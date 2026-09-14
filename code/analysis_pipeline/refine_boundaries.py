"""Anatomy-constrained mask closing and topology-preserving categorical smoothing.

The v3 baseline and all mesh coordinates remain protected. Physical distances
use midthickness; no arithmetic smoothing of integer label IDs is performed.
"""
import argparse
import numpy as np
import pandas as pd
import nibabel as nib
from scipy import sparse
from scipy.sparse.csgraph import connected_components,dijkstra
from scipy.spatial import cKDTree
from scipy.ndimage import map_coordinates
from atlas_pipeline import *

RADIUS=1.0
HOLE_AREA=2.0
SIGMA=.6
SMOOTH_PASSES=2
MAX_AREA_LOSS=.05


def graph(v,f):
    edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    length=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1)
    a=sparse.csr_matrix((np.r_[length,length],(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(len(v),len(v)))
    return a,edges,length


def restrict(a,mask):
    c=a.tocoo();ok=mask[c.row]&mask[c.col]
    return sparse.csr_matrix((c.data[ok],(c.row[ok],c.col[ok])),shape=a.shape)


def cc(a,labels,key):
    ids=np.flatnonzero(labels==key)
    if not len(ids):return []
    n,parts=connected_components(a[ids][:,ids],directed=False)
    return sorted([ids[parts==i] for i in range(n)],key=lambda z:(-len(z),int(z.min())))


def neighbors(a,i):return a.indices[a.indptr[i]:a.indptr[i+1]]


def source_parts(labels,direct,sphere,h):
    src=values(ROOT/f'source/{h}.MBNA124.32k.label.gii').astype(int)
    parts=values(ROOT/f'source/{h}.components.32k.label.gii').astype(int)
    sp,_=mesh(ROOT/f'source/{h}.sphere.32k.surf.gii')
    info=pd.read_csv(ROOT/'source/components.tsv',sep='\t').query('hemisphere==@h')
    lut=np.zeros(int(info.component_id.max())+1,int);lut[info.component_id]=info.source_id
    result=np.zeros(len(labels),int);ok=(labels>0)&(direct>0)&(direct<len(lut))
    ok &= lut[np.clip(direct,0,len(lut)-1)]==labels;result[ok]=direct[ok]
    for roi in range(1,125):
        need=np.flatnonzero((labels==roi)&(result==0));ids=np.flatnonzero(src==roi)
        if len(need):result[need]=parts[ids[cKDTree(sp[ids]).query(sphere[need])[1]]]
    return result


def gm_votes(out,h,rep):
    suffix='' if rep=='native' else '.32k';folder='native' if rep=='native' else 'fsLR32k'
    w,f=mesh(out/f'{folder}/{h}.white{suffix}.surf.gii');p,pf=mesh(out/f'{folder}/{h}.pial{suffix}.surf.gii')
    assert np.array_equal(f,pf)
    im=nib.load(out/'anat/age_aseg.nii.gz');data=np.asanyarray(im.dataobj);inv=np.linalg.inv(im.affine)
    samples=np.array([map_coordinates(data,nib.affines.apply_affine(inv,(1-t)*w+t*p).T,order=0,mode='constant',cval=0) for t in [.25,.5,.75]])
    return (samples==(3 if h=='L' else 42)).sum(0)


def close_candidates(a,cortex,area):
    distance=dijkstra(a,indices=np.flatnonzero(cortex),min_only=True,limit=RADIUS)
    dilated=distance<=RADIUS
    eroded_distance=dijkstra(a,indices=np.flatnonzero(~dilated),min_only=True,limit=RADIUS)
    closing=(eroded_distance>RADIUS)&~cortex&(distance<=RADIUS)
    small_holes=np.zeros(len(cortex),bool)
    background=cc(a,(~cortex).astype(int),1)
    if background:
        largest=max(background,key=lambda x:area[x].sum())
        for ids in background:
            if ids is largest:continue
            if area[ids].sum()<=HOLE_AREA:small_holes[ids]=True
    return closing|small_holes,closing,small_holes,distance


def topology_restore(labels,mask,base,base_mask,a,sphere,direct,h):
    """Reject resampling-created extras by restoring their exact v3 state."""
    result=labels.copy();mask=mask.copy();events=[]
    baseparts=source_parts(base,direct,sphere,h)
    for iteration in range(12):
        changed=False;parts=source_parts(result,direct,sphere,h)
        for roi in range(1,125):
            comps=cc(a,result,roi)
            present=np.unique(baseparts[base==roi]);present=present[present>0]
            for part in present:
                if not np.any(parts==part):
                    ids=np.flatnonzero(baseparts==part)
                    for i in ids:
                        if result[i]!=base[i]:events.append({'vertex':int(i),'old_label':int(result[i]),'new_label':int(base[i]),'method':'restore_v3_source_part'})
                    result[ids]=base[ids];mask[ids]=base_mask[ids];changed=True
            comps=cc(a,result,roi)
            if len(comps)<=len(cc(a,base,roi)):continue
            keep=set()
            for part in present:
                scores=[(int(np.count_nonzero(baseparts[ids]==part)),len(ids)) for ids in comps]
                if scores:keep.add(max(range(len(scores)),key=lambda j:scores[j]))
            for j,ids in enumerate(comps):
                if j in keep:continue
                altered=ids[result[ids]!=base[ids]]
                for i in altered:events.append({'vertex':int(i),'old_label':int(result[i]),'new_label':int(base[i]),'method':'reject_resampling_extra_restore_v3'})
                result[altered]=base[altered];mask[altered]=base_mask[altered]
                changed |= len(altered)>0
        if not changed:break
    for roi in range(1,125):
        if len(cc(a,result,roi))>len(cc(a,base,roi)):
            # Restore all changes involving a problematic ROI, then verify globally.
            ids=np.flatnonzero((result!=base)&((result==roi)|(base==roi)))
            for i in ids:events.append({'vertex':int(i),'old_label':int(result[i]),'new_label':int(base[i]),'method':'restore_v3_problematic_ROI_neighborhood'})
            result[ids]=base[ids];mask[ids]=base_mask[ids]
    return result,mask,events


def add_supported(labels,mask,candidates,reference,direct,sphere,h,a,bad):
    result=labels.copy();mask=mask.copy();events=[]
    parts=source_parts(result,direct,sphere,h);refparts=source_parts(reference,direct,sphere,h)
    pending=set(map(int,np.flatnonzero(candidates&(result==0)&(reference>0)&~bad)))
    for iteration in range(15):
        n=0
        for i in sorted(pending):
            roi=int(reference[i]);ns=neighbors(a,i);donors=ns[(result[ns]==roi)&~bad[ns]]
            if not len(donors):continue
            observed=np.unique(parts[ns[result[ns]==roi]]);observed=observed[observed>0]
            if len(observed)!=1 or observed[0]!=refparts[i]:continue
            events.append({'vertex':i,'old_label':0,'new_label':roi,'method':'anatomy_supported_cortex_or_hole_fill'})
            result[i]=roi;mask[i]=True;parts[i]=refparts[i];pending.remove(i);n+=1
        if not n:break
    return result,mask,events


def smooth_labels(labels,mask,bad,a,area,reference,direct,sphere,h,edges,length):
    result=labels.copy();base=labels.copy();events=[];parts=source_parts(base,direct,sphere,h)
    protected=np.zeros(len(labels),bool)
    initial_area=np.bincount(base,weights=area,minlength=125);loss=np.zeros(125)
    for roi in range(1,125):
        if np.count_nonzero(base==roi)<=30:protected[base==roi]=True
        for ids in cc(a,base,roi):
            protected[ids[0]]=True
        for part in np.unique(parts[base==roi]):
            ids=np.flatnonzero(parts==part)
            if len(ids)<=20:protected[ids]=True
            elif len(ids):protected[ids[0]]=True
    original_boundary=np.zeros(len(labels),bool)
    use=(base[edges[:,0]]>0)&(base[edges[:,1]]>0)&(base[edges[:,0]]!=base[edges[:,1]])
    original_boundary[np.unique(edges[use])]=True
    def boundary_length(x):return float(length[(x[edges[:,0]]>0)&(x[edges[:,1]]>0)&(x[edges[:,0]]!=x[edges[:,1]])].sum())
    before_length=boundary_length(base)
    for iteration in range(SMOOTH_PASSES):
        proposals=[]
        for i in np.flatnonzero(original_boundary&mask&~bad&~protected):
            old=int(result[i]);ns=neighbors(a,i);dist=a.data[a.indptr[i]:a.indptr[i+1]]
            good=(result[ns]>0)&mask[ns]&~bad[ns]
            ns=ns[good];dist=dist[good]
            if not len(ns):continue
            weights=area[ns]*np.exp(-.5*(dist/SIGMA)**2)
            score=np.bincount(result[ns],weights=weights,minlength=125);score[old]+=area[i]
            target=int(score.argmax());total=float(score.sum())
            if target==old or target==0 or score[target]<.55*total or score[target]-score[old]<.15*total:continue
            if target not in reference[np.r_[i,ns]]:continue
            proposals.append((float(score[target]-score[old]),int(i),old,target))
        changed=0
        for advantage,i,old,target in sorted(proposals,key=lambda x:(-x[0],x[1])):
            if result[i]!=old or loss[old]+area[i]>.05*initial_area[old]:continue
            ns=neighbors(a,i);dist=a.data[a.indptr[i]:a.indptr[i+1]]
            oldns=ns[result[ns]==old];newns=ns[result[ns]==target]
            if not len(oldns) or not len(newns):continue
            if connected_components(a[oldns][:,oldns],directed=False,return_labels=False)!=1:continue
            if connected_components(a[newns][:,newns],directed=False,return_labels=False)!=1:continue
            target_parts=np.unique(parts[newns]);target_parts=target_parts[target_parts>0]
            if len(target_parts)!=1:continue
            delta=float(dist[result[ns]==old].sum()-dist[result[ns]==target].sum())
            if delta>=-1e-8:continue
            result[i]=target;parts[i]=int(target_parts[0]);loss[old]+=area[i];changed+=1
            events.append({'vertex':i,'old_label':old,'new_label':target,'method':'topology_preserving_boundary_vote','pass':iteration,'boundary_length_change_mm':delta})
        if not changed:break
    assert boundary_length(result)<=before_length+1e-7
    return result,events,{'before_mm':before_length,'after_mm':boundary_length(result),'sigma_mm':SIGMA,'passes':SMOOTH_PASSES,'max_loss_fraction_per_ROI':MAX_AREA_LOSS}


def process_rep(group,h,rep):
    out=ROOT/'groups'/group;is_native=rep=='native';folder='native' if is_native else 'fsLR32k';suffix='' if is_native else '.32k'
    v,f=mesh(out/f'{folder}/{h}.midthickness{suffix}.surf.gii');a,edges,length=graph(v,f);area=vertex_areas(v,f)
    base=values(out/f'baseline/{h}.MBNA124.{rep}.label.gii').astype(int);basemask=values(out/f'baseline/{h}.cortex.{rep}.shape.gii')>0
    sphere,_=mesh(out/(f'registration/{h}.sphere.MSMSulc.surf.gii' if is_native else f'fsLR32k/{h}.sphere.32k.surf.gii'))
    votes=gm_votes(out,h,rep)
    if is_native:
        direct=values(out/f'registration/{h}.components.native.label.gii').astype(int)
        reference=values(out/f'registration/{h}.adult_labels_before_mask.label.gii').astype(int)
        bad=values(out/f'registration/{h}.full_geometry_flags.shape.gii')>0
        start=base.copy();mask=basemask.copy();resampling_events=[];fraction=np.ones(len(v))
    else:
        direct=values(ROOT/f'source/{h}.components.32k.label.gii').astype(int)
        reference=values(ROOT/f'source/{h}.MBNA124.32k.label.gii').astype(int)
        bad=values(out/f'registration/{h}.geometry_32k.shape.gii')>0
        fraction=values(out/f'registration/{h}.v4_cortex_fraction.32k.shape.gii')
        raw=values(out/f'registration/{h}.v4_native_to_32k.label.gii').astype(int)
        # Restrict all new support, including resampling support, to current anatomy.
        mask=basemask|((fraction>=.5)&(reference>0)&(votes>=2)&~bad)
        start=np.where(mask,raw,0);start[bad]=base[bad]
        start,mask,resampling_events=topology_restore(start,mask,base,basemask,a,sphere,direct,h)
    proposal,closing,holes,distance=close_candidates(a,mask,area)
    eligible=proposal&(reference>0)&(votes>=2)&~bad
    if not is_native:eligible &= fraction>0
    # Remaining in-mask zeros need the same anatomy/source support and a bounded path.
    donors=np.flatnonzero(mask&(start>0)&~bad)
    dist_to_label=dijkstra(restrict(a,(mask|eligible)&~bad),indices=donors,min_only=True,limit=RADIUS)
    holes_inside=mask&(start==0)&~bad&(votes>=2)&(reference>0)&np.isfinite(dist_to_label)
    eligible |= holes_inside
    filled,newmask,fill_events=add_supported(start,mask,eligible,reference,direct,sphere,h,a,bad)
    final,smooth_events,boundary=smooth_labels(filled,newmask,bad,a,area,reference,direct,sphere,h,edges,length)
    final,newmask,reject_events=topology_restore(final,newmask,base,basemask,a,sphere,direct,h)
    for roi in range(1,125):
        assert 0<len(cc(a,final,roi))<=len(cc(a,base,roi)),(group,h,rep,roi)
    oldparts=source_parts(base,direct,sphere,h);parts=source_parts(final,direct,sphere,h)
    assert set(oldparts[oldparts>0]).issubset(set(parts[parts>0])),(group,h,rep,'source part lost')
    assert np.array_equal(final[bad],base[bad]),'Flagged labels changed'
    assert not np.any(final[~newmask]) and np.all(newmask[basemask]),'v3 support was removed'
    prefix=f'{h}.{rep}'
    write_labels(out/f'labels/{h}.MBNA124.{rep}.label.gii',final,h,ROOT/f'source/{h}.MBNA124.32k.label.gii')
    write_labels(out/f'labels/{h}.{rep}_before_smoothing.label.gii',filled,h,ROOT/f'source/{h}.MBNA124.32k.label.gii')
    write_metric(out/f'{folder}/{h}.cortex{suffix}.shape.gii',newmask,h,'v4 anatomy-supported cortex')
    write_metric(out/f'{folder}/{h}.unassigned{suffix}.shape.gii',newmask&(final==0),h,'v4 remaining unassigned')
    geom=values(out/f'registration/{h}.full_geometry_flags.shape.gii') if is_native else bad
    write_metric(out/f'{folder}/{h}.geometry_flags{suffix}.shape.gii',np.where(newmask,geom,0),h,'Unchanged geometry caution in final support')
    write_metric(out/f'{folder}/{h}.boundary_added{suffix}.shape.gii',newmask&~basemask,h,'v4 new support versus v3')
    flags=np.zeros(len(v),np.uint8);flags[(final!=base)&(base==0)]|=1;flags[(final!=base)&(base>0)]|=2;flags[newmask&~basemask]|=4
    write_metric(out/f'{folder}/{h}.repair_flags{suffix}.shape.gii',flags,h,'v4:1=filled label,2=changed label,4=expanded cortex')
    write_metric(out/f'registration/{h}.v4_source_parts.{rep}.shape.gii',parts,h,'Adult source component identity')
    write_metric(out/f'registration/{h}.v4_GM_votes.{rep}.shape.gii',votes,h,'Ipsilateral GM samples out of 3')
    pids=np.flatnonzero(proposal|holes_inside)
    table(out/f'metrics/{prefix}.mask_candidates.tsv',[{'vertex':int(i),'source_id':int(reference[i]),'closing':bool(closing[i]),'closed_hole':bool(holes[i]),
        'GM_votes':int(votes[i]),'geometry_flag':bool(bad[i]),'native_cortex_fraction':float(fraction[i]),'eligible':bool(eligible[i]),
        'accepted_new_support':bool(newmask[i] and not basemask[i]),'final_label':int(final[i]),'x_mm':v[i,0],'y_mm':v[i,1],'z_mm':v[i,2]} for i in pids])
    ids=np.flatnonzero((final!=base)|(newmask!=basemask))
    table(out/f'metrics/{prefix}.vertex_changes.tsv',[{'vertex':int(i),'old_label':int(base[i]),'new_label':int(final[i]),'old_cortex':bool(basemask[i]),'new_cortex':bool(newmask[i]),
        'repair_flags':int(flags[i]),'GM_votes':int(votes[i]),'geometry_flag':bool(bad[i]),'x_mm':v[i,0],'y_mm':v[i,1],'z_mm':v[i,2]} for i in ids])
    for name,data in [('fill',fill_events),('smooth',smooth_events),('resample_rejections',resampling_events),('final_rejections',reject_events)]:
        # Records within a file have one schema; absent operations retain an explicit empty table.
        table(out/f'metrics/{prefix}.{name}.tsv',data)
    roirows=[]
    for roi in range(1,125):
        roirows.append({'group':group,'hemisphere':h,'representation':rep,'source_id':roi,'v3_vertices':int((base==roi).sum()),'v4_vertices':int((final==roi).sum()),
            'v3_area_mm2':float(area[base==roi].sum()),'v4_area_mm2':float(area[final==roi].sum()),
            'v3_components':len(cc(a,base,roi)),'v4_components':len(cc(a,final,roi)),'geometry_flagged_vertices':int(((final==roi)&bad).sum())})
    table(out/f'metrics/{prefix}.ROIs.tsv',roirows)
    summary={'group':group,'hemisphere':h,'representation':rep,'ROIs':124,'changed_vertices':len(ids),'new_cortex_vertices':int((newmask&~basemask).sum()),
        'newly_assigned_vertices':int(((base==0)&(final>0)).sum()),'changed_nonzero_labels':int(((base>0)&(final!=base)).sum()),
        'unassigned_vertices':int((newmask&(final==0)).sum()),'geometry_flagged_candidates_retained':int((proposal&bad).sum()),'CGRSr_vertices':int((final==122).sum()),
        'source_parts_lost':0,'extra_components':0,'boundary_smoothing':boundary}
    dump(out/f'metrics/{prefix}.summary.json',summary);print(summary,flush=True)


def group_stage(group):
    out=ROOT/'groups'/group
    for h in HEMIS:
        process_rep(group,h,'native')
        command(out/f'logs/{h}_v4_native_to_32k',['-label-resample',out/f'labels/{h}.MBNA124.native.label.gii',out/f'registration/{h}.sphere.MSMSulc.surf.gii',out/f'fsLR32k/{h}.sphere.32k.surf.gii',
            'ADAP_BARY_AREA',out/f'registration/{h}.v4_native_to_32k.label.gii','-area-surfs',out/f'native/{h}.midthickness.surf.gii',out/f'fsLR32k/{h}.midthickness.32k.surf.gii','-current-roi',out/f'native/{h}.cortex.shape.gii'])
        command(out/f'logs/{h}_v4_cortex_fraction',['-metric-resample',out/f'native/{h}.cortex.shape.gii',out/f'registration/{h}.sphere.MSMSulc.surf.gii',out/f'fsLR32k/{h}.sphere.32k.surf.gii',
            'BARYCENTRIC',out/f'registration/{h}.v4_cortex_fraction.32k.shape.gii'])
        process_rep(group,h,'32k')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--group',choices=GROUPS,required=True);a=p.parse_args();group_stage(a.group)
