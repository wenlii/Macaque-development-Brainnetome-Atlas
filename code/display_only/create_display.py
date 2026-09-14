"""COSMETIC DISPLAY ONLY. Never use these outputs for quantitative analysis.

Create a separate, self-contained surface presentation from the v4 analysis
package. Geometry-flag exclusions are deliberately overridden for appearance.
No analysis label, mask, volume, or source geometry is overwritten.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,copy,csv,gc,hashlib,json,shutil
import numpy as np
import nibabel as nib
from scipy import sparse
from scipy.sparse.csgraph import connected_components,dijkstra

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
BASE=PROJECT/'macaque_age_atlas_v4_minimal'
ARCHIVE=PROJECT/'result/age_native_MBNA_v4'
GROUPS=['Group01_34p0-42p5','Group02_42p5-51p0','Group03_51p0-59p5','Group04_59p5-68p0','Group05_68p0-76p5','Group06_76p5-inf']
CLOSE_MM=2.0
MASK_SMOOTH_STEPS=6
LABEL_SMOOTH_STEPS=8
SUBDIVISIONS=2


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''):h.update(block)
    return h.hexdigest()


def dump(p,d):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def table(p,rows):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['no_records'],delimiter='\t');w.writeheader();w.writerows(rows)


def mesh(p):
    im=nib.load(p);return im.agg_data('pointset').astype(np.float32),im.agg_data('triangle').astype(np.int32)


def vals(p):return nib.load(p).agg_data().ravel()


def write_gifti(p,arrays,h,name,labeltable=None):
    im=nib.gifti.GiftiImage(darrays=arrays)
    if labeltable is not None:im.labeltable=copy.deepcopy(labeltable)
    im.meta['AnatomicalStructurePrimary']='CortexLeft' if h=='L' else 'CortexRight'
    im.meta['IntendedUse']='DISPLAY_ONLY_NOT_FOR_ANALYSIS'
    im.meta['Description']='Cosmetically filled and smoothed display derivative; anatomical validation does not apply.'
    im.meta['Name']=name
    nib.save(im,p)


def labels(p,x,h,lut):
    write_gifti(p,[nib.gifti.GiftiDataArray(x.astype(np.int32),intent='NIFTI_INTENT_LABEL')],h,'DISPLAY ONLY: cosmetically modified ROI labels',lut)


def metric(p,x,h,name):
    write_gifti(p,[nib.gifti.GiftiDataArray(x.astype(np.float32),intent='NIFTI_INTENT_SHAPE')],h,'DISPLAY ONLY: '+name)


def surface(p,v,f,h):
    write_gifti(p,[nib.gifti.GiftiDataArray(v.astype(np.float32),intent='NIFTI_INTENT_POINTSET'),nib.gifti.GiftiDataArray(f.astype(np.int32),intent='NIFTI_INTENT_TRIANGLE')],h,'DISPLAY ONLY: interpolated inflated mesh')


def adjacency(v,f):
    e=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
    length=np.linalg.norm(v[e[:,0]]-v[e[:,1]],axis=1)
    a=sparse.csr_matrix((np.r_[length,length],(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(v),len(v)))
    area=.5*np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)
    va=np.bincount(f.ravel(),weights=np.repeat(area/3,3),minlength=len(v))
    return a,e,length,va


def fill_background_islands(mask,a,area):
    ids=np.flatnonzero(~mask);n,cl=connected_components(a[ids][:,ids],directed=False)
    weights=np.bincount(cl,weights=area[ids]);keep=int(weights.argmax())
    out=mask.copy();out[ids[cl!=keep]]=True
    return out,n-1,int(np.count_nonzero(out&~mask))


def transition(a,area,mask=None):
    c=a.tocoo();w=(np.exp(-.5*(c.data/1.0)**2)*area[c.col]).astype(np.float32)
    if mask is not None:
        use=mask[c.row]&mask[c.col];c=sparse.coo_matrix((w[use],(c.row[use],c.col[use])),shape=a.shape)
    else:c=sparse.coo_matrix((w,(c.row,c.col)),shape=a.shape)
    mat=c.tocsr();total=np.asarray(mat.sum(1)).ravel();total[total==0]=1
    return sparse.diags(1/total)@mat


def subdivide(v,f):
    # Linear subdivision retains original points and the same piecewise surface.
    raw=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);edges,inv=np.unique(np.sort(raw,axis=1),axis=0,return_inverse=True)
    n=len(v);m=len(f);ab=n+inv[:m];bc=n+inv[m:2*m];ca=n+inv[2*m:]
    newv=np.concatenate([v,.5*(v[edges[:,0]]+v[edges[:,1]])]).astype(np.float32)
    newf=np.concatenate([np.c_[f[:,0],ab,ca],np.c_[ab,f[:,1],bc],np.c_[ca,bc,f[:,2]],np.c_[ab,bc,ca]]).astype(np.int32)
    rows=np.r_[np.arange(n),np.arange(n,n+len(edges)),np.arange(n,n+len(edges))]
    cols=np.r_[np.arange(n),edges[:,0],edges[:,1]]
    data=np.r_[np.ones(n),np.full(len(edges)*2,.5)].astype(np.float32)
    B=sparse.csr_matrix((data,(rows,cols)),shape=(len(newv),n))
    return newv,newf,B


def initialize():
    for d in ['provenance','figures','qc','groups']:(ROOT/d).mkdir(exist_ok=True)
    manifest=BASE/'manifest_sha256.tsv'
    rows=list(csv.DictReader(manifest.open(encoding='utf-8'),delimiter='\t'))
    for row in rows:assert sha(BASE/row['path'])==row['sha256'],row['path']
    shutil.copyfile(manifest,ROOT/'provenance/protected_v4_manifest.tsv')
    shutil.copyfile(BASE/'volume_labels.tsv',ROOT/'label_dictionary.tsv')
    shutil.copyfile(BASE/'age_groups.tsv',ROOT/'age_groups.tsv')
    dump(ROOT/'provenance/processing_contract.json',{
        'authorization':'User explicitly requested a cosmetically attractive separate version for display only, not downstream computation.',
        'base':str(BASE),'base_manifest_sha256':sha(manifest),'base_payload_files':len(rows),'groups':GROUPS,
        'intended_use':'DISPLAY_ONLY_NOT_FOR_ANALYSIS','geometry_exclusion_override':'Allowed only in this display derivative, including Group04 flagged holes',
        'mask_closing_mm':CLOSE_MM,'mask_smoothing_steps':MASK_SMOOTH_STEPS,'label_probability_diffusion_steps':LABEL_SMOOTH_STEPS,
        'hole_rule':'fill every excluded component except the largest-area medial-wall component before and after cosmetic smoothing',
        'medial_wall_core':'keep original background vertices at least 3 mm from cortex excluded',
        'ROI_protection':'preserve an interior seed for every original ROI component so all 124 ROI identities remain visible',
        'label_ids_averaged':False,'dense_display_subdivisions':SUBDIVISIONS,'dense_display_information':'linear geometry and categorical score interpolation; no extra MRI information',
        'analysis_volumes_generated':False,'analysis_files_modified':False,'quota_reset_cards_used':False,
        'figure_contract':{'backend':'Python/PyVista','claim':'cosmetic continuity and clean viewing only','archetype':'surface image plates','palette':'original adult ROI colors and unlit neutral medial wall','exports':'PNG with DISPLAY ONLY labels','no_anatomical_validation_claim':True}})


def process(group):
    out=ROOT/'groups'/group;out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for h in ['L','R']:
        source=BASE/'groups'/group
        v,f=mesh(source/f'{h}.midthickness.32k.surf.gii');inflated,fi=mesh(source/f'{h}.inflated.32k.surf.gii');assert np.array_equal(f,fi)
        li=nib.load(source/f'{h}.MBNA124.32k.label.gii');base=li.agg_data().astype(int)
        a,edges,length,area=adjacency(v,f);base_mask=base>0
        mask,hole_components,hole_vertices=fill_background_islands(base_mask,a,area)
        dist=dijkstra(a,indices=np.flatnonzero(mask),min_only=True,limit=CLOSE_MM)
        dilated=dist<=CLOSE_MM;to_out=dijkstra(a,indices=np.flatnonzero(~dilated),min_only=True,limit=CLOSE_MM)
        mask |= to_out>CLOSE_MM
        # Cosmetic signed-distance smoothing rounds the mask boundary.
        inside=dijkstra(a,indices=np.flatnonzero(~mask),min_only=True,limit=4.)
        outside=dijkstra(a,indices=np.flatnonzero(mask),min_only=True,limit=4.)
        signed=np.where(mask,np.minimum(inside,4.),-np.minimum(outside,4.)).astype(np.float32)
        T=transition(a,area)
        for _ in range(MASK_SMOOTH_STEPS):signed=.5*signed+.5*(T@signed)
        # Protect one visible seed for each existing ROI component, including tiny parts.
        bnodes=np.unique(edges[(base[edges[:,0]]!=base[edges[:,1]])])
        bd=dijkstra(a,indices=bnodes,min_only=True)
        seeds=[]
        for roi in range(1,125):
            ids=np.flatnonzero(base==roi);n,cl=connected_components(a[ids][:,ids],directed=False)
            for j in range(n):
                component=ids[cl==j];seed=int(component[np.argmax(bd[component])]);seeds.append(seed)
        seeds=np.array(seeds,int)
        mask=signed>0;mask[seeds]=True
        old_distance=dijkstra(a,indices=np.flatnonzero(base_mask),min_only=True,limit=3.)
        core=(~base_mask)&(old_distance>=3.);mask[core]=False
        mask,nmore,nverts=fill_background_islands(mask,a,area)
        signed=(np.abs(signed)+.05)*np.where(mask,1.,-1.)
        # Initial cosmetic fill prefers the common adult reference where available.
        donors=np.flatnonzero(base_mask)
        distance,pred,owner=dijkstra(a,indices=donors,min_only=True,return_predecessors=True)
        filled=base.copy();fill=(base==0)&mask
        filled[fill]=base[owner[fill]]
        adult=vals(ARCHIVE/f'source/{h}.MBNA124.32k.label.gii').astype(int)
        good=fill&(adult>0);filled[good]=adult[good];filled[~mask]=0
        score=np.zeros((len(v),124),np.float32);ids=np.flatnonzero(mask);score[ids,filled[ids]-1]=1.
        tm=transition(a,area,mask)
        for _ in range(LABEL_SMOOTH_STEPS):
            score=.5*score+.5*(tm@score)
            score[seeds]=0.;score[seeds,base[seeds]-1]=1.
        display=np.where(mask,score.argmax(1)+1,0).astype(np.int32)
        assert set(display[display>0])==set(range(1,125))
        z=np.flatnonzero(display==0);assert connected_components(a[z][:,z],directed=False,return_labels=False)==1
        labels(out/f'{h}.MBNA124.DISPLAY_ONLY.32k.label.gii',display,h,li.labeltable)
        surface(out/f'{h}.inflated.DISPLAY_ONLY.32k.surf.gii',inflated,f,h)
        metric(out/f'{h}.cortex.DISPLAY_ONLY.32k.shape.gii',mask,h,'cosmetic display coverage')
        geom=vals(ARCHIVE/'groups'/group/f'registration/{h}.geometry_32k.shape.gii')>0
        flags=np.zeros(len(v),np.uint8);flags[(base==0)&(display>0)]|=1;flags[(base>0)&(display!=base)]|=2;flags[(base==0)&(display>0)&geom]|=4
        metric(out/f'{h}.cosmetic_changes.DISPLAY_ONLY.32k.shape.gii',flags,h,'1=cosmetic fill;2=other label change;4=fill over prior geometry caution')
        changed=np.flatnonzero(display!=base)
        table(out/f'{h}.cosmetic_vertex_changes.tsv',[{'vertex':int(i),'v4_analysis_label':int(base[i]),'display_label':int(display[i]),'cosmetic_flags':int(flags[i])} for i in changed])
        # A dense presentation mesh reduces stair steps without changing MRI resolution.
        vd,fd=inflated,f;B=sparse.eye(len(v),dtype=np.float32,format='csr')
        for _ in range(SUBDIVISIONS):vd,fd,b=subdivide(vd,fd);B=b@B
        md=np.asarray(B@signed).ravel()>0
        dense=np.zeros(len(vd),np.int32)
        for start in range(0,len(vd),40000):
            stop=min(start+40000,len(vd));pr=B[start:stop]@score
            dense[start:stop]=np.where(md[start:stop],pr.argmax(1)+1,0)
        assert set(dense[dense>0])==set(range(1,125))
        assert np.array_equal(dense[:len(v)],display)
        labels(out/f'{h}.MBNA124.DISPLAY_ONLY.dense.label.gii',dense,h,li.labeltable)
        surface(out/f'{h}.inflated.DISPLAY_ONLY.dense.surf.gii',vd,fd,h)
        row={'group':group,'hemisphere':h,'source_vertices':len(v),'dense_display_vertices':len(vd),'dense_faces':len(fd),'ROIs_retained':124,
            'v4_zero_vertices':int((base==0).sum()),'display_zero_vertices_32k':int((display==0).sum()),'display_background_components_32k':1,
            'cosmetically_filled_vertices':int(((base==0)&(display>0)).sum()),'other_changed_vertices':int(((base>0)&(display!=base)).sum()),
            'filled_over_prior_geometry_caution':int(((flags&4)>0).sum()),'initial_enclosed_holes_filled':hole_components,'initial_hole_vertices_filled':hole_vertices,
            'intended_use':'DISPLAY_ONLY_NOT_FOR_ANALYSIS'}
        rows.append(row);print(row,flush=True)
        del B,vd,fd,score,a,tm,T;gc.collect()
    dump(out/'display_only.json',{'intended_use':'DISPLAY_ONLY_NOT_FOR_ANALYSIS','group':group,'analysis_package':'macaque_age_atlas_v4_minimal','summary':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['init','group']);p.add_argument('--group',choices=GROUPS);args=p.parse_args()
    initialize() if args.stage=='init' else process(args.group)
