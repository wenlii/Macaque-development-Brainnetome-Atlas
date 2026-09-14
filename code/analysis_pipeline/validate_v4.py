"""Output-based validation of v4; separate technical integrity from anatomy."""
import argparse,json,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import nibabel as nib
from atlas_pipeline import *
from refine_boundaries import graph,cc,source_parts,gm_votes


def readrows(p):return pd.read_csv(p,sep='\t').to_dict('records')


def specs(out,rep):
    folder,suffix=('native','') if rep=='native' else ('fsLR32k','.32k')
    root=ET.Element('CaretSpecFile',Version='1.0');ET.SubElement(root,'MetaData')
    entries=[('Invalid','VOLUME','anat/T1w.nii.gz'),('Invalid','VOLUME','volumes/MBNA248_T1w.nii.gz'),('Invalid','VOLUME','volumes/overlap_uncertainty_T1w.nii.gz')]
    for h in HEMIS:
        structure='CortexLeft' if h=='L' else 'CortexRight'
        for kind in ['white','pial','midthickness','inflated']:entries.append((structure,'SURFACE',f'{folder}/{h}.{kind}{suffix}.surf.gii'))
        entries.append((structure,'LABEL',f'labels/{h}.MBNA124.{rep}.label.gii'))
        for kind in ['cortex','geometry_flags','boundary_added','repair_flags','unassigned']:entries.append((structure,'METRIC',f'{folder}/{h}.{kind}{suffix}.shape.gii'))
    for st,typ,name in entries:
        assert (out/name).is_file(),name;ET.SubElement(root,'DataFile',Structure=st,DataFileType=typ,Selected='true').text=name
    ET.indent(root,space='  ');ET.ElementTree(root).write(out/f'atlas_{rep}.wb.spec',encoding='utf-8',xml_declaration=True)


def surfaces():
    checks=[];summaries=[];rois=[];parts_table=[]
    def ck(condition,name,group,h,rep):
        checks.append({'group':group,'hemisphere':h,'representation':rep,'check':name,'pass':bool(condition)});assert condition,(name,group,h,rep)
    for group in GROUPS:
        out=ROOT/'groups'/group
        for h,fs in HEMIS.items():
            for rep in ['native','32k']:
                folder,suffix=('native','') if rep=='native' else ('fsLR32k','.32k')
                base=values(out/f'baseline/{h}.MBNA124.{rep}.label.gii').astype(int)
                final=values(out/f'labels/{h}.MBNA124.{rep}.label.gii').astype(int)
                before=values(out/f'labels/{h}.{rep}_before_smoothing.label.gii').astype(int)
                mask=values(out/f'{folder}/{h}.cortex{suffix}.shape.gii')>0
                oldmask=values(out/f'baseline/{h}.cortex.{rep}.shape.gii')>0
                bad=values(out/(f'registration/{h}.full_geometry_flags.shape.gii' if rep=='native' else f'registration/{h}.geometry_32k.shape.gii'))>0
                v,f=mesh(out/f'{folder}/{h}.midthickness{suffix}.surf.gii');a,edges,length=graph(v,f);area=vertex_areas(v,f)
                ck(set(final[final>0])==set(range(1,125)),'124 ROIs',group,h,rep)
                ck(np.all(mask[oldmask]) and not np.any(final[~mask]),'v3 cortex retained and labels inside mask',group,h,rep)
                ck(np.array_equal(final[bad],base[bad]),'geometric caution labels not overwritten',group,h,rep)
                votes=gm_votes(out,h,rep);added=mask&~oldmask
                reference=values(out/f'registration/{h}.adult_labels_before_mask.label.gii').astype(int) if rep=='native' else values(ROOT/f'source/{h}.MBNA124.32k.label.gii').astype(int)
                ck(np.all(votes[added]>=2) and not np.any(bad[added]) and np.all(reference[added]>0),'new support has age GM and adult support with no geometry flag',group,h,rep)
                ck(np.array_equal(added,values(out/f'{folder}/{h}.boundary_added{suffix}.shape.gii')>0),'added support metric exact',group,h,rep)
                if rep=='32k':ck(np.all(values(out/f'registration/{h}.v4_cortex_fraction.32k.shape.gii')[added]>0),'32k additions have native support',group,h,rep)
                direct=values(out/f'registration/{h}.components.native.label.gii').astype(int) if rep=='native' else values(ROOT/f'source/{h}.components.32k.label.gii').astype(int)
                sphere,_=mesh(out/(f'registration/{h}.sphere.MSMSulc.surf.gii' if rep=='native' else f'fsLR32k/{h}.sphere.32k.surf.gii'))
                bp=source_parts(base,direct,sphere,h);fp=source_parts(final,direct,sphere,h)
                ck(set(bp[bp>0]).issubset(set(fp[fp>0])),'all v3 source parts retained',group,h,rep)
                for part in np.unique(bp[bp>0]):parts_table.append({'group':group,'hemisphere':h,'representation':rep,'source_component':int(part),'v3_vertices':int((bp==part).sum()),'v4_vertices':int((fp==part).sum())})
                replay=before.copy();smooth=readrows(out/f'metrics/{h}.{rep}.smooth.tsv');rejections=readrows(out/f'metrics/{h}.{rep}.final_rejections.tsv')
                loss=np.zeros(125)
                for event in smooth:
                    i=int(event['vertex']);old=int(event['old_label']);new=int(event['new_label'])
                    ck(replay[i]==old and not bad[i] and old>0 and new>0,'smoothing edit chain and geometry rule',group,h,rep)
                    loss[old]+=area[i];replay[i]=new
                for event in rejections:
                    i=int(event['vertex']);ck(replay[i]==event['old_label'],'rejection replay',group,h,rep);replay[i]=event['new_label']
                ck(np.array_equal(replay,final),'final labels reproduce from recorded smoothing',group,h,rep)
                initial_area=np.bincount(before,weights=area,minlength=125)
                ck(np.all(loss<=.05*initial_area+1e-8),'smoothing removes at most 5 percent of each ROI area',group,h,rep)
                changes=readrows(out/f'metrics/{h}.{rep}.vertex_changes.tsv')
                ck(set(int(x['vertex']) for x in changes)==set(np.flatnonzero((base!=final)|(mask!=oldmask))),'all final vertex changes recorded',group,h,rep)
                for roi in range(1,125):
                    n0=len(cc(a,base,roi));n1=len(cc(a,final,roi));ck(0<n1<=n0,'no additional ROI component',group,h,rep)
                    rois.append({'group':group,'hemisphere':h,'representation':rep,'source_id':roi,'v3_vertices':int((base==roi).sum()),'v4_vertices':int((final==roi).sum()),
                        'v3_components':n0,'v4_components':n1,'v3_area_mm2':float(area[base==roi].sum()),'v4_area_mm2':float(area[final==roi].sum()),'geometry_flagged_vertices':int(((final==roi)&bad).sum())})
                ck(h!='R' or len(cc(a,final,98))==2,'right V1.rm retains its two source components',group,h,rep)
                for kind in ['white','pial','midthickness','inflated']:
                    p=Path(f'{folder}/{h}.{kind}{suffix}.surf.gii');ck(sha(out/p)==sha(BASELINE/'groups'/group/p),'unchanged geometry '+kind,group,h,rep)
                source_table=nib.load(ROOT/f'source/{h}.MBNA124.32k.label.gii').labeltable.to_xml()
                ck(nib.load(out/f'labels/{h}.MBNA124.{rep}.label.gii').labeltable.to_xml()==source_table,'adult label dictionary unchanged',group,h,rep)
                d=json.loads((out/f'metrics/{h}.{rep}.summary.json').read_text());d.update(d.pop('boundary_smoothing'));summaries.append(d)
    # Exact candidate highlighted by the user; do not silently fill known flagged vertices.
    out=ROOT/'groups/Group04_59p5-68p0';ids=np.array([2055,2106,2107,2157,2206,2207])
    lab=values(out/'labels/R.MBNA124.32k.label.gii');mask=values(out/'fsLR32k/R.cortex.32k.shape.gii')
    assert np.all(lab[ids]==0) and np.all(mask[ids]==0)
    dump(ROOT/'qc/Group04_flagged_hole.json',{'group':'Group04_59p5-68p0','hemisphere':'R','vertices_32k':ids.tolist(),'all_remain_excluded':True,'reason':'all six carry pre-existing geometry caution; no new anatomical correction was performed'})
    table(ROOT/'metrics/surface_summary.tsv',summaries);table(ROOT/'metrics/ROI_changes.tsv',rois);table(ROOT/'metrics/source_component_preservation.tsv',parts_table)
    table(ROOT/'provenance/surface_checks.tsv',checks);dump(ROOT/'provenance/surface_validation.json',{'status':'PASS','checks':len(checks),'represented_source_parts_lost':0,'unchanged_geometry':True})
    print('Surface verification PASS',len(checks),flush=True)


def volumes():
    records=[];checks=[]
    def ck(ok,name,g):checks.append({'group':g,'check':name,'pass':bool(ok)});assert ok,(g,name)
    for group in GROUPS:
        out=ROOT/'groups'/group
        def load(name):return np.asanyarray(nib.load(out/'volumes'/name).dataobj).astype(int)
        merged=load('MBNA248_T1w.nii.gz');l=load('L.MBNA124_T1w.nii.gz');r=load('R.MBNA124_T1w.nii.gz')
        rl=load('L.MBNA124_resolved_T1w.nii.gz');rr=load('R.MBNA124_resolved_T1w.nii.gz');u=load('overlap_uncertainty_T1w.nii.gz')>0
        conflict=(l>0)&(r>0);choice=load('overlap_chosen_hemisphere_T1w.nii.gz')
        ck(set(np.unique(merged[merged>0]))==set(range(1,249)),'248 volume labels',group)
        ck(not np.any((rl>0)&(rr>0)) and np.array_equal(merged,np.where(rl>0,rl,np.where(rr>0,rr+124,0))),'disjoint hemisphere products reproduce merged volume',group)
        natural=np.where(l>0,l,np.where(r>0,r+124,0));ck(np.array_equal(merged[~conflict],natural[~conflict]),'ownership alters only overlap voxels',group)
        ck(np.all(np.isin(choice[conflict],[1,2])) and not np.any(choice[~conflict]),'one owner per overlap voxel',group)
        rows=readrows(out/'metrics/hemisphere_adjudication.tsv');ck(len(rows)==int(conflict.sum()),'one decision record per conflict',group)
        for x in rows:
            key=tuple(int(x[z]) for z in ['i','j','k']);ck(merged[key]==x['new_label'] and bool(u[key])==x['uncertainty_flag'],'adjudication record matches output',group)
        raw=nib.load(out/'anat/T1w.nii.gz');vol=nib.load(out/'volumes/MBNA248_T1w.nii.gz')
        ck(vol.shape==raw.shape and np.array_equal(vol.affine,raw.affine),'original template grid preserved',group)
        base=np.asanyarray(nib.load(out/'baseline/MBNA248_T1w.nii.gz').dataobj);ids=np.argwhere(base!=merged);xyz=nib.affines.apply_affine(raw.affine,ids)
        table(out/'metrics/voxel_changes.tsv',[{'i':int(i[0]),'j':int(i[1]),'k':int(i[2]),'x_mm':x[0],'y_mm':x[1],'z_mm':x[2],'v3_label':int(base[tuple(i)]),'v4_label':int(merged[tuple(i)]),'uncertain_owner':bool(u[tuple(i)])} for i,x in zip(ids,xyz)])
        rec=json.loads((out/'metrics/volume_summary.json').read_text());ck(rec['changed_volume_voxels_vs_v3']==len(ids),'summary changed voxel count exact',group);records.append(rec)
        for rep in ['native','32k']:
            specs(out,rep);command(out/f'logs/read_spec_{rep}',['-file-information',out/f'atlas_{rep}.wb.spec'])
    table(ROOT/'metrics/volume_summary.tsv',records);table(ROOT/'provenance/volume_checks.tsv',checks)
    dump(ROOT/'provenance/volume_validation.json',{'status':'PASS','checks':len(checks),'all_six_volumes_have_248_ROIs':True});print('Volume verification PASS',len(checks),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['surfaces','volumes']);a=p.parse_args();{'surfaces':surfaces,'volumes':volumes}[a.stage]()
