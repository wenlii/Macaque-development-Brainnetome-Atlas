"""Rasterize v4 labels and adjudicate bilateral conflicts in age-template space."""
import os
os.environ.setdefault("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS","4")
os.environ.setdefault("OMP_NUM_THREADS","4")
import argparse,json,shutil
import numpy as np
import pandas as pd
import nibabel as nib
import ants
from scipy.spatial import cKDTree
from atlas_pipeline import *

ADULT_T1 = PROJECT / "result/step1_v1/bridge/space-atlas_desc-affineResampled_T1w.nii.gz"


def hemisphere_prior(group):
    out=ROOT/'groups'/group
    priorfile=out/'registration/adult_signed_midline_distance.crop.nii.gz'
    if priorfile.exists() and (out/'metrics/hemisphere_prior_QC.json').exists(): return
    freeze(ADULT_T1)
    if group.startswith('Group03'):
        chain_path=PROJECT/'result/step2_v1/transforms/syn_chain.json'
        chain=json.loads(chain_path.read_text());paths=chain['fwd_image_transforms_moving_to_fixed']
    else:
        chain_path=PROJECT/'result/step5_v1/groups'/group/'transforms/chain.json'
        chain=json.loads(chain_path.read_text());paths=chain['forward_image']
        prep_path=PROJECT/'result/step5_v1/groups'/group/'prepared/preparation.json';freeze(prep_path)
        prep=json.loads(prep_path.read_text())
        assert sha(ADULT_T1)==prep['moving']['sha256']
        assert sha(out/'anat/T1w.nii.gz')==prep['fixed']['sha256']
    freeze(chain_path)
    xfms=out/'registration/volume_xfms';xfms.mkdir(exist_ok=True)
    localpaths=[]
    for path in paths:
        freeze(path);dest=xfms/Path(path).name;shutil.copyfile(path,dest);localpaths.append(str(dest))
    dump(xfms/'chain.json',{'source_chain':str(chain_path),'forward_image':localpaths,'flags':[False]*len(localpaths),'moving_reference':str(ADULT_T1)})
    left,_=mesh(ROOT/'source/L.midthickness.32k.surf.gii');right,_=mesh(ROOT/'source/R.midthickness.32k.surf.gii')
    mirror=right.copy();mirror[:,0]*=-1
    assert np.max(np.linalg.norm(left-mirror,axis=1))<1e-5
    sourcefile=ROOT/'source/adult_signed_midline_distance.nii.gz'
    if not sourcefile.exists():
        ref=nib.load(ADULT_T1);a=ref.affine;n=ref.shape
        x=(a[0,0]*np.arange(n[0],dtype=np.float32)[:,None,None]+a[0,1]*np.arange(n[1],dtype=np.float32)[None,:,None]
           +a[0,2]*np.arange(n[2],dtype=np.float32)[None,None,:]+a[0,3])
        write_volume(sourcefile,x,ref,label=False)
    moving=ants.image_read(str(sourcefile));fixed=ants.image_read(str(out/'anat/T1w_reconstruction_crop.nii.gz'))
    warped=ants.apply_transforms(fixed,moving,localpaths,whichtoinvert=[False]*len(localpaths),interpolator='linear')
    coverage=ants.apply_transforms(fixed,moving.new_image_like(np.ones(moving.shape,np.float32)),localpaths,
        whichtoinvert=[False]*len(localpaths),interpolator='linear')
    ants.image_write(warped,str(priorfile));ants.image_write(coverage,str(out/'registration/adult_prior_coverage.crop.nii.gz'))
    prior=warped.numpy();valid=coverage.numpy()>.99
    asegim=nib.load(out/'anat/age_aseg.nii.gz');aseg=np.asanyarray(asegim.dataobj)
    age=np.where(np.isin(aseg,[2,3]),1,np.where(np.isin(aseg,[41,42]),2,0))
    hemisphere=np.where(valid,np.where(prior>0,2,1),0).astype(np.int16)
    write_volume(out/'registration/adult_hemisphere_prior.crop.nii.gz',hemisphere,asegim)
    cortex=np.isin(aseg,[3,42]) & valid
    agreement=float(np.mean(hemisphere[cortex]==age[cortex]))
    # Independent API path checks image-versus-point transform conventions.
    indices=np.argwhere(cortex);rng=np.random.default_rng(20260912)
    sample=indices[rng.choice(len(indices),min(500,len(indices)),replace=False)]
    xyz=nib.affines.apply_affine(asegim.affine,sample);lps=xyz.copy();lps[:,:2]*=-1
    transformed=ants.apply_transforms_to_points(3,pd.DataFrame(lps,columns=['x','y','z']),localpaths,
                                               whichtoinvert=[False]*len(localpaths))
    point_distance=-transformed['x'].values
    error=np.abs(point_distance-prior[tuple(sample.T)])
    assert np.percentile(error,95)<.1,'Image/point direction mismatch'
    assert agreement>.95,'Adult hemisphere prior disagrees with age anatomy; inspect before use'
    dump(out/'metrics/hemisphere_prior_QC.json',{'source_plane':'adult atlas RAS x=0, confirmed by exact mirrored left/right geometry',
        'age_cortical_aseg_agreement':agreement,'image_point_difference_p95_mm':float(np.percentile(error,95)),
        'image_point_difference_max_mm':float(error.max()),'cortical_coverage_fraction':float(valid[np.isin(aseg,[3,42])].mean()),
        'confidence_margin_mm':.25,'geometry_reference':'existing adult-to-age volumetric registration; no new registration fit',
        'interpretation':'hemisphere prior, not independent cortical parcellation truth'})
    print(group,'hemisphere prior agreement',agreement,flush=True)


def volumes(group):
    out=ROOT/'groups'/group
    hemisphere_prior(group)
    crop=nib.load(out/'anat/T1w_reconstruction_crop.nii.gz');raw=nib.load(out/'anat/T1w.nii.gz')
    offset=np.rint((np.linalg.inv(raw.affine)@crop.affine)[:3,3]).astype(int)
    assert np.allclose((np.linalg.inv(raw.affine)@crop.affine)[:3,:3],np.eye(3),atol=1e-5)
    slices=tuple(slice(int(x),int(x+n)) for x,n in zip(offset,crop.shape))
    arrays={};trees={}
    for h in HEMIS:
        dest=out/f'volumes/{h}.MBNA124.reconstruction_crop.nii.gz'
        command(out/f'logs/{h}_label_to_volume',['-label-to-volume-mapping',out/f'labels/{h}.MBNA124.native.label.gii',
            out/f'native/{h}.midthickness.surf.gii',out/'anat/T1w_reconstruction_crop.nii.gz',dest,'-ribbon-constrained',
            out/f'native/{h}.white.surf.gii',out/f'native/{h}.pial.surf.gii','-voxel-subdiv','3'])
        arrays[h]=np.asanyarray(nib.load(dest).dataobj).astype(np.int16)
        full=np.zeros(raw.shape,np.int16);full[slices]=arrays[h];write_volume(out/f'volumes/{h}.MBNA124_T1w.nii.gz',full,raw)
        v,_=mesh(out/f'native/{h}.midthickness.surf.gii')
        lab=values(out/f'labels/{h}.MBNA124.native.label.gii')
        bad=values(out/f'registration/{h}.full_geometry_flags.shape.gii')>0
        trees[h]=cKDTree(v[(lab>0)&~bad])
    l,r=arrays['L'],arrays['R'];conflict=(l>0)&(r>0)
    merged=np.where(l>0,l,np.where(r>0,r+124,0)).astype(np.int16)
    prior=np.asanyarray(nib.load(out/'registration/adult_signed_midline_distance.crop.nii.gz').dataobj)
    coverage=np.asanyarray(nib.load(out/'registration/adult_prior_coverage.crop.nii.gz').dataobj)>.99
    aseg=np.asanyarray(nib.load(out/'anat/age_aseg.nii.gz').dataobj)
    coords=np.argwhere(conflict);world=nib.affines.apply_affine(crop.affine,coords)
    dl=trees['L'].query(world)[0] if len(world) else np.array([])
    dr=trees['R'].query(world)[0] if len(world) else np.array([])
    records=[];uncertain=np.zeros(crop.shape,bool);chosen=np.zeros(crop.shape,np.int16)
    for j,(idx,xyz) in enumerate(zip(coords,world)):
        key=tuple(idx);signed=float(prior[key]);covered=bool(coverage[key]);adult=2 if signed>0 else 1
        tissue=int(aseg[key]);age=1 if tissue in [2,3] else 2 if tissue in [41,42] else 0
        nearest=1 if dl[j]<dr[j] else 2;near=abs(signed)<.25
        if covered and not near and age in [0,adult]:
            side=adult;rule='adult_prior_supported';flag=age==0
        elif age and (not covered or near or nearest==age):
            side=age;rule='age_aseg_supported';flag=not covered or age!=adult
        elif covered and not near:
            side=adult;rule='adult_prior_with_age_disagreement';flag=True
        else:
            side=nearest;rule='nearest_reliable_hemisphere_surface';flag=True
            if abs(dl[j]-dr[j])<.1 and covered:side=adult;rule='adult_prior_near_surface_tie'
        selected=int(l[key]) if side==1 else int(r[key])+124
        merged[key]=selected;chosen[key]=side;uncertain[key]=flag
        raw_idx=idx+offset
        records.append({'i':int(raw_idx[0]),'j':int(raw_idx[1]),'k':int(raw_idx[2]),'x_mm':xyz[0],'y_mm':xyz[1],'z_mm':xyz[2],
            'left_candidate':int(l[key]),'right_candidate':int(r[key])+124,'adult_signed_distance_mm':signed,'adult_prior_covered':covered,
            'age_aseg_label':tissue,'age_hemisphere':age,'left_surface_distance_mm':float(dl[j]),'right_surface_distance_mm':float(dr[j]),
            'chosen_hemisphere':'L' if side==1 else 'R','new_label':selected,'rule':rule,'uncertainty_flag':flag})
    for name,array in [('MBNA248_T1w.nii.gz',merged),('overlap_before_resolution_T1w.nii.gz',conflict),
                       ('overlap_chosen_hemisphere_T1w.nii.gz',chosen),('overlap_uncertainty_T1w.nii.gz',uncertain),
                       ('assigned_cortical_ribbon_T1w.nii.gz',merged>0)]:
        full=np.zeros(raw.shape,np.int16);full[slices]=array;write_volume(out/'volumes'/name,full,raw)
    # Disjoint per-hemisphere products explicitly encode the adjudicated ownership.
    for h,array in [('L',np.where(merged<=124,merged,0)),('R',np.where(merged>=125,merged-124,0))]:
        full=np.zeros(raw.shape,np.int16);full[slices]=array;write_volume(out/f'volumes/{h}.MBNA124_resolved_T1w.nii.gz',full,raw)
    table(out/'metrics/hemisphere_adjudication.tsv',records)
    lut=pd.read_csv(ROOT/'source/volume_labels.tsv',sep='\t');counts=np.bincount(merged.ravel(),minlength=249)
    voxvol=abs(float(np.linalg.det(raw.affine[:3,:3])))
    table(out/'metrics/volume_parcels.tsv',[{'group':group,'volume_id':int(x.volume_id),'hemisphere':x.hemisphere,
        'source_id':int(x.source_id),'name':x.name,'voxels':int(counts[int(x.volume_id)]),'volume_mm3':float(counts[int(x.volume_id)]*voxvol)} for x in lut.itertuples() if x.volume_id>0])
    before=np.asanyarray(nib.load(out/'baseline/MBNA248_T1w.nii.gz').dataobj);after=np.asanyarray(nib.load(out/'volumes/MBNA248_T1w.nii.gz').dataobj)
    summary={'group':group,'volume_ROIs':int(np.count_nonzero(counts[1:])),'CGRSr_R_voxels':int(counts[246]),
        'overlap_voxels_before':int(conflict.sum()),'overlap_voxels_resolved':len(records),'unresolved_overlap_voxels':0,
        'adjudication_uncertain_voxels':int(uncertain.sum()),'changed_volume_voxels_vs_v3':int(np.count_nonzero(before!=after)),
        'primary_grid':list(raw.shape),'voxel_sizes_mm':list(map(float,raw.header.get_zooms())),
        'surface_geometry_modified':False,'contralateral_fill_used':False,'anatomical_release':'HOLD_SOURCE_GEOMETRY_AND_LOCAL_REVIEW'}
    dump(out/'metrics/volume_summary.json',summary)
    dump(out/'volumes/MBNA248_T1w.json',{**summary,'labels_table':'../../../source/volume_labels.tsv','background_or_unassigned':0,
        'left_ids':[1,124],'right_ids':[125,248],'conflict_resolution':'adult prior plus age anatomy; detailed audit and uncertainty mask retained',
        'native_surface_label_method':'v4 anatomy-constrained local mask closing, categorical boundary smoothing, source-part protection'})
    print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prior','volumes']);p.add_argument('--group',choices=GROUPS,required=True);a=p.parse_args()
    if a.stage=='prior':hemisphere_prior(a.group)
    else:volumes(a.group)
