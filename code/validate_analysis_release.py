"""Validate a relocated minimal atlas without access to the source project.

Default: hashes, file inventory, internal references, GIFTI topology and labels,
NIfTI grids, ROI counts, and compact QC tables. --hash-only needs only Python.
This validator reads the package and prints JSON; it never edits package files.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(4194304),b''):h.update(block)
    return h.hexdigest()


def rows(path):
    with path.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def validate(root, hash_only=False):
    root=root.resolve();checks=0
    def ck(ok,message):
        nonlocal checks
        checks+=1
        if not ok:raise AssertionError(message)
    def internal(base,relative):
        p=(base/relative).resolve()
        ck(p.is_relative_to(root),'Reference escapes release: '+relative)
        ck(p.is_file(),'Missing reference: '+relative)
        return p
    manifest=rows(root/'manifest_sha256.tsv');seen=set()
    for row in manifest:
        ck(row['path'] not in seen,'Duplicate manifest entry');seen.add(row['path'])
        p=internal(root,row['path'])
        ck(p.stat().st_size==int(row['bytes']) and sha(p)==row['sha256'],'Hash mismatch: '+row['path'])
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='manifest_sha256.tsv'}
    ck(seen==actual,'Manifest inventory does not match payload')
    description=json.loads((root/'dataset_description.json').read_text(encoding='utf-8'))
    ages=rows(root/'age_groups.tsv');groups=[x['group'] for x in ages]
    ck(len(groups)==6 and len(set(groups))==6,'Expected six unique age groups')
    ck(groups==description['groups'],'Dataset description group mismatch')
    ck(description['anatomical_release']=='HOLD_SOURCE_GEOMETRY_AND_LOCAL_REVIEW','Anatomical status changed')
    lut=rows(root/'volume_labels.tsv');positive={int(x['volume_id']):x for x in lut if int(x['volume_id'])>0}
    ck(set(positive)==set(range(1,249)),'Volume label dictionary must have 248 labels')
    group_qc={x['group']:x for x in rows(root/'qc/group_summary.tsv')}
    roi_qc=rows(root/'qc/roi_qc.tsv');roi_lookup={(x['group'],x['hemisphere'],int(x['source_id'])):x for x in roi_qc}
    ck(len(roi_qc)==1488 and len(roi_lookup)==1488,'ROI QC must have 6 x 248 unique rows')
    for group in groups:
        folder=root/'groups'/group
        spec=ET.parse(folder/'atlas.wb.spec').getroot()
        for item in spec.findall('DataFile'):internal(folder,item.text)
        sidecar=json.loads((folder/'MBNA248_T1w.json').read_text(encoding='utf-8'))
        ck(sidecar['group']==group,'Sidecar group mismatch')
        for field in ['labels_table','reference_template','ownership_uncertainty_mask','group_qc_table','roi_qc_table']:
            internal(folder,sidecar[field])
        ck(sidecar['volume_ROIs']==248,'Primary volume is incomplete')
    if hash_only:
        return {'status':'PASS','mode':'hashes_inventory_and_references','checks':checks,'manifest_payload_files':len(manifest),'groups':6,'root':str(root)}
    try:
        import numpy as np
        import nibabel as nib
    except ImportError as e:
        raise RuntimeError('Full checks require numpy and nibabel; install requirements.txt or use --hash-only.') from e
    for group in groups:
        folder=root/'groups'/group;q=group_qc[group]
        t1=nib.load(folder/'T1w.nii.gz');volume=nib.load(folder/'MBNA248_T1w.nii.gz')
        uncertain=nib.load(folder/'overlap_uncertainty_T1w.nii.gz')
        data=np.asanyarray(volume.dataobj);u=np.asanyarray(uncertain.dataobj)
        voxel_counts=np.bincount(data.ravel().astype(int),minlength=249)
        ck(volume.shape==t1.shape==uncertain.shape==(256,256,256),'Unexpected volume grid')
        ck(np.array_equal(volume.affine,t1.affine) and np.array_equal(volume.affine,uncertain.affine),'Volume affine mismatch')
        ck(np.allclose(volume.header.get_zooms(),[.5,.5,.5],atol=1e-5),'Expected original 0.5-mm grid')
        ck(np.issubdtype(data.dtype,np.integer) and set(np.unique(data[data>0]))==set(range(1,249)),'Expected 248 integer volume labels')
        ck(np.isin(u,[0,1]).all() and not np.any((u>0)&(data==0)),'Invalid ownership uncertainty mask')
        ck(int(np.count_nonzero(u))==int(q['uncertain_ownership_voxels']),'Uncertainty count mismatch')
        for form in ['qform','sform']:
            ck(np.array_equal(getattr(volume,'get_'+form)(),getattr(t1,'get_'+form)()),'NIfTI '+form+' changed')
        for h,offset in [('L',0),('R',124)]:
            mid=nib.load(folder/f'{h}.midthickness.32k.surf.gii');infl=nib.load(folder/f'{h}.inflated.32k.surf.gii')
            vv=mid.agg_data('pointset');faces=mid.agg_data('triangle');iv=infl.agg_data('pointset')
            ck(vv.shape==iv.shape==(32492,3) and np.isfinite(vv).all() and np.isfinite(iv).all(),'Invalid 32k surface')
            ck(np.array_equal(faces,infl.agg_data('triangle')) and faces.min()>=0 and faces.max()<32492,'Incompatible surface topology')
            labelim=nib.load(folder/f'{h}.MBNA124.32k.label.gii');label=labelim.agg_data()
            cortex=nib.load(folder/f'{h}.cortex.32k.shape.gii').agg_data();flags=nib.load(folder/f'{h}.geometry_flags.32k.shape.gii').agg_data()
            ck(label.shape==cortex.shape==flags.shape==(32492,),'Labels and metrics need matching 32k vertex order')
            ck(np.issubdtype(label.dtype,np.integer) and set(np.unique(label[label>0]))==set(range(1,125)),'Expected 124 integer hemisphere labels')
            ck(np.isin(cortex,[0,1]).all() and not np.any(label[cortex==0]),'Labels outside valid cortex')
            ck(not np.any(flags[cortex==0]),'Geometry flag outside published cortex')
            ck(int(np.count_nonzero((cortex>0)&(label==0)))==int(q[f'{h}_unassigned_vertices_32k']),'Unassigned vertex count mismatch')
            labeltable={x.key:x for x in labelim.labeltable.labels}
            vertex_counts=np.bincount(label.astype(int),minlength=125)
            flagged_counts=np.bincount(label[flags>0].astype(int),minlength=125)
            for source_id in range(1,125):
                roi=roi_lookup[group,h,source_id];entry=positive[source_id+offset]
                ck(entry['hemisphere']==h and int(entry['source_id'])==source_id,'Hemisphere numbering mismatch')
                ck(labeltable[source_id].label==entry['name'],'Label name differs from shared dictionary')
                ck(np.allclose([labeltable[source_id].red,labeltable[source_id].green,labeltable[source_id].blue],[float(entry[x]) for x in ['red','green','blue']],atol=1e-6),'Label colors differ')
                ck(int(vertex_counts[source_id])==int(roi['vertices_32k']),'ROI vertex count mismatch')
                ck(int(flagged_counts[source_id])==int(roi['geometry_flagged_vertices_32k']),'ROI geometry QC mismatch')
                ck(int(voxel_counts[source_id+offset])==int(roi['volume_voxels']),'ROI voxel count mismatch')
        ck(int(np.count_nonzero(data==246))==int(q['right_CGRSr_voxels']),'CGRSr count mismatch')
    return {'status':'PASS','mode':'full','checks':checks,'manifest_payload_files':len(manifest),'groups':6,
            'ROIs_per_group':{'left_32k':124,'right_32k':124,'volume':248},'anatomical_release':description['anatomical_release'],'root':str(root)}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',nargs='?',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--hash-only',action='store_true');args=p.parse_args()
    try:print(json.dumps(validate(args.root,args.hash_only),indent=2,ensure_ascii=False))
    except Exception as e:
        print(json.dumps({'status':'FAIL','error':str(e)},indent=2,ensure_ascii=False));sys.exit(1)
