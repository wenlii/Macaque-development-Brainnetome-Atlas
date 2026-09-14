"""Document, verify and seal this cosmetic DISPLAY ONLY companion package."""
import argparse,ast,csv,json,shutil,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from datetime import datetime,timezone
import numpy as np
import nibabel as nib
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from create_display import ROOT,BASE,ARCHIVE,GROUPS,sha,dump,table,mesh,vals


def rows(p):
    with p.open(encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))


def document():
    summaries=[]
    for group in GROUPS:
        out=ROOT/'groups'/group;summaries+=json.loads((out/'display_only.json').read_text())['summary']
        for res in ['dense','32k']:
            spec=ET.Element('CaretSpecFile',Version='1.0');ET.SubElement(spec,'MetaData')
            for h in ['L','R']:
                structure='CortexLeft' if h=='L' else 'CortexRight'
                entries=[('SURFACE',f'{h}.inflated.DISPLAY_ONLY.{res}.surf.gii'),('LABEL',f'{h}.MBNA124.DISPLAY_ONLY.{res}.label.gii')]
                if res=='32k':entries += [('METRIC',f'{h}.cortex.DISPLAY_ONLY.32k.shape.gii'),('METRIC',f'{h}.cosmetic_changes.DISPLAY_ONLY.32k.shape.gii')]
                for kind,path in entries:
                    assert (out/path).is_file();ET.SubElement(spec,'DataFile',Structure=structure,DataFileType=kind,Selected='true').text=path
            ET.indent(spec,space='  ');ET.ElementTree(spec).write(out/f'DISPLAY_ONLY_{res}.wb.spec',encoding='utf-8',xml_declaration=True)
        (out/'DISPLAY_ONLY_NOT_FOR_ANALYSIS.txt').write_text('DISPLAY ONLY. Cosmetic labels and masks: NOT for ROI extraction, morphometry, registration, statistics, or any downstream computation.\n仅供展示，不用于任何后续定量计算。分析请使用独立的 v4 分析包。\n',encoding='utf-8')
    table(ROOT/'qc/cosmetic_summary.tsv',summaries)
    copied=[]
    protected_archive=rows(ARCHIVE/'provenance/output_manifest_sha256.tsv');index={x['path']:x for x in protected_archive}
    for h in ['L','R']:
        for path in [f'source/{h}.MBNA124.32k.label.gii']+[f'groups/{g}/registration/{h}.geometry_32k.shape.gii' for g in GROUPS]:
            p=ARCHIVE/path;assert sha(p)==index[path]['sha256'],path
            copied.append({'path':str(p),'sha256':index[path]['sha256'],'role':'adult display-fill reference or original full-mesh geometry-caution audit'})
    table(ROOT/'provenance/additional_source_inputs.tsv',copied)
    source_readme='''# DISPLAY ONLY — cosmetic macaque atlas companion

This is a **cosmetically filled and smoothed visualization derivative**, created at the user's explicit request. It is **not an anatomical correction and must not be used for any downstream quantitative analysis**. Use the separate `macaque_age_atlas_v4_minimal` package for the analysis inputs. Its files and labels were not overwritten.

## View the result

- Open `figures/overview_DISPLAY_ONLY.png` for all six groups.
- Each `figures/Group*_DISPLAY_ONLY.png` contains left/right lateral and medial views.
- `figures/Group04_comparison_DISPLAY_ONLY.png` compares the unchanged v4 labels with this cosmetic result.
- `figures/Group04_detail_DISPLAY_ONLY.png` shows the formerly perforated regions at the same camera position and scale before and after processing.
- Open a group's `DISPLAY_ONLY_dense.wb.spec` in Connectome Workbench for the smooth presentation mesh. `DISPLAY_ONLY_32k.wb.spec` is the corresponding 32k alternative.

![Cosmetic presentation overview](figures/overview_DISPLAY_ONLY.png)

Every GIFTI and Workbench entry point contains `DISPLAY_ONLY` in its filename, and GIFTI metadata declares `IntendedUse=DISPLAY_ONLY_NOT_FOR_ANALYSIS`. The main presentation mesh has **519,842 vertices per hemisphere**; it is a two-step linear subdivision of the source 32,492-vertex inflated mesh. This is an interpolated display tessellation, not a new standard surface space or additional MRI resolution. Original 32k points retain their coordinates in the prefix of the dense mesh. Label IDs and the adult name/color dictionary are retained.

## Cosmetic modifications

The procedure fills enclosed background components while retaining the largest medial-wall component; closes narrow gaps with a 2-mm graph radius measured on source midthickness; smooths the display-mask signed distance; and protects the deep main medial-wall core. Filled locations initially use neighboring labels and, where available, the common adult source label. Eight iterations of area-weighted categorical score diffusion make label borders less jagged. A seed from each existing ROI component preserves all 124 label identities in each hemisphere. Integer IDs are not averaged.

The dense labels use interpolated categorical scores and the smoothed display coverage field. This makes boundaries finer than applying a single color to an entire original 32k face. Gray medial wall is rendered without lighting; colored cortex is gently shaded. Group04's previously retained grey holes are filled **including positions with original geometric caution flags**, solely to create a clean illustration.

These edits can change apparent ROI size, location, support and component structure. They cannot be treated as measurements, independent age-specific boundaries, or proof that original geometry is repaired. No T1-space volume atlas is produced for this display-only version. The scientific v4 limitations and anatomical HOLD status remain unchanged.

## Audit and regeneration

`qc/cosmetic_summary.tsv` and per-hemisphere `cosmetic_vertex_changes.tsv` record changes relative to v4. Cosmetic flag bits on the 32k grid are 1=filled formerly unlabeled point, 2=other label change, 4=fill over a former geometry-caution location. These are edit flags, not confidence values.

The source analysis-package hash manifest is protected under `provenance/`. `dataset_description.json` and the processing contract explicitly prohibit quantitative use. The companion source-code and input hashes permit inspection of how the images were made; no attached screenshot was painted or retouched.

With numpy, scipy, nibabel, pyvista and matplotlib installed, run `python code/render_display.py group --group Group04_59p5-68p0` or `python code/render_display.py overview` to render the existing display files. Rebuilding the cosmetic data with `code/create_display.py` requires the separate original v4 packages at the locations documented in the processing contract. Only Windows was tested here.

`code/finalize_display.py verify` checks ROI identities, mesh/label compatibility, the absence of isolated gray holes, explicit usage tags and protected v4 hashes when the source package is present. `manifest_sha256.tsv` seals this display package, excluding itself and Python bytecode caches.

The original project license notice and citation metadata are retained from the analysis package, with a display-only title/version. No public upload or quota reset card was used.
'''
    (ROOT/'README.md').write_text(source_readme,encoding='utf-8')
    zh='''# 仅供展示：已填洞和平滑的版本

**此版本只用于图片、幻灯片和界面展示，不用于任何后续计算。** 它主动修改了显示标签、覆盖范围和边界，包括原本带几何标记的位置；不代表这些区域已经得到解剖修复。原 v4 分析包保持不变。

## 从哪里看

- [六组总览](figures/overview_DISPLAY_ONLY.png)
- [Group04 前后对比](figures/Group04_comparison_DISPLAY_ONLY.png)
- [Group04 局部放大](figures/Group04_detail_DISPLAY_ONLY.png)
- `figures/` 下另有六组各自的四视图。
- 在 Workbench 中打开每组 `DISPLAY_ONLY_dense.wb.spec`。如需较小网格，可打开 `DISPLAY_ONLY_32k.wb.spec`。

圈出的独立灰洞已在展示版中填平，狭窄凹口与彩色标签边界经过外观平滑，主体内侧壁仍保留为灰色。左右各 124 个标签保留。主展示网格每侧 519,842 个顶点，由 32k 网格插值得到，只提高显示细腻度，没有增加 MRI 信息。

## 使用边界

所有显示标签／表面／入口名称均带 `DISPLAY_ONLY`，文件元数据也写明用途。**不能用于 ROI 信号提取、面积或厚度计算、连接矩阵、配准、统计或其它定量分析。** 分析仍使用独立的 `macaque_age_atlas_v4_minimal`，不将本展示版的标签或掩膜混入其中。

灰洞填补和边界平滑属于视觉加工。原来的几何可疑区域未被修复，分析版仍保留原有 HOLD 状态。本展示版不提供模板空间体积图谱，以免与分析结果混用。

[改动汇总](qc/cosmetic_summary.tsv) 和每组逐顶点表记录具体修改。32k 改动标记为：1=原零标签被补齐，2=其它标签变化，4=补齐了原先带几何可疑标记的位置，按位相加。完整方法及重绘方式见 [README](README.md)。

原 v4 分析标签、体积、QC 和封存校验值均保持不变；未对外上传，未使用额度重置卡。
'''
    (ROOT/'使用说明_仅供展示.md').write_text(zh,encoding='utf-8')
    (ROOT/'DISPLAY_ONLY_NOT_FOR_ANALYSIS.txt').write_text('DISPLAY ONLY / NOT FOR ANY DOWNSTREAM COMPUTATION\n仅供展示，不得用于后续计算。原 v4 分析文件保持不变。\n',encoding='utf-8')
    shutil.copyfile(BASE/'LICENSE.md',ROOT/'LICENSE.md')
    cff=(BASE/'CITATION.cff').read_text(encoding='utf-8').replace('version: v4-minimal','version: v4-display-only')
    cff=cff.replace('title: Macaque Age-Template MBNA Atlas - Minimal 32k Surface and Volume Package','title: Macaque Atlas Cosmetic Surface Presentation - DISPLAY ONLY')
    (ROOT/'CITATION.cff').write_text(cff,encoding='utf-8')
    (ROOT/'requirements.txt').write_text('numpy\nscipy\nnibabel\npyvista\nmatplotlib\n',encoding='utf-8')
    dump(ROOT/'dataset_description.json',{'Name':'Macaque v4 cosmetic surface presentation - DISPLAY ONLY','Version':'v4-display-only','IntendedUse':'DISPLAY_ONLY_NOT_FOR_ANALYSIS',
        'AllowedUse':['illustration','slides','visual demonstration'],'QuantitativeUseAllowed':False,'CosmeticLabels':True,'AnalysisSource':'macaque_age_atlas_v4_minimal',
        'AnalysisFilesModified':False,'AnatomicalCorrectionClaim':False,'OriginalGeometryCautionsOverriddenForDisplay':True,'VolumeAtlasIncluded':False,
        'VerticesPerHemisphere':{'32k':32492,'dense_display':519842},'Groups':GROUPS,'Panels':'all image filenames and headers label the display-only purpose',
        'InterpolationAddsMRIInformation':False,'SourceAnalysisStatus':'HOLD_SOURCE_GEOMETRY_AND_LOCAL_REVIEW'})
    print('Display documentation and Workbench specs written',flush=True)


def background_components(labels,f):
    ids=np.flatnonzero(labels==0);mapping=np.full(len(labels),-1,np.int32);mapping[ids]=np.arange(len(ids))
    pairs=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);pairs=pairs[(labels[pairs]==0).all(1)]
    e=mapping[pairs];a=sparse.csr_matrix((np.ones(len(e)*2,np.uint8),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(ids),len(ids)))
    return int(connected_components(a,directed=False,return_labels=False))


def verify(write_report=False):
    checks=[]
    def ck(ok,name,g='',h='',res=''):
        checks.append({'group':g,'hemisphere':h,'resolution':res,'check':name,'pass':bool(ok)});assert ok,(name,g,h,res)
    for group in GROUPS:
        out=ROOT/'groups'/group
        for h in ['L','R']:
            coarse=vals(out/f'{h}.MBNA124.DISPLAY_ONLY.32k.label.gii').astype(int)
            mask=vals(out/f'{h}.cortex.DISPLAY_ONLY.32k.shape.gii')>0
            ck(np.array_equal(coarse>0,mask),'32k display mask agrees with labels',group,h,'32k')
            for res,npoints in [('32k',32492),('dense',519842)]:
                v,f=mesh(out/f'{h}.inflated.DISPLAY_ONLY.{res}.surf.gii');im=nib.load(out/f'{h}.MBNA124.DISPLAY_ONLY.{res}.label.gii');label=im.agg_data().astype(int)
                ck(len(v)==len(label)==npoints and f.min()>=0 and f.max()<npoints,'mesh and label lengths match',group,h,res)
                ck(set(label[label>0])==set(range(1,125)),'all 124 labels retained',group,h,res)
                ck(im.meta.get('IntendedUse')=='DISPLAY_ONLY_NOT_FOR_ANALYSIS','explicit non-analysis metadata',group,h,res)
                ck(background_components(label,f)==1,'no isolated gray holes: only the main medial-wall component remains',group,h,res)
                if res=='dense':ck(np.array_equal(label[:32492],coarse),'dense mesh retains the original display vertices',group,h,res)
            for res in ['32k','dense']:
                spec=ET.parse(out/f'DISPLAY_ONLY_{res}.wb.spec')
                for element in spec.findall('DataFile'):ck((out/element.text).is_file() and 'DISPLAY_ONLY' in element.text,'relative display-only entry exists',group,h,res)
        print('Verified DISPLAY ONLY',group,flush=True)
    flagged=np.array([2055,2106,2107,2157,2206,2207])
    ck(np.all(vals(ROOT/'groups/Group04_59p5-68p0/R.MBNA124.DISPLAY_ONLY.32k.label.gii')[flagged]>0),'Group04 originally flagged six-point hole cosmetically filled')
    sources_available=BASE.exists()
    if BASE.exists():
        for row in rows(ROOT/'provenance/protected_v4_manifest.tsv'):ck(sha(BASE/row['path'])==row['sha256'],'protected v4 analysis file unchanged: '+row['path'])
    if (ROOT/'provenance/additional_source_inputs.tsv').exists():
        for row in rows(ROOT/'provenance/additional_source_inputs.tsv'):
            p=Path(row['path'])
            if p.exists():ck(sha(p)==row['sha256'],'additional source unchanged')
            else:sources_available=False
    ck(not list(ROOT.rglob('*.nii*')),'no analysis volume included')
    for p in (ROOT/'code').glob('*.py'):ast.parse(p.read_text(encoding='utf-8'))
    result={'status':'PASS','checks':len(checks),'all_hemispheres_have_124_labels':True,'background_components_per_hemisphere':1,'Group04_target_holes_filled_for_display':True,
            'intended_use':'DISPLAY_ONLY_NOT_FOR_ANALYSIS','analysis_files_modified':False,'source_integrity_check':'PASS' if sources_available else 'external source unavailable'}
    if write_report:
        table(ROOT/'qc/technical_checks.tsv',checks);dump(ROOT/'qc/technical_validation.json',result)
    print(result,flush=True)


def seal():
    validation=json.loads((ROOT/'qc/technical_validation.json').read_text());visual=json.loads((ROOT/'qc/visual_review.json').read_text())
    assert validation['status']=='PASS' and visual['review_complete']
    for x in visual['files']:assert sha(ROOT/x['path'])==x['sha256']
    dump(ROOT/'provenance/final_receipt.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),'completion':'COMPLETE_FOR_ALL_SIX_GROUPS','intended_use':'DISPLAY_ONLY_NOT_FOR_ANALYSIS',
        'technical_validation':validation,'images_reviewed':len(visual['files']),'analysis_v4_preserved':True,'quota_reset_cards_used':False,'external_upload':False})
    manifest=ROOT/'manifest_sha256.tsv';files=[p for p in sorted(ROOT.rglob('*')) if p.is_file() and p!=manifest and '__pycache__' not in p.parts]
    records=[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in files];table(manifest,records)
    for row in records:assert sha(ROOT/row['path'])==row['sha256']
    print({'sealed_files':len(records),'total_files_including_manifest':len(records)+1,'payload_bytes':sum(x['bytes'] for x in records),'root':str(ROOT)},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['document','verify','seal']);p.add_argument('--write-report',action='store_true');a=p.parse_args()
    {'document':document,'verify':lambda:verify(a.write_report),'seal':seal}[a.stage]()
