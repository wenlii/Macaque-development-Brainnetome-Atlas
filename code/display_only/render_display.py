"""DISPLAY ONLY: polished surface plates from cosmetic GIFTI labels.

Figure outline: six-group overview; four views per group; a Group04 before/after
and magnified confirmation. All panels are drawn from meshes, not retouched PNGs.
Original adult ROI colors are retained; the main medial wall is flat neutral gray.
"""
import argparse,gc
from pathlib import Path
import numpy as np
import nibabel as nib
import pyvista as pv
import matplotlib as mpl
from create_display import ROOT,BASE,GROUPS,mesh,vals

mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':10,'pdf.fonttype':42,'svg.fonttype':'none'})


def make_mesh(group,h,analysis=False):
    if analysis:
        out=BASE/'groups'/group
        v,f=mesh(out/f'{h}.inflated.32k.surf.gii');im=nib.load(out/f'{h}.MBNA124.32k.label.gii')
    else:
        out=ROOT/'groups'/group
        v,f=mesh(out/f'{h}.inflated.DISPLAY_ONLY.dense.surf.gii');im=nib.load(out/f'{h}.MBNA124.DISPLAY_ONLY.dense.label.gii')
    label=im.agg_data().astype(int);pal=np.full((125,3),.83)
    for item in im.labeltable.labels:
        if item.key:pal[item.key]=[item.red,item.green,item.blue]
    lv=label[f];face_label=np.where(lv[:,0]==lv[:,1],lv[:,0],np.where(lv[:,1]==lv[:,2],lv[:,1],lv[:,0]))
    poly=pv.PolyData(v,np.c_[np.full(len(f),3,np.int32),f].ravel())
    poly.cell_data['RGB']=(pal[face_label]*255).astype(np.uint8)
    poly.compute_normals(cell_normals=False,point_normals=True,inplace=True)
    color=poly.extract_cells(np.flatnonzero(face_label>0));wall=poly.extract_cells(np.flatnonzero(face_label==0))
    return {'v':v,'faces':f,'labels':label,'colored':color,'wall':wall}


def panel(plot,item,h,view,title,detail=False,focus=None):
    plot.add_mesh(item['colored'],scalars='RGB',rgb=True,preference='cell',smooth_shading=True,ambient=.68,diffuse=.32,specular=0,show_scalar_bar=False)
    plot.add_mesh(item['wall'],color=[.83,.83,.83],lighting=False,show_scalar_bar=False)
    v=item['v'];center=(v.min(0)+v.max(0))/2
    direction=-1 if h=='L' else 1
    if view=='medial':direction*=-1
    if detail:
        # Original 32k points occupy the prefix in the subdivided display mesh.
        points=focus if focus is not None else v[np.isin(item['labels'],[122,123])]
        center=(points.min(0)+points.max(0))/2
    plot.camera_position=[center+np.array([direction*180.,0,0]),center,[0,0,1]]
    plot.enable_parallel_projection()
    viewport=plot.renderer.GetViewport();width,height=plot.window_size
    aspect=width*(viewport[2]-viewport[0])/(height*(viewport[3]-viewport[1]))
    if detail:
        plot.camera.parallel_scale=max(np.ptp(points[:,1])/(2*aspect),np.ptp(points[:,2])/2)*1.5
    else:plot.camera.parallel_scale=max(np.ptp(v[:,1])/(2*aspect),np.ptp(v[:,2])/2)*1.16
    actor=plot.add_text(title,font_size=12,color='#222222',font='arial',position='upper_left')
    # Keep the usage label readable on close-up images.
    actor.GetTextProperty().SetBackgroundColor(1,1,1);actor.GetTextProperty().SetBackgroundOpacity(1)


def setup(shape,size):
    p=pv.Plotter(shape=shape,off_screen=True,window_size=size,border=False)
    # Avoid the SSAA framebuffer pass: it corrupts multiple viewports on this host.
    p.set_background('white',all_renderers=True);p.ren_win.SetMultiSamples(4)
    return p


def group_plate(group):
    data={h:make_mesh(group,h) for h in ['L','R']}
    p=setup((1,4),(3600,950))
    for col,(h,view) in enumerate([('L','lateral'),('L','medial'),('R','medial'),('R','lateral')]):
        p.subplot(0,col);panel(p,data[h],h,view,f'{group.split("_")[0]} | {h} {view}\nDISPLAY ONLY')
    p.screenshot(ROOT/'figures'/f'{group}_DISPLAY_ONLY.png');p.close();del data;gc.collect()
    print('Rendered display plate',group,flush=True)


def group04_comparison():
    group=GROUPS[3]
    before={h:make_mesh(group,h,True) for h in ['L','R']};after={h:make_mesh(group,h) for h in ['L','R']}
    for detail in [False,True]:
        p=setup((2,2),(2200,1550))
        for row,h in enumerate(['L','R']):
            focus=before[h]['v'][np.isin(before[h]['labels'],[122,123])]
            p.subplot(row,0);panel(p,before[h],h,'medial',f'{h} medial | unchanged v4 analysis',detail,focus)
            p.subplot(row,1);panel(p,after[h],h,'medial',f'{h} medial | cosmetic DISPLAY ONLY',detail,focus)
        name='Group04_detail_DISPLAY_ONLY.png' if detail else 'Group04_comparison_DISPLAY_ONLY.png'
        p.screenshot(ROOT/'figures'/name);p.close()
    print('Rendered Group04 display comparisons',flush=True)


def overview():
    p=setup((6,4),(3600,4200))
    # Keep one shared mesh pair per group for the two viewing directions.
    keep=[]
    for row,group in enumerate(GROUPS):
        pair={h:make_mesh(group,h) for h in ['L','R']};keep.append(pair)
        for col,(h,view) in enumerate([('L','lateral'),('L','medial'),('R','medial'),('R','lateral')]):
            p.subplot(row,col)
            title=f'{group}\n{h} {view} | DISPLAY ONLY' if col==0 else f'{h} {view}'
            panel(p,pair[h],h,view,title)
    p.screenshot(ROOT/'figures/overview_DISPLAY_ONLY.png');p.close();del keep;gc.collect()
    print('Rendered six-group DISPLAY ONLY overview',flush=True)


# PNG is the approved presentation output; vector exports are not generated.
# fig.savefig('display_only.pdf', bbox_inches='tight')
# fig.savefig('display_only.svg', bbox_inches='tight')
# For matplotlib-based companion panels: fig.savefig('display_only.png', dpi=600).

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['group','group04','overview']);p.add_argument('--group',choices=GROUPS);a=p.parse_args()
    {'group':lambda:group_plate(a.group),'group04':group04_comparison,'overview':overview}[a.stage]()
