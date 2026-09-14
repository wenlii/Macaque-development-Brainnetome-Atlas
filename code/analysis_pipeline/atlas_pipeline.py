"""Shared I/O for the approved v4 label repair; all v3 inputs remain protected."""
from pathlib import Path
import os, json, csv, hashlib, subprocess, copy
from datetime import datetime, timezone
import numpy as np
import nibabel as nib

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
BASELINE = PROJECT / "result/age_native_MBNA_v3"
GROUPS = ["Group01_34p0-42p5", "Group02_42p5-51p0", "Group03_51p0-59p5", "Group04_59p5-68p0", "Group05_68p0-76p5", "Group06_76p5-inf"]
HEMIS = {"L": "lh", "R": "rh"}
WB = Path(r"D:\software\workbench\workbench\bin_windows64\wb_command.exe")

def dump(path, obj):
    def cv(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x))
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=cv)+"\n", encoding="utf-8")

def table(path, rows, columns=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    columns = columns or (list(rows[0]) if rows else ["no_records"])
    with path.open("w",newline="",encoding="utf-8") as fp:
        w=csv.DictWriter(fp,fieldnames=columns,delimiter="\t");w.writeheader();w.writerows(rows)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fp:
        for b in iter(lambda:fp.read(4194304),b""):h.update(b)
    return h.hexdigest()

def freeze(path):
    p=ROOT/"provenance/new_input_manifest.json";m=json.loads(p.read_text()) if p.exists() else {}
    key=str(Path(path).resolve());rec={"bytes":Path(path).stat().st_size,"sha256":sha(path)}
    if key in m: assert m[key]==rec, f"Input changed: {path}"
    m[key]=rec;dump(p,m)

def command(logstem,args):
    logstem=Path(logstem);logstem.parent.mkdir(parents=True,exist_ok=True)
    argv=[str(WB)]+[str(x) for x in args]
    p=subprocess.run(argv,capture_output=True,text=True,errors="replace",env=dict(os.environ,OMP_NUM_THREADS="4"))
    logstem.with_suffix(".log").write_text(p.stdout+p.stderr,encoding="utf-8")
    dump(logstem.with_suffix(".json"),{"argv":argv,"exit_code":p.returncode,"utc":datetime.now(timezone.utc).isoformat()})
    if p.returncode:raise RuntimeError(str(logstem)+" failed")

def mesh(path):
    im=nib.load(path);return im.agg_data("pointset").astype(float),im.agg_data("triangle").astype(np.int32)

def values(path):return np.asarray(nib.load(path).agg_data()).ravel()

def write_metric(path,data,h,name):
    im=nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.asarray(data,np.float32),intent="NIFTI_INTENT_SHAPE")])
    im.meta["AnatomicalStructurePrimary"]="CortexLeft" if h=="L" else "CortexRight";im.darrays[0].meta["Name"]=name;nib.save(im,path)

def write_labels(path,data,h,source=None,labeltable=None):
    im=nib.gifti.GiftiImage(darrays=[nib.gifti.GiftiDataArray(np.asarray(data,np.int32),intent="NIFTI_INTENT_LABEL")])
    im.labeltable=copy.deepcopy(labeltable if labeltable is not None else nib.load(source).labeltable)
    im.meta["AnatomicalStructurePrimary"]="CortexLeft" if h=="L" else "CortexRight"
    im.darrays[0].meta["Name"]="Source-informed MBNA boundary repair v4";nib.save(im,path)

def write_volume(path,data,ref,label=True):
    hdr=ref.header.copy();hdr.set_data_dtype(np.int16 if label else np.float32)
    hdr.set_intent("label" if label else "none")
    im=nib.Nifti1Image(np.asarray(data,np.int16 if label else np.float32),ref.affine,hdr)
    im.set_qform(ref.get_qform(),int(ref.header["qform_code"]));im.set_sform(ref.get_sform(),int(ref.header["sform_code"]));nib.save(im,path)

def vertex_areas(v,f):
    area=.5*np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)
    return np.bincount(f.ravel(),weights=np.repeat(area/3,3),minlength=len(v))

def signed_sphere_triangles(v,f):return np.einsum("ij,ij->i",v[f[:,0]],np.cross(v[f[:,1]],v[f[:,2]]))
