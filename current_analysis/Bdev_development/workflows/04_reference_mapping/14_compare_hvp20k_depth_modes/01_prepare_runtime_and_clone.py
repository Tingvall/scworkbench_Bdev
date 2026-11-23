#!/usr/bin/env python3
from __future__ import annotations

import argparse, hashlib, json, subprocess
from pathlib import Path
import yaml

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()

def main() -> None:
    p = argparse.ArgumentParser(); p.add_argument("--config", required=True, type=Path)
    config_path = p.parse_args().config.resolve(); cfg = yaml.safe_load(config_path.read_text())
    rt = cfg["runtime"]; source = Path(cfg["inputs"]["source_bundle"]["path"])
    if sha256(source / "manifest.json") != cfg["inputs"]["source_bundle"]["manifest_sha256"]: raise RuntimeError("source manifest changed")
    if sha256(Path(rt["sif"])) != rt["sif_sha256"]: raise RuntimeError("SIF changed")
    run_lib = Path(rt["run_library"]); run_lib.mkdir(parents=True, exist_ok=True)
    results = Path(cfg["storage"]["results_root"])
    for name in ("01_tables","02_figures","03_reports","04_logs","05_objects","06_manifests"): (results/name).mkdir(parents=True, exist_ok=True)
    Path(cfg["storage"]["sessions"]).mkdir(parents=True, exist_ok=True)
    packages=[]
    for key in ("scworkbench","scwbATAC","scwbLatent"):
        rec=rt[key]; path=Path(rec["source_path"])
        if git(path,"rev-parse","HEAD") != rec["source_commit"]: raise RuntimeError(f"{key} commit changed")
        if git(path,"status","--porcelain"): raise RuntimeError(f"{key} source is dirty")
        manifest=path/"inst"/"scworkbench"/"operation-pack.json"
        if "manifest_sha256" in rec and sha256(manifest) != rec["manifest_sha256"]: raise RuntimeError(f"{key} manifest changed")
        packages.append((key,path))
    rlibs=f"{run_lib}:{rt['shared_development_library']}"
    for _,path in packages:
        subprocess.run([rt["apptainer"],"exec","--cleanenv","--bind","/mnt/sata2:/mnt/sata2","--bind","/home/jtg/msi_projects/scworkbench_code:/scworkbench_code","--env",f"R_LIBS_USER={rlibs}",rt["sif"],"R","CMD","INSTALL",f"--library={run_lib}",f"/scworkbench_code/{path.name}"],check=True)
    clone=Path(cfg["storage"]["selection_clone"])
    if clone.exists():
        if sha256(clone/"manifest.json") != cfg["inputs"]["source_bundle"]["manifest_sha256"]: raise RuntimeError("existing clone is not pristine")
    else:
        clone.parent.mkdir(parents=True,exist_ok=True)
        subprocess.run(["cp","--recursive","--reflink=auto","--",str(source),str(clone)],check=True)
    receipt={"schema_version":"bdev.hvp20k_runtime_clone.v1","status":"complete","run_identity":cfg["run_identity"],"config_path":str(config_path),"config_sha256":sha256(config_path),"source_manifest_sha256":sha256(source/"manifest.json"),"sif_sha256":sha256(Path(rt["sif"])),"packages":[{"package":k,"commit":git(v,"rev-parse","HEAD")} for k,v in packages]}
    (results/"06_manifests"/"01_runtime_and_clone_receipt.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")

if __name__ == "__main__": main()
