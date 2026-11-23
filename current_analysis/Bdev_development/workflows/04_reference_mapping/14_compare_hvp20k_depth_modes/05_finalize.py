#!/usr/bin/env python3
import argparse,hashlib,json
from pathlib import Path
import yaml
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def load(p):
 x=yaml.safe_load(p.read_text());b=x.pop('base_config',None)
 if b is None:return x
 y=load((p.parent/b).resolve())
 for k,v in x.items():
  if isinstance(v,dict) and isinstance(y.get(k),dict):y[k].update(v)
  else:y[k]=v
 return y
p=argparse.ArgumentParser();p.add_argument('--config',required=True,type=Path);cp=p.parse_args().config.resolve();c=load(cp);r=Path(c['storage']['results_root']);out=r/'06_manifests'/'05_run_inventory.json';files=[{'path':str(x.relative_to(r)),'bytes':x.stat().st_size,'sha256':sha(x)} for x in sorted(r.rglob('*')) if x.is_file() and x!=out];out.write_text(json.dumps({'schema_version':'bdev.hvp20k_inventory.v1','status':'complete','run_identity':c['run_identity'],'config_sha256':sha(cp),'current_default_changed':False,'files':files},indent=2,sort_keys=True)+'\n')
