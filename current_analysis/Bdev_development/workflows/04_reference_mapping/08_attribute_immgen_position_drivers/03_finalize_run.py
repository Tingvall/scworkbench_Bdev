#!/usr/bin/env python3
from __future__ import annotations
import hashlib, json
from pathlib import Path
import sys, yaml

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

cfg_path=Path(sys.argv[2]).resolve(); assert sys.argv[1]=='--config'
cfg=yaml.safe_load(cfg_path.read_text()); root=Path(cfg['storage']['results_root'])
summary=json.loads((root/'03_reports/02_position_driver_summary.json').read_text())
required=[
 '01_tables/01_position_driver_latent_coordinates.tsv', '01_tables/02_provisional_label_diagnostic_cells.tsv.gz',
 '01_tables/04_peak_driver_candidates.tsv.gz', '01_tables/09_library_label_position_diagnostic.tsv',
 '01_tables/10_library_latent_dimension_drivers.tsv', '01_tables/13_technical_bias_associations.tsv',
 '01_tables/14_recurring_peak_driver_candidates.tsv.gz', '01_tables/15_peak_axis_prevalence_audit.tsv.gz',
 '01_tables/16_latent_dimension_technical_associations.tsv', '02_figures/01_latent_dimension_pull.png',
 '02_figures/02_expected_vs_observed_position.png', '02_figures/03_technical_bias_associations.png']
for rel in required: assert (root/rel).is_file(), rel
assert summary['libraries_evaluated']==28 and summary['reference_cells']==1382
files=[]
for p in sorted(x for x in root.rglob('*') if x.is_file()):
    if p.name=='03_run_file_inventory.json': continue
    files.append({'path':str(p.relative_to(root)), 'bytes':p.stat().st_size, 'sha256':digest(p)})
out=root/'06_manifests/03_run_file_inventory.json'; out.write_text(json.dumps({'schema_version':'bdev.run_file_inventory.v1','run_identity':cfg['run_identity'],'files':files},indent=2)+'\n')
print('status',summary['status']); print('files',len(files)); print('inventory',out)
