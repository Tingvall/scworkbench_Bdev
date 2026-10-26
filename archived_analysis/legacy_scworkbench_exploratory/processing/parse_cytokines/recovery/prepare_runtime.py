#!/usr/bin/env python3
"""Run inside the supplied SIF; create an isolated, explicitly patched engine."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

repo, destination = map(Path, sys.argv[1:])
source = Path('/opt/scworkbench')
target = destination / 'source'
if target.exists():
    raise SystemExit(f'Refusing to overwrite recovery source: {target}')
original = source / 'R/bundle_init.R'
def sha256(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()
expected = '8979a3300a017bfc0b3fd0e78239c4b622a9ee3049f6c79491c5c85f8a995828'
if sha256(original) != expected:
    raise SystemExit('Source identity differs from the inspected scWorkbench build')
row_source_sha256 = '6c3bd455e140684ef282f2d947b44121e6f8cf32ffebe2534f02ce4b1bdfbe97'
if sha256(source / 'R/provenance_json.R') != row_source_sha256:
    raise SystemExit('Row-key source identity differs from the inspected build')
shutil.copytree(source, target, ignore=shutil.ignore_patterns(
    '.git', '.codex', '.Rhistory', '.tmp_exports_overview.csv', 'archive', 'third_party'))
patch = repo / 'processing/parse_cytokines/recovery/registry_json.patch'
helper = repo / 'src/scworkbench_registry_json.R'
subprocess.run(['patch', '--batch', '--forward', '-p1', '-i', str(patch)],
               cwd=target, check=True)
shutil.copy2(helper, target / 'R/manifest_registry_json.R')
row_patch = repo / 'processing/parse_cytokines/recovery/object_row_keys.patch'
row_helper = repo / 'src/scworkbench_object_row_keys.R'
subprocess.run(['patch', '--batch', '--forward', '-p1', '-i', str(row_patch)],
               cwd=target, check=True)
shutil.copy2(row_helper, target / 'R/object_row_keys_character.R')
readonly_patch = repo / 'processing/parse_cytokines/recovery/feature_set_readonly.patch'
registry_materialization_source_sha256 = sha256(source / 'R/registry_materialization.R')
subprocess.run(['patch', '--batch', '--forward', '-p1', '-i', str(readonly_patch)],
               cwd=target, check=True)
artifact_patch = repo / 'processing/parse_cytokines/recovery/artifact_listing.patch'
artifact_source_sha256 = sha256(source / 'R/artifacts.R')
if artifact_source_sha256 != 'fcbc64ce02c3dac1ba98dc822c518c1781c71cda3ad33080f5f09a3ef379712b':
    raise SystemExit('Artifact listing source differs from the inspected build')
subprocess.run(['patch', '--batch', '--forward', '-p1', '-i', str(artifact_patch)],
               cwd=target, check=True)
description = target / 'DESCRIPTION'
lines = description.read_text().splitlines(keepends=True)
if lines.count('Version: 0.1.0\n') != 1:
    raise SystemExit('Unexpected source package version')
description.write_text(''.join('Version: 0.1.0.9004\n' if line == 'Version: 0.1.0\n'
                               else line for line in lines))
(destination / 'library').mkdir()
(destination / 'build_provenance.json').write_text(json.dumps({
    'base_source': '/opt/scworkbench (inside supplied SIF)',
    'base_bundle_init_sha256': expected,
    'base_provenance_json_sha256': row_source_sha256,
    'package_version': '0.1.0.9004',
    'patch_sha256': sha256(patch),
    'helper_sha256': sha256(helper),
    'row_key_patch_sha256': sha256(row_patch),
    'row_key_helper_sha256': sha256(row_helper),
    'readonly_patch_sha256': sha256(readonly_patch),
    'artifact_listing_patch_sha256': sha256(artifact_patch),
    'base_artifacts_sha256': artifact_source_sha256,
    'base_registry_materialization_sha256': registry_materialization_source_sha256,
    'patched_bundle_init_sha256': sha256(target / 'R/bundle_init.R'),
    'identity_scheme_changed': False,
    'bundle_payload_migration_required': False,
}, indent=2) + '\n')
print(target)
