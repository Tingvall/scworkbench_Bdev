#!/usr/bin/env python3
"""Copy a finalized bundle and retain byte-level preservation evidence."""
import hashlib
import json
from pathlib import Path
import shutil
import sys


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def inventory(root):
    rows = []
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise RuntimeError(f'Symlink requires explicit containment review: {path}')
        if not path.is_file():
            continue
        before = path.stat()
        digest = sha256(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError(f'Source changed during inventory: {path}')
        rows.append({'path': str(path.relative_to(root)), 'bytes': before.st_size,
                     'sha256': digest})
    return rows


def main():
    mode, source_arg, target_arg, report_arg = sys.argv[1:]
    source, target, report = map(Path, (source_arg, target_arg, report_arg))
    evidence = report / 'source_bundle_files.json'
    if mode == 'copy':
        if target.exists():
            raise SystemExit(f'Refusing to overwrite: {target}')
        lock = source.parent / f'.{source.name}.scworkbench-writer-lock'
        if lock.exists() or any((source / 'db').glob('*-wal')) or any((source / 'db').glob('*-shm')):
            raise SystemExit('Source has a writer lock or unresolved SQLite sidecars')
        print('Hashing original bundle', flush=True)
        original = inventory(source)
        evidence.write_text(json.dumps(original, indent=2) + '\n')
        print('Copying original bundle', flush=True)
        shutil.copytree(source, target, copy_function=shutil.copy2)
        print('Verifying every copied file', flush=True)
        if inventory(target) != original:
            raise RuntimeError('Copied bundle bytes differ from source inventory')
        (report / 'copy_validation.json').write_text(json.dumps({
            'source': str(source), 'destination': str(target), 'files': len(original),
            'bytes': sum(row['bytes'] for row in original), 'all_files_identical': True,
        }, indent=2) + '\n')
    elif mode == 'verify-source':
        unchanged = inventory(source) == json.loads(evidence.read_text())
        (report / 'source_preservation.json').write_text(json.dumps({
            'source': str(source), 'every_file_sha256_unchanged': unchanged,
        }, indent=2) + '\n')
        if not unchanged:
            raise RuntimeError('Original bundle changed during recovery')
        summary = report / 'REPORT.md'
        if summary.exists():
            with summary.open('a') as handle:
                handle.write('\nFinal preservation check passed: every original bundle file has its initial SHA-256.\n')
    else:
        raise SystemExit(f'Unknown mode: {mode}')
    print('PASS:', mode, flush=True)


if __name__ == '__main__':
    main()
