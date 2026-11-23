# Build and qualify the Bdev candidate SIF

This workflow extends the immutable `20261029` baseline with the exact clean
source commits pinned in config 01. It preserves the inherited PyBench runtime,
installs scWorkbench plus eight operation packs, and never overwrites an
existing image.

The candidate is not accepted merely because it builds. Qualification runs in
clean processes inside the candidate and requires exact installed versions and
manifest hashes, collision-free activation of all eight packs, the installed
core sparse typed-slice regression, all seven installed scwbATAC contracts,
working inherited Python scientific imports, and a recorded SIF SHA-256.

Run in order:

1. `01_preflight.py --config <config>`
2. `02_build_candidate_sif.sh`
3. `03_qualify_candidate_sif.py --config <config>`

Implementation remains locked until the qualified checksum is explicitly
accepted and pinned in the implementation runtime manifest.
