# ImmGen ATAC native-L4 trajectory mapping

## Outcome

Controlled development run `immgen_atac_l4_20261110T202521Z` completed the
real-ImmGen portion of Step 4. All 28 frozen GSE100738/SRP110978 paired-end
ATAC libraries representing 16 ImmGen B/Stem&Prog populations were reprocessed
against mm10 and counted on the exact ordered 191,433-feature Bdev peak axis.

The preprocessing receipt is complete. Every FASTQ passed its declared MD5,
every peak-count file contains the exact feature axis, and every library passed
the minimum retained-fragment gate. Retained deduplicated,
blacklist-filtered, Tn5-shifted fragments ranged from 2,323,295 to 9,542,315
per library. The declared archived mm10 Bowtie2 index and blacklist remained
read-only.

## Mapping contract

The mapping run used 1,382 Bdev reference cells and the qualified deterministic
10-dimensional Bdev-only PeakVI recipe. Its reference latent checksum was
reproduced exactly as
`3bb607f73151f8cc5ecc50c101da289eecfed50fdc43efd510c9db97b24eb818`.
The retained model projected all 28 ImmGen libraries without refitting, and
`model_bound_latent_union` authenticated the combined reference/query space.

Support calibration and transfer used native fitted Slingshot `L4`, primary
k 30, and the predeclared q95/q99 support thresholds. Categorical Bdev cluster
transfer was intentionally disabled.

| Support call | Libraries |
|---|---:|
| supported | 25 |
| borderline | 3 |
| out-of-reference | 0 |

The three borderline libraries are `GSM2692242` (`B.FrE.BM`), `GSM2692325`
(`LTHSC.34+.BM`), and `GSM2692316` (`MPP3.48+.BM`). All 28 libraries are
supported or borderline and therefore retain reportable native-L4 estimates;
the transfer contract would have withheld accepted L4 for any
out-of-reference library.

## Traceability and disposition

The source-bundle manifest remained
`cc8b92e8b55bb1b12b7f818c3bc272193c2139abec81b9777d48bfe2dd022a05`
before and after the run. The finalized mapping-bundle manifest is
`1fb8773e02cb7f8a9d85a27b37823381ba902599bc442de4bd657a85bd9d7934`.
The run inventory contains 48 checksummed result files and has SHA-256
`caacd244497cc3aed180cd9129683b725084dc1ba367acab69c1ef8a05c1a13f`.

This is a development qualification result, not a paper-analysis result. No
replacement SIF was built, `Bdev_implementation` remains locked, and these
outputs must not be promoted as implementation results. The mapping operation
also required a transient 2.0-GiB sparse-to-dense allocation during support
calibration; future runtime qualification should retain sufficient memory or
remove that allocation within the owning operation pack before the replacement
SIF is frozen.

Scientific acceptance was reopened on 2026-11-11 after visualization on the
approved `umap.atac` showed implausible placement, especially for later-stage
populations. The support labels and transferred L4 values remain reproducible
technical outputs, but they are not accepted biological results and must not be
used downstream. Step 4 remains unresolved pending investigation.
