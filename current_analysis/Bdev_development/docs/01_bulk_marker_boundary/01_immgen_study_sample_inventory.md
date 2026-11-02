# ImmGen cis-regulatory atlas sample inventory

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** derive ImmGen lineage gene sets outside scWorkbench  
**Status:** study cohort, GENCODE vM25 matrix, and RNA library manifest frozen;
contrast design not yet frozen  
**Reviewed:** 2026-10-31

## Scope

The ImmGen input for this Bdev step is based on Yoshida et al., *The
cis-Regulatory Atlas of the Mouse Immune System* (Cell, 2019; DOI
10.1016/j.cell.2018.12.036). The paper's Table S1 is the authoritative cohort
definition. Bulk RNA analysis remains outside scWorkbench; only versioned
generated gene sets will be imported into scWorkbench for use on single-cell
data.

Primary sources:

- article: <https://doi.org/10.1016/j.cell.2018.12.036>;
- Table S1: <https://ars.els-cdn.com/content/image/1-s2.0-S0092867418316507-mmc1.xlsx>;
- Table S2 expression matrix: <https://ars.els-cdn.com/content/image/1-s2.0-S0092867418316507-mmc2.csv>;
- ATAC-seq: GEO `GSE100738`;
- corresponding low-input RNA-seq: GEO `GSE109125`, under umbrella
  `GSE127267`.

The downloaded Table S1 workbook had SHA-256
`948f8dd40f06e76db322cfa30d86303e18eb35743ce7a9a9350196ac77cb7947`.
It contains 181 ATAC libraries across 90 population labels and 169 RNA
libraries across 86 population labels.

The published Table S2 matrix had SHA-256
`1941954b3b18652624d375ecf14ed65c543f4f7b4be009134f87abb7aa67b335`
and exactly 86 population columns. Those columns are the matched subset below;
the four ATAC-only Table S1 populations are absent. This resolves the apparent
difference between the 90 raw ATAC population labels in Table S1 and the 86
cell populations described by the paper's joint analysis.

## Matching rule

A population is called RNA/ATAC matched when the Table S1 population identity
is present in both the `SortedPopulations` ATAC sheet and the `Read Statistics`
RNA sheet. This is a **population-level assay match**. The study states that
the corresponding populations were isolated in parallel, but Table S1 does
not provide a one-to-one pairing key between individual RNA and ATAC
libraries. Replicate numbers must therefore not be treated as paired across
assays.

One documented naming alias was reconciled:

- ATAC `FRC.CD140a+.Madcam-.CD35-.SLN` = RNA `FRC.SLN`.

There are also source-level naming variations that do not change the Table S1
population match. For example, Table S1 uses `T8.IEL.LCMV.d7.SI`, whereas the
RNA GEO titles and Table S2 use `T8.IEL.LCMV.d7.Gut`. GEO RNA replicate
numbers and some shortened population labels are not uniformly identical to
the Table S1 `sample.name` values. The final input manifest must therefore
resolve GSM accessions with Table S1 metadata and library statistics rather
than rely on direct string equality alone.

The full RNA cohort maps to an ATAC population after applying this alias. Four
ATAC populations have no RNA counterpart in Table S1 or `GSE109125`:

- `preB.FrD.BM`;
- `DC.103+11b+.SI`;
- `DC.103+11b-.SI`;
- `MF.LP.SI`.

Thus, 86 of 90 ATAC populations have a population-matched RNA assay. The four
additional Rorc thymocyte GSM records currently present in `GSE100738` are not
rows in the paper's Table S1 and are excluded from this study cohort.

## Lineage summary

The lineage labels below are taken directly from Table S1 rather than inferred
from sample names.

| Table S1 lineage | ATAC populations | ATAC libraries | RNA populations | RNA libraries | Matched populations | ATAC-only populations |
|---|---:|---:|---:|---:|---:|---|
| Stem&Prog | 5 | 5 | 5 | 10 | 5 | none |
| B | 19 | 39 | 18 | 33 | 18 | `preB.FrD.BM` |
| abT | 15 | 30 | 15 | 34 | 15 | none |
| T.act | 10 | 20 | 10 | 19 | 10 | none |
| gdT | 7 | 14 | 7 | 13 | 7 | none |
| ILC | 10 | 20 | 10 | 17 | 10 | none |
| monocyte | 2 | 6 | 2 | 4 | 2 | none |
| myeloid | 12 | 27 | 11 | 23 | 11 | `MF.LP.SI` |
| DC | 5 | 10 | 3 | 6 | 3 | `DC.103+11b+.SI`, `DC.103+11b-.SI` |
| stroma | 5 | 10 | 5 | 10 | 5 | none |
| **Total** | **90** | **181** | **86** | **169** | **86** | **4** |

## Populations by lineage

### Stem and progenitor

Matched: `LTHSC.34-.BM`, `LTHSC.34+.BM`, `STHSC.150-.BM`, `MPP4.135+.BM`,
and `MPP3.48+.BM`.

### B lineage

Matched: `proB.CLP.BM`, `proB.FrA.BM`, `proB.FrBC.BM`, `B.FrE.BM`, `B1b.PC`,
`B.T1.Sp`, `B.T2.Sp`, `B.T3.Sp`, `B.MZ.Sp`, `B.Fo.Sp`, `B.GC.CB.Sp`,
`B.GC.CC.Sp`, `B.PB.Sp`, `B.PC.Sp`, `B.PC.BM`, `B.mem.Sp`, `B.Sp`, and
`B.Fem.Sp`.

ATAC only: `preB.FrD.BM`.

### Alpha-beta T lineage

Matched: `preT.DN1.Th`, `preT.DN2a.Th`, `preT.DN2b.Th`, `preT.DN3.Th`,
`T.DN4.Th`, `T.ISP.Th`, `T.DP.Th`, `T.4.Th`, `T.8.Th`, `T.4.Nve.Sp`,
`T.4.Nve.Fem.Sp`, `Treg.4.FP3+.Nrplo.Co`, `Treg.4.25hi.Sp`, `T.8.Nve.Sp`, and
`NKT.Sp`.

### Activated T/NKT lineage

Matched: `T.4.Sp.aCD3+CD40.18hr`, `T8.TN.P14.Sp`, `T8.IEL.LCMV.d7.SI`,
`T8.TE.LCMV.d7.Sp`, `T8.MP.LCMV.d7.Sp`, `T8.Tcm.LCMV.d180.Sp`,
`T8.Tem.LCMV.d180.Sp`, `NKT.Sp.LPS.3hr`, `NKT.Sp.LPS.18hr`, and
`NKT.Sp.LPS.3d`.

### Gamma-delta T lineage

Matched: `Tgd.g2+d17.24a+.Th`, `Tgd.g2+d17.LN`, `Tgd.g2+d1.24a+.Th`,
`Tgd.g2+d1.LN`, `Tgd.g1.1+d1.24a+.Th`, `Tgd.g1.1+d1.LN`, and `Tgd.Sp`.

### Innate lymphoid lineage

Matched: `NK.27+11b-.BM`, `NK.27+11b+.BM`, `NK.27-11b+.BM`,
`NK.27+11b-.Sp`, `NK.27+11b+.Sp`, `NK.27-11b+.Sp`, `ILC2.SI`,
`ILC3.NKp46-CCR6-.SI`, `ILC3.NKp46+.SI`, and `ILC3.CCR6+.SI`.

### Monocyte lineage

Matched: `Mo.6C+II-.Bl` and `Mo.6C-II-.Bl`.

### Other myeloid populations

Matched: `GN.BM`, `GN.Sp`, `GN.Thio.PC`, `MF.226+II+480lo.PC`,
`MF.ICAM+480hi.PC`, `MF.RP.Sp`, `MF.pIC.Alv.Lu`, `MF.Alv.Lu`,
`MF.microglia.CNS`, `MF.PC`, and `MF.Fem.PC`.

ATAC only: `MF.LP.SI`.

### Dendritic-cell lineage

Matched: `DC.4+.Sp`, `DC.8+.Sp`, and `DC.pDC.Sp`.

ATAC only: `DC.103+11b+.SI` and `DC.103+11b-.SI`.

### Stromal lineage

Matched: `FRC.CD140a+.Madcam-.CD35-.SLN` (RNA alias `FRC.SLN`), `IAP.SLN`,
`BEC.SLN`, `LEC.SLN`, and `Ep.MEChi.Th`.

## Replicate and design cautions

Fifteen RNA populations have only one library in the paper's Table S1:

- B: `B.Fem.Sp`, `B.PC.BM`, `B.PC.Sp`, `B.Sp`, `B.T1.Sp`, `B.T3.Sp`, and
  `B.mem.Sp`;
- ILC: `ILC3.CCR6+.SI`, `ILC3.NKp46-CCR6-.SI`, and `NK.27-11b+.Sp`;
- activated T/NKT: `NKT.Sp.LPS.3d` and `T8.Tem.LCMV.d180.Sp`;
- alpha-beta T: `NKT.Sp`;
- gamma-delta T: `Tgd.g2+d17.24a+.Th`;
- myeloid: `MF.pIC.Alv.Lu`.

Different populations within one lineage are not interchangeable biological
replicates. Lineage, tissue, activation state, contributing laboratory, sex,
and fine cell identity overlap in this atlas. Marker discovery must retain the
individual RNA libraries and an explicit sample-design table; it must not run
on population means or treat subtype diversity as replicate sampling without a
reviewed model.

## Consequence for the Bdev workflow

The 86 RNA populations define the available bulk-expression cohort. The four
ATAC-only populations cannot contribute to RNA-derived gene sets. The matched
subset should be retained in the metadata because it will also define which
ImmGen populations can later be evaluated in both expression and chromatin
space.

The exact GENCODE vM25 RNA cohort, its library-to-Table S1 mappings, and the
selected allowlist are frozen in the capability inputs. The statistical
contract and its annotation, state-gene, and software preflight were
subsequently frozen and validated; marker-model execution remains unrun.

## Approved steady-state discovery cohort

**Status:** cohort boundary, statistical method, selection thresholds, and
analysis preflight approved or frozen on 2026-10-31.

Use the 64 populations below as the primary marker-discovery cohort. This
cohort retains untreated, non-reporter reference populations, excludes
explicit sex-specific samples from model fitting, and contains only
populations with both RNA and ATAC data in the study. Paper Table S1 lists 127
RNA libraries for these populations; the selected GENCODE vM25 matrix contains
126 library columns after applying the population allowlist. The vM25 library
manifest, rather than the paper replicate count, will define the executable
cohort. The retained populations are:

- Stem and progenitor (5 populations, 10 vM25 libraries): `LTHSC.34-.BM`,
  `LTHSC.34+.BM`, `STHSC.150-.BM`, `MPP4.135+.BM`, and `MPP3.48+.BM`;
- B lineage (11 populations, 21 vM25 libraries): `proB.CLP.BM`, `proB.FrA.BM`,
  `proB.FrBC.BM`, `B.FrE.BM`, `B1b.PC`, `B.T1.Sp`, `B.T2.Sp`, `B.T3.Sp`,
  `B.MZ.Sp`, `B.Fo.Sp`, and `B.Sp`;
- alpha-beta T lineage (13 populations, 24 vM25 libraries): `preT.DN1.Th`,
  `preT.DN2a.Th`, `preT.DN2b.Th`, `preT.DN3.Th`, `T.DN4.Th`, `T.ISP.Th`,
  `T.DP.Th`, `T.4.Th`, `T.8.Th`, `T.4.Nve.Sp`, `T.8.Nve.Sp`,
  `Treg.4.25hi.Sp`, and `NKT.Sp`;
- gamma-delta T lineage (7 populations, 13 vM25 libraries): all seven Table S1
  `gdT` populations;
- innate lymphoid lineage (10 populations, 18 vM25 libraries): all ten Table S1
  `ILC` populations;
- monocyte lineage (2 populations, 4 vM25 libraries): both Table S1 monocyte
  populations;
- other myeloid populations (8 populations, 16 vM25 libraries): `GN.BM`, `GN.Sp`,
  `MF.226+II+480lo.PC`, `MF.ICAM+480hi.PC`, `MF.RP.Sp`, `MF.Alv.Lu`,
  `MF.microglia.CNS`, and `MF.PC`;
- dendritic-cell lineage (3 populations, 8 vM25 libraries): all three RNA-profiled
  Table S1 DC populations; and
- stromal lineage (5 populations, 12 vM25 libraries): all five Table S1 stromal
  populations.

Do not use the `T.act` group for steady-state marker discovery.

### Proposed exclusions

Exclude 14 explicitly perturbed populations (27 RNA libraries):

- NP-KLH/alum day 14: `B.GC.CB.Sp`, `B.GC.CC.Sp`, and `B.mem.Sp`;
- thioglycollate-induced: `GN.Thio.PC`;
- poly(I:C)-treated: `MF.pIC.Alv.Lu`;
- anti-CD3 plus anti-CD40: `T.4.Sp.aCD3+CD40.18hr`;
- LCMV infection: `T8.IEL.LCMV.d7.SI`, `T8.TE.LCMV.d7.Sp`,
  `T8.MP.LCMV.d7.Sp`, `T8.Tcm.LCMV.d180.Sp`, and
  `T8.Tem.LCMV.d180.Sp`; and
- LPS treatment: `NKT.Sp.LPS.3hr`, `NKT.Sp.LPS.18hr`, and `NKT.Sp.LPS.3d`.

Exclude five engineered reporter or transgenic populations (9 RNA libraries)
from primary discovery:

- Prdm1-GFP: `B.PB.Sp`, `B.PC.Sp`, and `B.PC.BM`;
- Foxp3-IRES-GFP: `Treg.4.FP3+.Nrplo.Co`; and
- P14 TCR transgenic: `T8.TN.P14.Sp`.

Retain these populations only as secondary checks of whether a steady-state
lineage gene set generalizes to terminal differentiation or engineered
backgrounds.

Exclude three explicitly female populations (6 RNA libraries) from discovery:

- `B.Fem.Sp`;
- `T.4.Nve.Fem.Sp`; and
- `MF.Fem.PC`.

Use them as held-out comparisons with `B.Sp`, `T.4.Nve.Sp`, and `MF.PC`,
respectively. Candidate lineage genes should be rejected if their apparent
specificity is driven by these male/female contrasts. Y-linked genes and
`Xist` should be excluded a priori. Other genes showing a strong sex effect in
the three available matched comparisons should also be removed. This prevents
known sex signals from entering the gene sets, but it does not make all
lineages fully sex-balanced; that limitation must be reported.

### Proposed marker-selection principle

Fine populations are biological strata, not exchangeable replicates. Derive
each lineage gene set from genes that are consistently elevated across
multiple retained populations within that lineage and against approved
competitor lineages. Preserve library-level counts, account for the fine
population structure, and require leave-one-population-out stability. Use the
excluded sex, treatment, and engineered populations only for sensitivity and
generalization checks, never as discovery replicates.
