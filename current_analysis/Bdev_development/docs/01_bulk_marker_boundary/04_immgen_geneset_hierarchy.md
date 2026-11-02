# ImmGen hierarchical gene-set and contrast design

**Capability:** `01_bulk_marker_boundary`  
**Bdev step:** define population roles and contrast neighborhoods  
**Status:** approved and validated  
**Approved:** 2026-10-31

## Decision

Use a hierarchical gene-set design with four distinct products:

1. superlineage sets for progenitor, lymphoid, myeloid/DC, and stromal
   separation;
2. broad Table S1 lineage sets;
3. B-lineage branch sets; and
4. local B developmental-stage sets defined against declared neighboring
   populations.

This hierarchy is frozen in
[`01_immgen_geneset_hierarchy.yaml`](../../config/01_bulk_marker_boundary/01_immgen_geneset_hierarchy.yaml).
It defines biological roles and contrasts only. It does not authorize or select
a differential-expression engine, design formula, or marker threshold.

## Broad-population roles

Do not remove broad populations from the frozen input cohort. Assign them
roles that prevent imprecise gates from becoming fine-stage targets:

| Population | Role |
|---|---|
| `B.Sp` | broad B-lineage validation only |
| `Tgd.Sp` | broad gamma-delta lineage validation only |
| `MF.PC` | broad macrophage validation only |

These populations do not contribute to discovery of broad-lineage or
fine-stage genes, but will test whether accepted sets generalize to less
specific gates.

Retain `NKT.Sp`, `GN.BM`, and `GN.Sp` in discovery because they preserve the
only approved steady-state NKT branch and the only neutrophil/granulocyte
branch. They are not treated as fine B-development targets.

After applying these role overrides, 61 of the 64 approved populations
contribute to broad discovery. All nine Table S1 lineages remain represented:

| Lineage | Discovery populations |
|---|---:|
| Stem&Prog | 5 |
| B | 10 |
| abT | 13 |
| gdT | 6 |
| ILC | 10 |
| monocyte | 2 |
| myeloid | 7 |
| DC | 3 |
| stroma | 5 |

Fine populations are equal-weight biological strata for broad-set discovery;
their library counts do not determine their biological weight.

## B-lineage hierarchy

The B discovery populations are grouped as:

| Branch | Populations |
|---|---|
| Early lymphoid progenitor | `proB.CLP.BM`, `proB.FrA.BM` |
| Committed bone-marrow B | `proB.FrBC.BM`, `B.FrE.BM` |
| Splenic transitional B | `B.T1.Sp`, `B.T2.Sp`, `B.T3.Sp` |
| Mature conventional B | `B.Fo.Sp`, `B.MZ.Sp` |
| B1b branch | `B1b.PC` |

The B1b branch can produce a candidate signature, but its single population
cannot establish cross-population branch stability.

The directed biological graph is:

```text
proB.CLP → proB.FrA → proB.FrBC ──gap──> B.FrE ──BM/spleen──> B.T1
                                                                  ├──> B.T2 ──> B.Fo
                                                                  │          └──> B.MZ
                                                                  └──> B.T3
```

The `proB.FrBC` to `B.FrE` edge spans the missing RNA-profiled `preB.FrD`
stage. The `B.FrE` to `B.T1` edge spans bone marrow and spleen. Neither edge
may be interpreted as a clean direct transition without explicit confounder
screening.

Local marker contrasts use a slightly wider symmetric neighbor graph. In
particular, `B.T2`, `B.T3`, follicular, and marginal-zone populations are also
compared with phenotypically close branch neighbors. Mouse studies support T1
as a precursor of T2 and T3 and the subsequent separation of follicular and
marginal-zone fates; the graph nevertheless treats disputed or indirect
relationships as contrast neighborhoods rather than proof of lineage ancestry.

## Validation result

The hierarchy validator passed and confirmed:

- all 64 selected populations are assigned consistently;
- the three broad-validation populations are excluded only from discovery,
  not removed from validation;
- all nine broad lineages remain represented by 61 discovery populations;
- unique NKT and granulocyte branches are retained;
- the five B branch groups partition the ten discovery-eligible B populations;
- nine fine B targets have symmetric local neighborhoods; and
- eight developmental edges and five branch contrasts have unique stable IDs.

## Next step

Freeze the statistical and marker-selection contract. It must define balanced
population weighting, count filtering and normalization, local and broad
models, covariate/confounding policy, effect and false-discovery requirements,
stability resampling, gene exclusions, gene-set sizes, and deterministic
ranking. No differential-expression analysis should run before that contract
is approved.

## Scientific sources

- Yoshida et al., *The cis-Regulatory Atlas of the Mouse Immune System*, Cell
  2019: <https://doi.org/10.1016/j.cell.2018.12.036>.
- Rowland et al., *BAFF-R signaling aids the differentiation of immature B
  cells into transitional B cells following tonic BCR signaling*, Journal of
  Immunology 2010: <https://pmc.ncbi.nlm.nih.gov/articles/PMC2950883/>.
- Kleiman et al., *Distinct Transcriptomic Features are Associated with
  Transitional and Mature B-Cell Populations in the Mouse Spleen*, Frontiers
  in Immunology 2015: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4324157/>.
