# User narrow47 native margin-0.01 reduction

This diagnostic threshold-sensitivity step reuses the pinned calibrated signal
from workflow 23 and runs the native scWorkbench `cell_vector_family_reduce`
operation with minimum score 0.80 and minimum margin 0.01. It does not rerun
UCell, use lineage gates, promote score artifacts, or mutate the bundle.
