"""Fast computations of cohomology operations on simplicial data."""

from fastop.cohomology import (
    Mod2Cohomology,
    Mod2CohomologyElement,
    PrimeFieldCohomology,
    PrimeFieldCohomologyElement,
)
from fastop.delta_complex import DeltaComplex
from fastop import spaces
from fastop.simplicial import SimplicialComplex

__version__ = "0.1.0.dev0"

__all__ = [
    "Mod2Cohomology",
    "Mod2CohomologyElement",
    "PrimeFieldCohomology",
    "PrimeFieldCohomologyElement",
    "DeltaComplex",
    "SimplicialComplex",
    "spaces",
    "__version__",
]
