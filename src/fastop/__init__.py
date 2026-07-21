"""Fast computations of cohomology operations on simplicial data."""

from fastop.cohomology import (
    Mod2Cohomology,
    Mod2CohomologyElement,
    PrimeFieldCohomology,
    PrimeFieldCohomologyElement,
)
from fastop import spaces
from fastop.delta_complex import DeltaComplex
from fastop.group_action import FiniteGroupAction
from fastop.simplicial import SimplicialComplex
from fastop.simplicial_set import (
    SimplexReference,
    SimplicialSet,
    SymmetricPowerSimplicialSet,
)

__version__ = "0.1.0a1"

__all__ = [
    "Mod2Cohomology",
    "Mod2CohomologyElement",
    "PrimeFieldCohomology",
    "PrimeFieldCohomologyElement",
    "DeltaComplex",
    "FiniteGroupAction",
    "SimplicialComplex",
    "SimplexReference",
    "SimplicialSet",
    "SymmetricPowerSimplicialSet",
    "spaces",
    "__version__",
]
