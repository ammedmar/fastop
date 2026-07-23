"""Fast computations of cohomology operations on simplicial data."""

from fastop.cohomology import PrimeFieldCohomology, PrimeFieldCohomologyElement
from fastop import spaces
from fastop.delta_complex import DeltaComplex
from fastop.group_action import FiniteGroupAction
from fastop.simplicial import SimplicialComplex
from fastop.simplicial_set import (
    SimplexReference,
    SimplicialSet,
)

__version__ = "0.1.0a1"

__all__ = [
    "PrimeFieldCohomology",
    "PrimeFieldCohomologyElement",
    "DeltaComplex",
    "FiniteGroupAction",
    "SimplicialComplex",
    "SimplexReference",
    "SimplicialSet",
    "spaces",
    "__version__",
]
