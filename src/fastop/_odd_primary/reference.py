"""Temporary oddp-backed reference implementation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._odd_primary.indices import OperationIndex
from fastop._prime_field import Vector

if TYPE_CHECKING:
    from fastop.simplicial import Simplex, SimplicialComplex


def cochain_operation_vector(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    index: OperationIndex,
    target_face_to_index: dict["Simplex", int],
    *,
    algorithm: str,
) -> Vector:
    """Apply the current oddp implementation and return a target-degree vector."""
    try:
        from oddp import Steenrod
    except ImportError as exc:
        raise ImportError(
            "oddp is required for odd-primary Steenrod operations; "
            "install oddp or add it to PYTHONPATH"
        ) from exc

    result = Steenrod.cochain_operation(
        _complex_for_oddp(complex_),
        cochain,
        index.p,
        index.oddp_s,
        index.oddp_q,
        bockstein=index.bockstein,
        algorithm=algorithm,
    )
    return {
        target_face_to_index[simplex]: coefficient % index.p
        for simplex, coefficient in result.items()
        if coefficient % index.p
    }


def _complex_for_oddp(complex_: "SimplicialComplex") -> dict[int, set["Simplex"]]:
    return {
        degree: set(complex_.faces(degree))
        for degree in range(complex_.dimension + 1)
    }
