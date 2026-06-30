"""Temporary oddp-backed reference implementation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._linear_algebra import Vector
from fastop._universal import TensorTerm

if TYPE_CHECKING:
    from fastop.simplicial import Simplex, SimplicialComplex


def cochain_operation_vector_oddp(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    *,
    p: int,
    bockstein: bool,
    oddp_s: int,
    oddp_q: int,
    target_face_to_index: dict["Simplex", int],
    algorithm: str,
) -> Vector:
    """Apply oddp's cochain evaluator directly for comparison tests."""
    Steenrod = _load_steenrod()

    result = Steenrod.cochain_operation(
        _complex_for_oddp(complex_),
        cochain,
        p,
        oddp_s,
        oddp_q,
        bockstein=bockstein,
        algorithm=algorithm,
    )
    return {
        target_face_to_index[simplex]: coefficient % p
        for simplex, coefficient in result.items()
        if coefficient % p
    }


def universal_terms_oddp(
    *,
    p: int,
    bockstein: bool,
    oddp_s: int,
    oddp_q: int,
) -> dict[TensorTerm, int]:
    """Return oddp's universal tensor-chain terms."""
    Steenrod = _load_steenrod()
    return dict(
        Steenrod.chain_operations(
            p,
            oddp_s,
            oddp_q,
            bockstein=bockstein,
        )
    )


def _load_steenrod():
    try:
        from oddp import Steenrod
    except ImportError as exc:
        raise ImportError(
            "oddp is required for odd-primary Steenrod operations; "
            "install oddp or add it to PYTHONPATH"
        ) from exc
    return Steenrod


def _complex_for_oddp(complex_: "SimplicialComplex") -> dict[int, set["Simplex"]]:
    return {
        degree: set(complex_.faces(degree))
        for degree in range(complex_.dimension + 1)
    }
