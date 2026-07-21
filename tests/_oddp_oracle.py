"""Development-only adapters for comparing fastop with the oddp oracle."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._linear_algebra import Vector
from fastop._universal import UniversalOperation

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
    """Apply oddp's cochain evaluator and return a fastop sparse vector."""
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


def universal_operation_oddp(
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    missing_vertices_per_factor: int,
    oddp_s: int,
    oddp_q: int,
) -> UniversalOperation:
    """Build universal tensor data using oddp as an external oracle."""
    Steenrod = _load_steenrod()
    return UniversalOperation.from_terms(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
        terms=dict(
            Steenrod.chain_operations(
                p,
                oddp_s,
                oddp_q,
                bockstein=bockstein,
            )
        ),
    )


def _load_steenrod():
    try:
        from oddp import Steenrod
    except ImportError as exc:
        raise ImportError(
            "oddp is required only for the optional oracle tests"
        ) from exc
    return Steenrod


def _complex_for_oddp(complex_: "SimplicialComplex") -> dict[int, set["Simplex"]]:
    return {
        degree: set(complex_.faces(degree))
        for degree in range(complex_.dimension + 1)
    }
