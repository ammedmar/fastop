"""Temporary oddp-backed reference implementation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._cochain_evaluation import (
    evaluate_all_targets,
    evaluate_source_focused,
    evaluate_target_omissions,
)
from fastop._linear_algebra import Vector
from fastop._universal import UniversalOperation, native_universal_operation

if TYPE_CHECKING:
    from fastop.simplicial import Simplex, SimplicialComplex


def cochain_operation_vector(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    missing_vertices_per_factor: int,
    oddp_s: int,
    oddp_q: int,
    target_face_to_index: dict["Simplex", int],
    algorithm: str,
) -> Vector:
    """Build universal data and evaluate it natively."""
    return cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal_operation(
            p=p,
            r=r,
            source_degree=source_degree,
            bockstein=bockstein,
            target_degree=target_degree,
            missing_vertices_per_factor=missing_vertices_per_factor,
            oddp_s=oddp_s,
            oddp_q=oddp_q,
        ),
        target_face_to_index,
        algorithm=algorithm,
    )


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


def universal_operation(
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
    """Build universal tensor data, falling back to oddp when needed."""
    native = native_universal_operation(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
    )
    if native is not None:
        return native

    Steenrod = _load_steenrod()
    tensor_chain = Steenrod.chain_operations(
        p,
        oddp_s,
        oddp_q,
        bockstein=bockstein,
    )
    return UniversalOperation.from_terms(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
        terms=dict(tensor_chain),
    )


def cochain_operation_vector_from_universal(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
    target_face_to_index: dict["Simplex", int],
    *,
    algorithm: str = "source_focused",
) -> Vector:
    """Evaluate universal data natively and return a target-degree vector."""
    target_faces = complex_.faces(universal.target_degree)
    if algorithm in {"all_targets", "direct", "all-target"}:
        result = evaluate_all_targets(target_faces, cochain, universal)
    elif algorithm in {"target_omissions", "target", "signature", "signatures"}:
        result = evaluate_target_omissions(
            target_faces,
            cochain,
            universal.signature_table(),
        )
    elif algorithm in {"source_mod_3", "prime-three"}:
        if universal.p != 3:
            raise ValueError("source_mod_3 is only available at p=3")
        result = evaluate_source_focused(
            target_faces,
            cochain,
            universal.signature_table(),
        )
    elif algorithm in {"source_focused", "support", "sparse"}:
        result = evaluate_source_focused(
            target_faces,
            cochain,
            universal.signature_table(),
        )
    else:
        raise ValueError(f"unknown odd-primary evaluation algorithm: {algorithm!r}")
    return {
        target_face_to_index[simplex]: coefficient
        for simplex, coefficient in result.items()
    }


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
