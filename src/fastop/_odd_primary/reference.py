"""Temporary oddp-backed reference implementation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastop._odd_primary.evaluate import evaluate_all_targets
from fastop._odd_primary.indices import OperationIndex
from fastop._odd_primary.universal import UniversalOperation
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
    """Build universal data and evaluate it natively."""
    return cochain_operation_vector_from_universal(
        complex_,
        cochain,
        universal_operation(index),
        target_face_to_index,
    )


def cochain_operation_vector_oddp(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    index: OperationIndex,
    target_face_to_index: dict["Simplex", int],
    *,
    algorithm: str,
) -> Vector:
    """Apply oddp's cochain evaluator directly for comparison tests."""
    Steenrod = _load_steenrod()

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


def universal_operation(index: OperationIndex) -> UniversalOperation:
    """Build universal tensor data using the current oddp implementation."""
    Steenrod = _load_steenrod()
    tensor_chain = Steenrod.chain_operations(
        index.p,
        index.oddp_s,
        index.oddp_q,
        bockstein=index.bockstein,
    )
    return UniversalOperation.from_terms(index, dict(tensor_chain))


def cochain_operation_vector_from_universal(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
    target_face_to_index: dict["Simplex", int],
) -> Vector:
    """Evaluate universal data natively and return a target-degree vector."""
    result = evaluate_all_targets(
        complex_.faces(universal.target_degree),
        cochain,
        universal,
    )
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
