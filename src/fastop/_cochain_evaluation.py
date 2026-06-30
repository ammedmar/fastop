"""Cochain-level Steenrod operation evaluators."""

from __future__ import annotations

from itertools import combinations, product
from typing import TYPE_CHECKING

from fastop._linear_algebra import Vector
from fastop._universal import SignatureTable, UniversalOperation, universal_operation

try:
    from fastop._native import evaluate_all_targets as _native_evaluate_all_targets
except ImportError:  # pragma: no cover - depends on optional extension build
    _native_evaluate_all_targets = None

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
    """Build universal data and evaluate it on a cochain."""
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


def cochain_operation_vector_from_universal(
    complex_: "SimplicialComplex",
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
    target_face_to_index: dict["Simplex", int],
    *,
    algorithm: str = "auto",
) -> Vector:
    """Evaluate universal data and return a target-degree vector."""
    target_faces = complex_.faces(universal.target_degree)
    if algorithm == "auto":
        algorithm = _auto_evaluation_algorithm(target_faces, cochain, universal)
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
        raise ValueError(f"unknown cochain evaluation algorithm: {algorithm!r}")
    return {
        target_face_to_index[simplex]: coefficient
        for simplex, coefficient in result.items()
    }


def _auto_evaluation_algorithm(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
) -> str:
    """Choose an evaluator using coarse work estimates."""
    support_size = sum(
        1
        for face, coefficient in cochain.items()
        if coefficient % universal.p and len(face) == universal.source_degree + 1
    )
    if support_size == 0:
        return "all_targets"

    source_work = support_size ** universal.p
    target_work = len(target_faces) * max(len(universal.terms), 1)

    if source_work < target_work // 4:
        return "source_focused"
    return "all_targets"


def evaluate_source_mod_2(
    target_length: int,
    cocycle_support: list["Simplex"] | tuple["Simplex", ...],
    target_simplices: set["Simplex"] | frozenset["Simplex"],
) -> set["Simplex"]:
    """Evaluate the current source-focused mod-2 square support rule."""
    answer: set["Simplex"] = set()
    for left, right in combinations(cocycle_support, 2):
        left_vertices = set(left)
        right_vertices = set(right)
        union = left_vertices | right_vertices
        if len(union) != target_length:
            continue

        simplex = tuple(sorted(union))
        if simplex not in target_simplices:
            continue

        left_only = left_vertices - right_vertices
        right_only = right_vertices - left_vertices
        symmetric_difference = sorted(left_only | right_only)
        indices = {}
        for vertex in symmetric_difference:
            indices[vertex] = (
                simplex.index(vertex) + symmetric_difference.index(vertex)
            ) % 2
        left_indices = {indices[vertex] for vertex in left_only}
        right_indices = {indices[vertex] for vertex in right_only}
        if (left_indices == {0} and right_indices == {1}) or (
            left_indices == {1} and right_indices == {0}
        ):
            if simplex in answer:
                answer.remove(simplex)
            else:
                answer.add(simplex)
    return answer


def evaluate_all_targets(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    universal: UniversalOperation,
) -> dict["Simplex", int]:
    """Evaluate every target simplex against every universal tensor term."""
    if _native_evaluate_all_targets is not None:
        return _native_evaluate_all_targets(
            target_faces,
            cochain,
            universal.p,
            universal.terms,
        )

    tensor_terms = tuple(universal.terms.items())
    factors = {
        factor
        for tensor, _ in tensor_terms
        for factor in tensor
    }
    factor_reference_count = sum(len(tensor) for tensor, _ in tensor_terms)
    if len(factors) < factor_reference_count:
        return _evaluate_all_targets_with_factor_cache(
            target_faces,
            cochain,
            universal.p,
            tensor_terms,
        )

    answer: dict["Simplex", int] = {}
    for target in target_faces:
        coefficient = 0
        for tensor, tensor_coefficient in tensor_terms:
            term_value = tensor_coefficient
            for factor in tensor:
                source = tuple(target[index] for index in factor)
                source_coefficient = cochain.get(source)
                if source_coefficient is None:
                    break
                term_value *= source_coefficient
            else:
                coefficient += term_value
        coefficient %= universal.p
        if coefficient:
            answer[target] = coefficient
    return answer


def _evaluate_all_targets_with_factor_cache(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    p: int,
    tensor_terms: tuple[tuple[tuple[tuple[int, ...], ...], int], ...],
) -> dict["Simplex", int]:
    answer: dict["Simplex", int] = {}
    missing = object()
    for target in target_faces:
        coefficient = 0
        factor_coefficients = {}
        for tensor, tensor_coefficient in tensor_terms:
            term_value = tensor_coefficient
            for factor in tensor:
                source_coefficient = factor_coefficients.get(factor, missing)
                if source_coefficient is missing:
                    source = tuple(target[index] for index in factor)
                    source_coefficient = cochain.get(source)
                    factor_coefficients[factor] = source_coefficient
                if source_coefficient is None:
                    break
                term_value *= source_coefficient
            else:
                coefficient += term_value
        coefficient %= p
        if coefficient:
            answer[target] = coefficient
    return answer


def evaluate_source_focused(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    signatures: SignatureTable,
) -> dict["Simplex", int]:
    """Evaluate by enumerating ordered tuples from the cochain support."""
    target_face_set = set(target_faces)
    support = tuple(
        (face, coefficient % signatures.p)
        for face, coefficient in sorted(cochain.items())
        if coefficient % signatures.p and len(face) == signatures.source_degree + 1
    )
    answer: dict["Simplex", int] = {}
    target_length = signatures.target_degree + 1

    for source_tuple in product(support, repeat=signatures.p):
        target = tuple(
            sorted({vertex for source, _ in source_tuple for vertex in source})
        )
        if len(target) != target_length or target not in target_face_set:
            continue

        pattern = tuple(
            _omitted_positions_in_simplex(target, source)
            for source, _ in source_tuple
        )
        if any(
            len(omitted) != signatures.missing_vertices_per_factor
            for omitted in pattern
        ):
            continue

        coefficient = signatures.coefficients.get(pattern)
        if coefficient is None:
            continue

        term_value = coefficient
        for _, source_coefficient in source_tuple:
            term_value *= source_coefficient
        term_value %= signatures.p
        if term_value:
            answer[target] = (answer.get(target, 0) + term_value) % signatures.p
            if answer[target] == 0:
                del answer[target]

    return answer


def evaluate_target_omissions(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    signatures: SignatureTable,
) -> dict["Simplex", int]:
    """Evaluate every target simplex using omission-pattern signatures."""
    answer: dict["Simplex", int] = {}
    for target in target_faces:
        coefficient = 0
        for pattern, pattern_coefficient in signatures.coefficients.items():
            term_value = pattern_coefficient
            for omitted_positions in pattern:
                source = _selected_face_from_omissions(target, omitted_positions)
                source_coefficient = cochain.get(source)
                if source_coefficient is None:
                    break
                term_value *= source_coefficient
            else:
                coefficient += term_value
        coefficient %= signatures.p
        if coefficient:
            answer[target] = coefficient
    return answer


def _omitted_positions_in_simplex(
    target: "Simplex",
    source: "Simplex",
) -> tuple[int, ...]:
    source_vertices = set(source)
    return tuple(
        index
        for index, vertex in enumerate(target)
        if vertex not in source_vertices
    )


def _selected_face_from_omissions(
    target: "Simplex",
    omitted_positions: tuple[int, ...],
) -> "Simplex":
    omitted = set(omitted_positions)
    return tuple(
        vertex
        for index, vertex in enumerate(target)
        if index not in omitted
    )
