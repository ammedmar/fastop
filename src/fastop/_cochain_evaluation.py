"""Cochain-level Steenrod operation evaluators."""

from __future__ import annotations

from itertools import combinations, product
from typing import TYPE_CHECKING

from fastop._universal import SignatureTable, UniversalOperation

if TYPE_CHECKING:
    from fastop.simplicial import Simplex


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
    answer: dict["Simplex", int] = {}
    for target in target_faces:
        coefficient = 0
        for tensor, tensor_coefficient in universal.terms.items():
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
