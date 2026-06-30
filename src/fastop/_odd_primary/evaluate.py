"""Native evaluators for odd-primary universal operation data."""

from __future__ import annotations

from itertools import product
from typing import TYPE_CHECKING

from fastop._odd_primary.universal import SignatureTable, UniversalOperation

if TYPE_CHECKING:
    from fastop.simplicial import Simplex


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
            product = tensor_coefficient
            for factor in tensor:
                source = tuple(target[index] for index in factor)
                source_coefficient = cochain.get(source)
                if source_coefficient is None:
                    break
                product *= source_coefficient
            else:
                coefficient += product
        coefficient %= universal.p
        if coefficient:
            answer[target] = coefficient
    return answer


def evaluate_sparse_support(
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
