"""Cochain-level Steenrod operation evaluators."""

from __future__ import annotations

from itertools import combinations, product
from typing import TYPE_CHECKING

from fastop._linear_algebra import Vector
from fastop._universal import (
    OmissionPattern,
    SignatureTable,
    UniversalOperation,
    universal_operation,
)

try:
    from fastop._native import evaluate_all_targets as _native_evaluate_all_targets
    from fastop._native import (
        evaluate_source_mod_3_covered as _native_evaluate_source_mod_3_covered,
    )
except ImportError:  # pragma: no cover - depends on optional extension build
    _native_evaluate_all_targets = None
    _native_evaluate_source_mod_3_covered = None

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
        result = evaluate_source_mod_3(
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
    if universal.p == 3:
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
    """Evaluate the source-focused mod-2 square support rule."""
    support_by_length: dict[int, list[tuple["Simplex", set[int]]]] = {}
    for face in cocycle_support:
        support_by_length.setdefault(len(face), []).append((face, set(face)))

    answer: set["Simplex"] = set()
    for source_length, support in support_by_length.items():
        overlap_length = 2 * source_length - target_length
        if overlap_length < 0 or overlap_length > source_length:
            continue

        buckets: dict["Simplex", list[tuple["Simplex", set[int]]]] = {}
        for face, vertices in support:
            for overlap in combinations(face, overlap_length):
                buckets.setdefault(overlap, []).append((face, vertices))

        for candidates in buckets.values():
            for left, right in combinations(candidates, 2):
                _, left_vertices = left
                _, right_vertices = right
                union = left_vertices | right_vertices
                if len(union) != target_length:
                    continue

                simplex = tuple(sorted(union))
                if simplex not in target_simplices:
                    continue

                if _mod2_square_pair_contributes(
                    simplex,
                    left_vertices,
                    right_vertices,
                ):
                    if simplex in answer:
                        answer.remove(simplex)
                    else:
                        answer.add(simplex)
    return answer


def _evaluate_source_mod_2_bruteforce(
    target_length: int,
    cocycle_support: list["Simplex"] | tuple["Simplex", ...],
    target_simplices: set["Simplex"] | frozenset["Simplex"],
) -> set["Simplex"]:
    """Evaluate the mod-2 square support rule by checking all source pairs."""
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

        if _mod2_square_pair_contributes(simplex, left_vertices, right_vertices):
            if simplex in answer:
                answer.remove(simplex)
            else:
                answer.add(simplex)
    return answer


def _mod2_square_pair_contributes(
    simplex: "Simplex",
    left_vertices: set[int],
    right_vertices: set[int],
) -> bool:
    """Return whether one source pair contributes to one mod-2 target."""
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
    return (left_indices == {0} and right_indices == {1}) or (
        left_indices == {1} and right_indices == {0}
    )


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


def evaluate_source_mod_3(
    target_faces: set["Simplex"] | frozenset["Simplex"],
    cochain: dict["Simplex", int],
    signatures: SignatureTable,
) -> dict["Simplex", int]:
    """Evaluate p=3 operations using pair-determined omission combinatorics."""
    if signatures.p != 3:
        raise ValueError("source_mod_3 is only available at p=3")

    support = tuple(
        (face, coefficient % 3)
        for face, coefficient in sorted(cochain.items())
        if coefficient % 3 and len(face) == signatures.source_degree + 1
    )
    if not support:
        return {}

    covered_patterns = {}
    uncovered_patterns = {}
    for pattern, coefficient in signatures.coefficients.items():
        if _pattern_covers_target(pattern, signatures.target_degree):
            covered_patterns[pattern] = coefficient
        else:
            uncovered_patterns[pattern] = coefficient

    answer: dict["Simplex", int] = {}
    if covered_patterns:
        covered_signature_table = SignatureTable(
            p=signatures.p,
            r=signatures.r,
            source_degree=signatures.source_degree,
            bockstein=signatures.bockstein,
            target_degree=signatures.target_degree,
            missing_vertices_per_factor=signatures.missing_vertices_per_factor,
            coefficients=covered_patterns,
        )
        if _native_evaluate_source_mod_3_covered is not None:
            answer.update(
                _native_evaluate_source_mod_3_covered(
                    target_faces,
                    support,
                    covered_signature_table.coefficients,
                    covered_signature_table.target_degree,
                )
            )
        else:
            _add_source_mod_3_covered_patterns(
                answer,
                set(target_faces),
                support,
                covered_signature_table,
            )

    if uncovered_patterns:
        fallback = evaluate_target_omissions(
            target_faces,
            cochain,
            SignatureTable(
                p=signatures.p,
                r=signatures.r,
                source_degree=signatures.source_degree,
                bockstein=signatures.bockstein,
                target_degree=signatures.target_degree,
                missing_vertices_per_factor=signatures.missing_vertices_per_factor,
                coefficients=uncovered_patterns,
            ),
        )
        for target, coefficient in fallback.items():
            answer[target] = (answer.get(target, 0) + coefficient) % 3
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


def _add_source_mod_3_covered_patterns(
    answer: dict["Simplex", int],
    target_face_set: set["Simplex"],
    support: tuple[tuple["Simplex", int], ...],
    signatures: SignatureTable,
) -> None:
    target_length = signatures.target_degree + 1
    positions = tuple(range(target_length))
    support_indexes: dict[
        tuple[int, ...],
        dict[tuple[int, ...], list[tuple["Simplex", int]]],
    ] = {}

    for pattern, pattern_coefficient in signatures.coefficients.items():
        selected_positions = tuple(
            tuple(position for position in positions if position not in omissions)
            for omissions in pattern
        )
        anchor_a, anchor_b, remaining = _best_source_mod_3_anchor_pair(
            selected_positions
        )
        positions_a = selected_positions[anchor_a]
        positions_b = selected_positions[anchor_b]
        positions_c = selected_positions[remaining]
        pair_cover = set(positions_a).union(positions_b)
        fixed_positions_c = tuple(
            position for position in positions_c if position in pair_cover
        )
        fixed_indices_c = tuple(
            index
            for index, position in enumerate(positions_c)
            if position in pair_cover
        )
        source_c_by_key = support_indexes.setdefault(
            fixed_indices_c,
            _source_index_by_positions(support, fixed_indices_c),
        )

        for source_a, coefficient_a in support:
            for source_b, coefficient_b in support:
                partial_target = _partial_target_from_pair(
                    source_a,
                    source_b,
                    positions_a,
                    positions_b,
                    target_length,
                )
                if partial_target is None:
                    continue

                key = tuple(partial_target[position] for position in fixed_positions_c)
                for source_c, coefficient_c in source_c_by_key.get(key, ()):
                    target = _target_from_remaining_source(
                        partial_target,
                        source_c,
                        positions_c,
                    )
                    if target is None or target not in target_face_set:
                        continue

                    source_coefficients = [0, 0, 0]
                    source_coefficients[anchor_a] = coefficient_a
                    source_coefficients[anchor_b] = coefficient_b
                    source_coefficients[remaining] = coefficient_c
                    term_value = pattern_coefficient
                    for coefficient in source_coefficients:
                        term_value *= coefficient
                    term_value %= 3
                    if term_value:
                        answer[target] = (answer.get(target, 0) + term_value) % 3
                        if answer[target] == 0:
                            del answer[target]


def _pattern_covers_target(
    pattern: OmissionPattern,
    target_degree: int,
) -> bool:
    positions = set(range(target_degree + 1))
    covered = set()
    for omissions in pattern:
        covered.update(positions.difference(omissions))
    return covered == positions


def _best_source_mod_3_anchor_pair(
    selected_positions: tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]],
) -> tuple[int, int, int]:
    pairs = ((0, 1, 2), (0, 2, 1), (1, 2, 0))
    return max(
        pairs,
        key=lambda pair: len(
            set(selected_positions[pair[0]]).union(selected_positions[pair[1]])
        ),
    )


def _source_index_by_positions(
    support: tuple[tuple["Simplex", int], ...],
    positions: tuple[int, ...],
) -> dict[tuple[int, ...], list[tuple["Simplex", int]]]:
    index: dict[tuple[int, ...], list[tuple["Simplex", int]]] = {}
    for source, coefficient in support:
        key = tuple(source[position] for position in positions)
        index.setdefault(key, []).append((source, coefficient))
    return index


def _partial_target_from_pair(
    source_a: "Simplex",
    source_b: "Simplex",
    positions_a: tuple[int, ...],
    positions_b: tuple[int, ...],
    target_length: int,
) -> list[int | None] | None:
    target: list[int | None] = [None] * target_length
    for index, position in enumerate(positions_a):
        target[position] = source_a[index]
    for index, position in enumerate(positions_b):
        vertex = source_b[index]
        if target[position] is not None and target[position] != vertex:
            return None
        target[position] = vertex
    if not _assigned_positions_are_increasing(target):
        return None
    return target


def _target_from_remaining_source(
    partial_target: list[int | None],
    source: "Simplex",
    positions: tuple[int, ...],
) -> "Simplex" | None:
    target = partial_target.copy()
    for index, position in enumerate(positions):
        vertex = source[index]
        if target[position] is not None and target[position] != vertex:
            return None
        target[position] = vertex
    if any(vertex is None for vertex in target):
        return None
    if not all(target[index] < target[index + 1] for index in range(len(target) - 1)):
        return None
    return tuple(target)


def _assigned_positions_are_increasing(target: list[int | None]) -> bool:
    last = None
    for vertex in target:
        if vertex is None:
            continue
        if last is not None and last >= vertex:
            return False
        last = vertex
    return True


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
