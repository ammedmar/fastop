"""Native construction of universal odd-primary operation formulas.

This module contains only the combinatorial formula-building path needed by
``fastop``.  It is derived from the MIT-licensed ``oddp`` project by Federico
Cantero-Morán and Anibal M. Medina-Mardones; see ``THIRD_PARTY_NOTICES.md``.
The implementation deliberately uses immutable faces and ordinary mappings
instead of importing ``oddp``'s resolution object hierarchy.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations, permutations
from math import factorial

Face = tuple[int, ...]
MilnorTensor = tuple[Face, ...]
TensorTerm = tuple[Face, ...]
Codegeneracy = dict[Face, dict[tuple[Face, Face], int]]


def computed_universal_terms(
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
) -> dict[TensorTerm, int]:
    """Compute the universal tensor formula for ``P^r`` or ``beta P^r``.

    The returned tensor factors are faces of the target simplex. Coefficients
    are reduced modulo ``p`` and zero terms are omitted.
    """
    if not _is_odd_prime(p):
        raise ValueError("p must be an odd prime")
    if not isinstance(r, int) or isinstance(r, bool) or r < 0:
        raise ValueError("r must be a nonnegative integer")
    if (
        not isinstance(source_degree, int)
        or isinstance(source_degree, bool)
        or source_degree < 0
    ):
        raise ValueError("source_degree must be a nonnegative integer")
    if not isinstance(bockstein, bool):
        raise TypeError("bockstein must be a boolean")

    missing = 2 * r * (p - 1) + int(bockstein)
    target_degree = source_degree + missing
    resolution_degree = missing * p

    periodic = _phi_dual(resolution_degree, p)
    milnor = _psi_dual_homogeneous(
        periodic,
        resolution_degree,
        p,
        target_degree,
    )
    tensor = _abc_dual(milnor, p, target_degree)
    return _alexander_dual(tensor, p, target_degree)


def _is_odd_prime(value: int) -> bool:
    if not isinstance(value, int) or isinstance(value, bool) or value < 3:
        return False
    return all(value % divisor for divisor in range(2, int(value**0.5) + 1))


def _add_mod(mapping: dict, key, coefficient: int, p: int) -> None:
    value = (mapping.get(key, 0) + coefficient) % p
    if value:
        mapping[key] = value
    else:
        mapping.pop(key, None)


def _sign_permutation(sequence) -> int:
    ordered = sorted(set(sequence))
    counts = {value: 0 for value in ordered}
    transpositions = 0
    for value in sequence:
        transpositions += sum(counts[other] for other in ordered if other > value)
        counts[value] += 1
    return (-1) ** (transpositions % 2)


def _sign_reorder(first, second) -> int:
    parity = sum(sum(other > value for other in first) for value in second)
    return (-1) ** (parity % 2)


def _sign_complement(face, subface) -> int:
    complement = tuple(sorted(set(face).difference(subface)))
    return _sign_reorder(complement, subface)


def _simplicial_action(face: Face, p: int, shift: int) -> tuple[Face, int]:
    """Apply the dual cyclic action to a face of the ``(p-1)``-simplex."""
    shift %= p
    nonwrapped = tuple(value - shift for value in face if value - shift >= 0)
    wrapped = tuple(value - shift + p for value in face if value - shift < 0)
    return nonwrapped + wrapped, (-1) ** (len(nonwrapped) * len(wrapped))


def _simplicial_alexander_dual(face: Face, dimension: int) -> tuple[Face, int]:
    vertices = range(dimension + 1)
    complement = tuple(value for value in vertices if value not in face)
    return complement, _sign_complement(vertices, face)


@lru_cache(maxsize=None)
def _partitions(face: Face) -> tuple[dict[tuple[Face, ...], int], ...]:
    """Return all ordered join decompositions of a nonempty face."""
    if len(face) == 1:
        return ({(face,): 1},)

    continuing: list[dict[tuple[Face, ...], int]] = [{(face,): 1}]
    finished: list[dict[tuple[Face, ...], int]] = [{}]
    for cardinality in range(1, len(face)):
        continuing.append({})
        finished.append({})
        for partition, coefficient in continuing[cardinality - 1].items():
            initial, end = partition[:-1], partition[-1]
            for first_size in range(1, len(end)):
                for first in combinations(end, first_size):
                    first_set = set(first)
                    second = tuple(value for value in end if value not in first_set)
                    new_partition = initial + (first, second)
                    new_coefficient = coefficient * _sign_reorder(first, second)
                    destination = (
                        continuing if first_size != len(end) - 1 else finished
                    )
                    destination[cardinality][new_partition] = new_coefficient
    for cardinality, complete in enumerate(finished):
        continuing[cardinality].update(complete)
    return tuple(continuing)


def _phi_dual(resolution_degree: int, p: int) -> dict[Face, int]:
    """Map the minimal-resolution generator into the periodic resolution."""
    if resolution_degree and resolution_degree % (p - 1) == 0:
        face_size = p - 1
    else:
        face_size = resolution_degree % (p - 1)

    result: dict[Face, int] = {}
    face = tuple(range(face_size))
    _add_mod(result, face, 1, p)
    last = tuple(
        range(
            max(p - 1 - (face_size % 2) - face_size, 0),
            max(p - 1 - (face_size % 2), 1),
        )
    )
    while face != last:
        for index, value in enumerate(reversed(face)):
            if value != last[-index - 1]:
                face = face[: -index - 1] + tuple(
                    range(value + 2, value + 3 + index)
                )
                _add_mod(result, face, 1, p)
                break

    if p != 3:
        quotient, remainder = divmod(resolution_degree, p - 1)
        half_factorial = factorial((p - 1) // 2)
        denominator = half_factorial ** (quotient + 1)
        numerator = factorial((p - 1 - remainder) // 2)
        scale = numerator * pow(denominator, -1, p) % p
        result = {
            face: coefficient * scale % p
            for face, coefficient in result.items()
            if coefficient * scale % p
        }
    return result


@lru_cache(maxsize=None)
def _straightening(p: int) -> dict[Face, int]:
    short: dict[Face, int] = {}
    full: dict[Face, int] = {}
    for size in range(1, p):
        for face in combinations(range(p), size):
            for shift in range(p):
                rotated, _ = _simplicial_action(face, p, shift)
                if rotated in short:
                    full[face] = (short[rotated] + shift) % p
                    break
            else:
                rotated, _ = _simplicial_action(face, p, face[0])
                for index in range(len(rotated)):
                    if rotated[index - 1] % p != rotated[index] - 1:
                        short[rotated] = rotated[index]
                        full[face] = (rotated[index] + face[0]) % p
                        break
    return full


@lru_cache(maxsize=None)
def _codegeneracy(p: int) -> Codegeneracy:
    straightening = _straightening(p)
    result: Codegeneracy = {}
    for omega_size in range(1, p):
        for tau_size in range(1, omega_size + 1):
            for omega in combinations(range(p), omega_size):
                for tau in combinations(omega, tau_size):
                    tau_set = set(tau)
                    remainder = tuple(value for value in omega if value not in tau_set)
                    first_sign = _sign_reorder(remainder, tau) * (-1) ** tau_size
                    for permuted_remainder in permutations(remainder):
                        second_sign = _sign_permutation(permuted_remainder)
                        subdivision_vertex = list(tau)
                        simplex_face = [straightening[tuple(subdivision_vertex)]]
                        for vertex in permuted_remainder:
                            subdivision_vertex.append(vertex)
                            simplex_vertex = straightening[
                                tuple(sorted(subdivision_vertex))
                            ]
                            if simplex_vertex in simplex_face:
                                break
                            simplex_face.append(simplex_vertex)
                        else:
                            third_sign = _sign_permutation(simplex_face)
                            source_face = tuple(sorted(simplex_face))
                            new_tau, fourth_sign = _simplicial_alexander_dual(
                                tau, p - 1
                            )
                            fifth_sign = (-1) ** omega_size
                            image = result.setdefault(source_face, {})
                            key = (omega, new_tau)
                            coefficient = (
                                first_sign
                                * second_sign
                                * third_sign
                                * fourth_sign
                                * fifth_sign
                            )
                            _add_mod(image, key, coefficient, p)
    return result


def _psi_dual_homogeneous(
    periodic: dict[Face, int],
    resolution_degree: int,
    p: int,
    target_degree: int,
) -> dict[MilnorTensor, int]:
    """Build only Milnor terms that can contribute homogeneous tensors."""
    result: dict[MilnorTensor, int] = {}
    if resolution_degree % p:
        return result

    weight = resolution_degree // p
    rule = tuple(-label % p for label in range(target_degree + 1))

    def emit(pairs, coefficient: int) -> None:
        slots: list[Face] = [()] * (target_degree + 1)
        for face, label in pairs:
            slots[label] = face
        _add_mod(result, tuple(slots), coefficient, p)

    if resolution_degree == 0:
        for coefficient in periodic.values():
            emit((), coefficient)
        return result

    rounds = (resolution_degree - 1) // (p - 1)

    def place(states, faces, base_mass: int, tail: int):
        mass = base_mass
        number_of_faces = len(faces)
        for index, face in enumerate(faces):
            mass += len(face)
            minimum_after = -(-(resolution_degree - mass) // (p - 1))
            after = max((number_of_faces - 1 - index) + tail, minimum_after)
            upper_label = target_degree + 1 - after
            new_states = []
            for pairs, last_label, weights in states:
                for label in range(last_label + 1, upper_label):
                    new_weights = list(weights)
                    for vertex in face:
                        residue = (vertex - rule[label]) % p
                        new_weights[residue] += 1
                        if new_weights[residue] > weight:
                            break
                    else:
                        new_states.append(
                            (
                                pairs + ((face, label),),
                                label,
                                tuple(new_weights),
                            )
                        )
            states = new_states
            if not states:
                return []
        return states

    codegeneracy = _codegeneracy(p)

    def recurse(
        states,
        active: Face,
        coefficient: int,
        remaining: int,
        base_mass: int,
    ):
        if remaining == 0:
            for pairs, last_label, weights in states:
                for label in range(last_label + 1, target_degree + 1):
                    new_weights = list(weights)
                    for vertex in active:
                        residue = (vertex - rule[label]) % p
                        new_weights[residue] += 1
                        if new_weights[residue] > weight:
                            break
                    else:
                        emit(pairs + ((active, label),), coefficient)
            return

        length = len(states[0][0]) + 1
        for (first, second), codegeneracy_coefficient in codegeneracy[
            active
        ].items():
            second_partitions = _partitions(second)
            maximum_size = min(
                len(second_partitions),
                target_degree + 1 - length - (remaining - 1),
            )
            for size in range(maximum_size):
                for partition, partition_coefficient in second_partitions[
                    size
                ].items():
                    batch = (first,) + partition[:-1]
                    new_states = place(states, batch, base_mass, remaining)
                    if new_states:
                        recurse(
                            new_states,
                            partition[-1],
                            coefficient
                            * codegeneracy_coefficient
                            * partition_coefficient
                            % p,
                            remaining - 1,
                            base_mass + sum(len(face) for face in batch),
                        )

    initial = [((), -1, (0,) * p)]
    for periodic_face, coefficient in periodic.items():
        periodic_partitions = _partitions(periodic_face)
        maximum_size = min(
            len(periodic_partitions), target_degree + 1 - rounds
        )
        for size in range(maximum_size):
            for partition, partition_coefficient in periodic_partitions[
                size
            ].items():
                batch = partition[:-1]
                states = place(initial, batch, 0, rounds + 1)
                if states:
                    recurse(
                        states,
                        partition[-1],
                        coefficient * partition_coefficient % p,
                        rounds,
                        sum(len(face) for face in batch),
                    )
    return result


def _abc_dual(
    milnor: dict[MilnorTensor, int], p: int, target_degree: int
) -> dict[TensorTerm, int]:
    result: dict[TensorTerm, int] = {}
    rule = tuple(-label % p for label in range(target_degree + 1))
    for tensor, coefficient in milnor.items():
        transformed: list[Face] = []
        action_sign = 1
        for label, factor in enumerate(tensor):
            new_factor, sign = _simplicial_action(factor, p, rule[label])
            transformed.append(new_factor)
            action_sign *= sign

        factors: list[list[int]] = [[] for _ in range(p)]
        signs = {vertex: 0 for vertex in range(p)}
        parity = 0
        for label, factor in enumerate(transformed):
            for vertex in factor:
                factors[vertex].append(label)
                for lower_vertex in range(vertex):
                    signs[lower_vertex] += 1
                parity += signs[vertex]
        new_tensor = tuple(tuple(factor) for factor in factors)
        _add_mod(
            result,
            new_tensor,
            coefficient * action_sign * (-1) ** parity,
            p,
        )
    return result


def _alexander_dual(
    tensor_terms: dict[TensorTerm, int], p: int, target_degree: int
) -> dict[TensorTerm, int]:
    result: dict[TensorTerm, int] = {}
    for tensor, coefficient in tensor_terms.items():
        new_tensor: list[Face] = []
        sign = 1
        parity = 0
        for factor in reversed(tensor):
            complement, complement_sign = _simplicial_alexander_dual(
                factor, target_degree
            )
            new_tensor.append(complement)
            sign *= (-1) ** parity * complement_sign
            parity += (target_degree + 1) * len(factor)
        _add_mod(result, tuple(new_tensor), coefficient * sign, p)
    return result
