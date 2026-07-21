"""Finite simplicial sets with normalized nondegenerate cells."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product
from math import comb
from typing import Iterable, TYPE_CHECKING

from fastop._group_action import (
    CellAction,
    FiniteGroupAction,
    cell_orbits,
    normalize_actions,
    validate_free_action,
    validate_permutations,
)
from fastop.cohomology import PrimeFieldCohomology

if TYPE_CHECKING:
    from fastop.delta_complex import DeltaComplex


@dataclass(frozen=True)
class SimplexReference:
    """A simplex as a nondegenerate cell composed with a surjection.

    ``operator`` is an order-preserving surjection from the vertices of this
    simplex to the local vertices of its underlying nondegenerate cell.
    Repeated adjacent values encode degeneracies.
    """

    degree: int
    cell: int
    operator: tuple[int, ...]

    @property
    def dimension(self) -> int:
        return len(self.operator) - 1

    @property
    def is_nondegenerate(self) -> bool:
        return self.operator == tuple(range(self.degree + 1))


SimplicialFace = int | SimplexReference
SimplicialFaceTable = tuple[tuple[tuple[SimplexReference, ...], ...], ...]


@dataclass(frozen=True)
class SimplicialSet:
    """A finite simplicial set stored through its nondegenerate simplices."""

    face_maps: SimplicialFaceTable
    _product_labels: tuple[tuple[tuple[SimplexReference, ...], ...], ...] | None = field(
        default=None,
        compare=False,
        repr=False,
    )
    _product_factors: tuple["SimplicialSet", ...] | None = field(
        default=None,
        compare=False,
        repr=False,
    )

    def __init__(
        self,
        face_maps: Iterable[Iterable[Iterable[SimplicialFace]]],
        *,
        check: bool = True,
    ):
        raw = tuple(
            tuple(tuple(faces) for faces in degree_faces)
            for degree_faces in face_maps
        )
        while raw and not raw[-1]:
            raw = raw[:-1]
        if not raw or not raw[0]:
            raise ValueError("a simplicial set needs at least one vertex")

        normalized = []
        for degree, degree_faces in enumerate(raw):
            expected_faces = 0 if degree == 0 else degree + 1
            normalized_degree = []
            for faces in degree_faces:
                if len(faces) != expected_faces:
                    raise ValueError(
                        f"each nondegenerate {degree}-simplex needs exactly "
                        f"{expected_faces} faces"
                    )
                normalized_degree.append(tuple(
                    self.nondegenerate(degree - 1, face)
                    if isinstance(face, int) and not isinstance(face, bool)
                    else face
                    for face in faces
                ))
            normalized.append(tuple(normalized_degree))

        table = tuple(normalized)
        self._validate_references(table)
        object.__setattr__(self, "face_maps", table)
        object.__setattr__(self, "_product_labels", None)
        object.__setattr__(self, "_product_factors", None)
        if check:
            self._validate_identities()

    @classmethod
    def from_face_maps(
        cls,
        face_maps: Iterable[Iterable[Iterable[SimplicialFace]]],
        *,
        check: bool = True,
    ) -> "SimplicialSet":
        """Create a finite simplicial set from faces of nondegenerate cells."""
        return cls(face_maps, check=check)

    @classmethod
    def from_delta_complex(cls, complex_: "DeltaComplex") -> "SimplicialSet":
        """Freely add degeneracies to a finite Delta-complex."""
        return cls(complex_.face_maps, check=False)

    @classmethod
    def from_sage(cls, complex_) -> "SimplicialSet":
        """Copy a finite Sage simplicial set without retaining Sage objects."""
        cells = {
            int(degree): tuple(values)
            for degree, values in complex_.cells().items()
            if int(degree) >= 0
        }
        if not cells or not cells.get(0):
            raise ValueError("the Sage simplicial set has no vertices")
        dimension = max(cells)
        indices = {
            degree: {simplex: index for index, simplex in enumerate(cells.get(degree, ()))}
            for degree in range(dimension + 1)
        }
        face_maps = []
        for degree in range(dimension + 1):
            degree_faces = []
            for simplex in cells.get(degree, ()):
                if degree == 0:
                    degree_faces.append(())
                    continue
                faces = []
                for face in complex_.faces(simplex):
                    underlying = face.nondegenerate()
                    underlying_degree = int(underlying.dimension())
                    degeneracies = {int(value) for value in face.degeneracies()}
                    faces.append(SimplexReference(
                        underlying_degree,
                        indices[underlying_degree][underlying],
                        _operator_from_degeneracies(
                            int(face.dimension()),
                            degeneracies,
                        ),
                    ))
                degree_faces.append(tuple(faces))
            face_maps.append(tuple(degree_faces))
        return cls(face_maps)

    @classmethod
    def minimal_sphere(cls, dimension: int) -> "SimplicialSet":
        """Return the one-vertex minimal simplicial model of ``S^dimension``."""
        if not isinstance(dimension, int) or isinstance(dimension, bool):
            raise TypeError("dimension must be an integer")
        if dimension < 1:
            raise ValueError("dimension must be positive")
        face_maps = [[()]] + [[] for _ in range(dimension)]
        if dimension == 1:
            boundary = (0, 0)
        else:
            degenerate_vertex = SimplexReference(
                0,
                0,
                tuple(0 for _ in range(dimension)),
            )
            boundary = tuple(degenerate_vertex for _ in range(dimension + 1))
        face_maps[dimension] = [boundary]
        return cls(face_maps)

    @staticmethod
    def nondegenerate(degree: int, cell: int) -> SimplexReference:
        """Return the identity reference to one nondegenerate cell."""
        return SimplexReference(degree, cell, tuple(range(degree + 1)))

    @property
    def dimension(self) -> int:
        """Return the largest nondegenerate simplex dimension."""
        return len(self.face_maps) - 1

    @property
    def supports_vertex_algorithms(self) -> bool:
        """Return whether cells carry globally comparable vertex sets."""
        return False

    def cells(self, dimension: int | None = None):
        """Return nondegenerate cell indices in one or every dimension."""
        if dimension is None:
            return {
                degree: range(len(degree_faces))
                for degree, degree_faces in enumerate(self.face_maps)
            }
        if dimension < 0 or dimension > self.dimension:
            return ()
        return range(len(self.face_maps[dimension]))

    def f_vector(self) -> tuple[int, ...]:
        """Return the number of nondegenerate cells in every dimension."""
        return tuple(len(degree_faces) for degree_faces in self.face_maps)

    def face(self, degree: int, cell: int, index: int) -> int | None:
        """Return a normalized boundary face, or ``None`` if it is degenerate."""
        reference = self.face_maps[degree][cell][index]
        return reference.cell if reference.is_nondegenerate else None

    def face_reference(
        self,
        reference: SimplexReference,
        index: int,
    ) -> SimplexReference:
        """Apply one face map without prematurely discarding degeneracies."""
        if index < 0 or index > reference.dimension:
            raise IndexError("face index is outside the simplex")
        if reference.is_nondegenerate:
            return self.face_maps[reference.degree][reference.cell][index]

        composite = reference.operator[:index] + reference.operator[index + 1 :]
        image = tuple(dict.fromkeys(composite))
        image_to_position = {value: position for position, value in enumerate(image)}
        quotient_operator = tuple(image_to_position[value] for value in composite)

        if image == tuple(range(reference.degree + 1)):
            restricted = self.nondegenerate(reference.degree, reference.cell)
        else:
            restricted = self._restrict_reference(
                reference.degree,
                reference.cell,
                image,
            )
        return SimplexReference(
            restricted.degree,
            restricted.cell,
            tuple(restricted.operator[value] for value in quotient_operator),
        )

    def restrict(
        self,
        degree: int,
        cell: int,
        positions: tuple[int, ...],
    ) -> int | None:
        """Restrict a cell and return its normalized cell, if nondegenerate."""
        if not positions:
            raise ValueError("a restriction needs at least one retained position")
        if tuple(sorted(set(positions))) != positions:
            raise ValueError("retained positions must be strictly increasing")
        if positions[-1] > degree:
            raise ValueError("retained position is outside the cell")
        reference = self._restrict_reference(degree, cell, positions)
        return reference.cell if reference.is_nondegenerate else None

    def cohomology(
        self,
        p: int = 2,
        *,
        reduced: bool = False,
        convention: int = 1,
    ) -> PrimeFieldCohomology:
        """Return normalized mod-``p`` cohomology with a chosen basis."""
        return PrimeFieldCohomology(self, p=p, reduced=reduced, convention=convention)

    def cartesian_product(self, *others: "SimplicialSet") -> "SimplicialSet":
        """Return the finite simplicial-set cartesian product."""
        if not others:
            return self
        factors = (self,) + tuple(others)
        max_dimension = sum(factor.dimension for factor in factors)
        base_options = [
            tuple(
                (degree, cell)
                for degree, degree_faces in enumerate(factor.face_maps)
                for cell in range(len(degree_faces))
            )
            for factor in factors
        ]
        labels_by_degree = []
        for total_degree in range(max_dimension + 1):
            labels = []
            for bases in product(*base_options):
                degrees = tuple(degree for degree, _ in bases)
                if max(degrees) > total_degree or sum(degrees) < total_degree:
                    continue
                operator_options = tuple(
                    _surjections(total_degree, degree)
                    for degree in degrees
                )
                for operators in product(*operator_options):
                    if _common_degeneracies(operators):
                        continue
                    labels.append(tuple(
                        SimplexReference(degree, cell, operator)
                        for (degree, cell), operator in zip(bases, operators)
                    ))
            labels_by_degree.append(tuple(labels))

        label_indices = tuple(
            {label: index for index, label in enumerate(labels)}
            for labels in labels_by_degree
        )
        product_faces = [tuple(() for _ in labels_by_degree[0])]
        for degree in range(1, max_dimension + 1):
            product_faces.append(tuple(
                tuple(
                    _normalize_product_reference(
                        tuple(
                            factor.face_reference(component, index)
                            for factor, component in zip(factors, label)
                        ),
                        label_indices,
                    )
                    for index in range(degree + 1)
                )
                for label in labels_by_degree[degree]
            ))

        answer = SimplicialSet(product_faces, check=False)
        object.__setattr__(answer, "_product_labels", tuple(labels_by_degree))
        object.__setattr__(answer, "_product_factors", factors)
        return answer

    def factor_permutation_action(self, permutation: Iterable[int]) -> CellAction:
        """Return the action permuting factors of a stored cartesian product."""
        if self._product_labels is None or self._product_factors is None:
            raise ValueError("factor permutations require a cartesian product")
        permutation = tuple(permutation)
        factor_count = len(self._product_factors)
        if set(permutation) != set(range(factor_count)):
            raise ValueError("factor action must be a permutation")
        for target, source in enumerate(permutation):
            if self._product_factors[target] != self._product_factors[source]:
                raise ValueError("only equal factors may be permuted")

        label_indices = tuple(
            {label: index for index, label in enumerate(labels)}
            for labels in self._product_labels
        )
        return tuple(
            tuple(
                label_indices[degree][tuple(label[index] for index in permutation)]
                for label in labels
            )
            for degree, labels in enumerate(self._product_labels)
        )

    def symmetric_power(
        self,
        power: int,
    ) -> "SimplicialSet | SymmetricPowerSimplicialSet":
        """Return the quotient of a cartesian power by factor permutations.

        Orbit representatives are built directly as unordered tuples.  This
        avoids materializing the much larger ordered cartesian power first.
        """
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError("power must be an integer")
        if power < 1:
            raise ValueError("power must be positive")
        if power == 1:
            return self
        return SymmetricPowerSimplicialSet(self, power)

    def symmetric_power_f_vector(self, power: int) -> tuple[int, ...]:
        """Count a symmetric power without constructing its face model."""
        _validate_symmetric_power(power)
        return _symmetric_power_f_vector(self.f_vector(), power)

    def quotient(
        self,
        generators: FiniteGroupAction | Iterable[Iterable[Iterable[int]]],
        *,
        require_free: bool = False,
    ) -> "SimplicialSet":
        """Return the quotient by strict finite simplicial automorphisms."""
        actions = normalize_actions(generators)
        if not actions:
            return self
        counts = self.f_vector()
        validate_permutations(actions, counts)
        for action in actions:
            self._validate_action(action)
        if require_free:
            validate_free_action(actions, counts)

        orbit_indices, orbit_representatives = cell_orbits(actions, counts)
        quotient_faces = [tuple(() for _ in orbit_representatives[0])]
        for degree in range(1, self.dimension + 1):
            quotient_faces.append(tuple(
                tuple(
                    SimplexReference(
                        reference.degree,
                        orbit_indices[reference.degree][reference.cell],
                        reference.operator,
                    )
                    for reference in self.face_maps[degree][representative]
                )
                for representative in orbit_representatives[degree]
            ))
        return SimplicialSet(quotient_faces)

    def _restrict_reference(
        self,
        degree: int,
        cell: int,
        positions: tuple[int, ...],
    ) -> SimplexReference:
        reference = self.nondegenerate(degree, cell)
        omitted = set(range(degree + 1)).difference(positions)
        for index in sorted(omitted, reverse=True):
            reference = self.face_reference(reference, index)
        return reference

    def _validate_references(self, table: SimplicialFaceTable) -> None:
        for degree, degree_faces in enumerate(table):
            for faces in degree_faces:
                for reference in faces:
                    if not isinstance(reference, SimplexReference):
                        raise TypeError("simplicial faces must be cell references")
                    if reference.dimension != degree - 1:
                        raise ValueError("simplicial face has the wrong dimension")
                    if reference.degree < 0 or reference.degree >= len(table):
                        raise ValueError("underlying simplex degree is outside the model")
                    if reference.cell < 0 or reference.cell >= len(table[reference.degree]):
                        raise ValueError("underlying simplex index is outside the model")
                    _validate_operator(reference)

    def _validate_identities(self) -> None:
        for degree in range(2, self.dimension + 1):
            for cell in self.cells(degree):
                reference = self.nondegenerate(degree, cell)
                for i in range(degree):
                    for j in range(i + 1, degree + 1):
                        left = self.face_reference(
                            self.face_reference(reference, j),
                            i,
                        )
                        right = self.face_reference(
                            self.face_reference(reference, i),
                            j - 1,
                        )
                        if left != right:
                            raise ValueError(
                                "face maps violate the simplicial identities "
                                f"at degree {degree}, cell {cell}, i={i}, j={j}"
                            )

    def _validate_action(self, action: CellAction) -> None:
        for degree in range(1, self.dimension + 1):
            for cell in self.cells(degree):
                for index, reference in enumerate(self.face_maps[degree][cell]):
                    acted_reference = SimplexReference(
                        reference.degree,
                        action[reference.degree][reference.cell],
                        reference.operator,
                    )
                    face_of_acted = self.face_maps[degree][action[degree][cell]][index]
                    if acted_reference != face_of_acted:
                        raise ValueError("cell action does not commute with face maps")


def _validate_operator(reference: SimplexReference) -> None:
    operator = reference.operator
    if not operator or operator[0] != 0 or operator[-1] != reference.degree:
        raise ValueError("simplex operator must be surjective")
    if any(
        right < left or right > left + 1
        for left, right in zip(operator, operator[1:])
    ):
        raise ValueError("simplex operator must be order preserving")


@lru_cache(maxsize=None)
def _surjections(total_degree: int, base_degree: int) -> tuple[tuple[int, ...], ...]:
    if base_degree < 0 or base_degree > total_degree:
        return ()
    answer = []
    for increments in combinations(range(total_degree), base_degree):
        increment_positions = set(increments)
        value = 0
        operator = [value]
        for position in range(total_degree):
            if position in increment_positions:
                value += 1
            operator.append(value)
        answer.append(tuple(operator))
    return tuple(answer)


def _degeneracy_positions(operator: tuple[int, ...]) -> set[int]:
    return {
        position
        for position in range(len(operator) - 1)
        if operator[position] == operator[position + 1]
    }


def _common_degeneracies(operators: tuple[tuple[int, ...], ...]) -> set[int]:
    if not operators:
        return set()
    common = _degeneracy_positions(operators[0])
    for operator in operators[1:]:
        common.intersection_update(_degeneracy_positions(operator))
    return common


def _operator_from_degeneracies(
    total_degree: int,
    degeneracies: set[int],
) -> tuple[int, ...]:
    value = 0
    operator = [value]
    for position in range(total_degree):
        if position not in degeneracies:
            value += 1
        operator.append(value)
    return tuple(operator)


def _normalize_product_reference(
    components: tuple[SimplexReference, ...],
    label_indices: tuple[dict[tuple[SimplexReference, ...], int], ...],
) -> SimplexReference:
    total_degree = components[0].dimension
    common = _common_degeneracies(tuple(component.operator for component in components))
    quotient_operator = _operator_from_degeneracies(total_degree, common)
    base_degree = total_degree - len(common)

    representatives = []
    previous = None
    for position, value in enumerate(quotient_operator):
        if value != previous:
            representatives.append(position)
            previous = value
    base_label = tuple(
        SimplexReference(
            component.degree,
            component.cell,
            tuple(component.operator[position] for position in representatives),
        )
        for component in components
    )
    return SimplexReference(
        base_degree,
        label_indices[base_degree][base_label],
        quotient_operator,
    )


def _symmetric_index_labels(
    references: tuple[SimplexReference, ...],
    power: int,
    *,
    total_degree: int,
    max_base_degree: int,
):
    """Yield unordered tuples with no common degeneracy.

    If ``r`` factors remain, they can remove at most
    ``r * max_base_degree`` common degeneracy positions.  Applying this
    bound during generation avoids visiting the overwhelmingly many tuples
    that can only remain degenerate, especially in fifth symmetric powers.
    """
    masks = tuple(
        sum(
            1 << position
            for position in _degeneracy_positions(reference.operator)
        )
        for reference in references
    )
    prefix: list[int] = []

    def generate(start: int, remaining: int, common_mask: int):
        if remaining == 0:
            if common_mask == 0:
                yield tuple(prefix)
            return
        if common_mask.bit_count() > remaining * max_base_degree:
            return
        if common_mask == 0:
            for suffix in combinations_with_replacement(
                range(start, len(references)),
                remaining,
            ):
                yield tuple(prefix) + suffix
            return

        for index in range(start, len(references)):
            next_mask = common_mask & masks[index]
            if next_mask.bit_count() > (remaining - 1) * max_base_degree:
                continue
            prefix.append(index)
            yield from generate(index, remaining - 1, next_mask)
            prefix.pop()

    yield from generate(0, power, (1 << total_degree) - 1)


class SymmetricPowerSimplicialSet:
    """A finite symmetric power with faces computed from compact cell labels."""

    def __init__(self, base, power: int):
        _validate_symmetric_power(power)
        self.base = base
        self.power = power
        self._dimension = power * base.dimension
        self._references: dict[int, tuple[SimplexReference, ...]] = {}
        self._reference_indices: dict[int, dict[SimplexReference, int]] = {}
        self._labels: dict[int, tuple[tuple[int, ...], ...]] = {}
        self._label_indices: dict[int, dict[tuple[int, ...], int]] = {}

    @staticmethod
    def nondegenerate(degree: int, cell: int) -> SimplexReference:
        """Return the identity reference to one nondegenerate cell."""
        return SimplexReference(degree, cell, tuple(range(degree + 1)))

    @property
    def dimension(self) -> int:
        """Return the dimension of the symmetric power."""
        return self._dimension

    @property
    def supports_vertex_algorithms(self) -> bool:
        """Return whether cells carry globally comparable vertex sets."""
        return False

    def cells(self, dimension: int | None = None):
        """Return nondegenerate cell indices in one or every dimension."""
        if dimension is None:
            return {
                degree: self.cells(degree)
                for degree in range(self.dimension + 1)
            }
        if dimension < 0 or dimension > self.dimension:
            return ()
        self._ensure_degree(dimension)
        return range(len(self._labels[dimension]))

    def f_vector(self) -> tuple[int, ...]:
        """Return the number of nondegenerate cells in every dimension."""
        return _symmetric_power_f_vector(self.base.f_vector(), self.power)

    def face(self, degree: int, cell: int, index: int) -> int | None:
        """Return a normalized boundary face, or ``None`` if degenerate."""
        reference = self.face_reference(self.nondegenerate(degree, cell), index)
        return reference.cell if reference.is_nondegenerate else None

    def face_reference(
        self,
        reference: SimplexReference,
        index: int,
    ) -> SimplexReference:
        """Apply one face map without prematurely discarding degeneracies."""
        if index < 0 or index > reference.dimension:
            raise IndexError("face index is outside the simplex")
        if reference.is_nondegenerate:
            self._ensure_degree(reference.degree)
            components = tuple(
                self.base.face_reference(
                    self._references[reference.degree][component],
                    index,
                )
                for component in self._labels[reference.degree][reference.cell]
            )
            return self._normalize_components(components)

        composite = reference.operator[:index] + reference.operator[index + 1 :]
        image = tuple(dict.fromkeys(composite))
        image_to_position = {value: position for position, value in enumerate(image)}
        quotient_operator = tuple(image_to_position[value] for value in composite)

        if image == tuple(range(reference.degree + 1)):
            restricted = self.nondegenerate(reference.degree, reference.cell)
        else:
            restricted = self._restrict_reference(
                reference.degree,
                reference.cell,
                image,
            )
        return SimplexReference(
            restricted.degree,
            restricted.cell,
            tuple(restricted.operator[value] for value in quotient_operator),
        )

    def restrict(
        self,
        degree: int,
        cell: int,
        positions: tuple[int, ...],
    ) -> int | None:
        """Restrict a cell and return its normalized cell, if nondegenerate."""
        if not positions:
            raise ValueError("a restriction needs at least one retained position")
        if tuple(sorted(set(positions))) != positions:
            raise ValueError("retained positions must be strictly increasing")
        if positions[-1] > degree:
            raise ValueError("retained position is outside the cell")
        reference = self._restrict_reference(degree, cell, positions)
        return reference.cell if reference.is_nondegenerate else None

    def cohomology(
        self,
        p: int = 2,
        *,
        reduced: bool = False,
        convention: int = 1,
    ) -> PrimeFieldCohomology:
        """Return normalized mod-``p`` cohomology with a chosen basis."""
        return PrimeFieldCohomology(self, p=p, reduced=reduced, convention=convention)

    def symmetric_power(self, power: int):
        """Return a further symmetric power of this finite simplicial set."""
        if not isinstance(power, int) or isinstance(power, bool):
            raise TypeError("power must be an integer")
        if power < 1:
            raise ValueError("power must be positive")
        if power == 1:
            return self
        return SymmetricPowerSimplicialSet(self, power)

    def symmetric_power_f_vector(self, power: int) -> tuple[int, ...]:
        """Count a further symmetric power without constructing it."""
        _validate_symmetric_power(power)
        return _symmetric_power_f_vector(self.f_vector(), power)

    def _restrict_reference(
        self,
        degree: int,
        cell: int,
        positions: tuple[int, ...],
    ) -> SimplexReference:
        self._ensure_degree(degree)
        components = tuple(
            _restrict_reference_on_model(
                self.base,
                self._references[degree][component],
                positions,
            )
            for component in self._labels[degree][cell]
        )
        return self._normalize_components(components)

    def _normalize_components(
        self,
        components: tuple[SimplexReference, ...],
    ) -> SimplexReference:
        total_degree = components[0].dimension
        common = _common_degeneracies(tuple(
            component.operator for component in components
        ))
        quotient_operator = _operator_from_degeneracies(total_degree, common)
        base_degree = total_degree - len(common)
        self._ensure_degree(base_degree)

        representatives = []
        previous = None
        for position, value in enumerate(quotient_operator):
            if value != previous:
                representatives.append(position)
                previous = value
        label = tuple(sorted(
            self._reference_indices[base_degree][SimplexReference(
                component.degree,
                component.cell,
                tuple(
                    component.operator[position]
                    for position in representatives
                ),
            )]
            for component in components
        ))
        return SimplexReference(
            base_degree,
            self._label_indices[base_degree][label],
            quotient_operator,
        )

    def _ensure_degree(self, total_degree: int) -> None:
        if total_degree in self._labels:
            return
        references = tuple(sorted(
            (
                SimplexReference(degree, cell, operator)
                for degree in range(min(self.base.dimension, total_degree) + 1)
                for cell in self.base.cells(degree)
                for operator in _surjections(total_degree, degree)
            ),
            key=_reference_key,
        ))
        labels = tuple(_symmetric_index_labels(
            references,
            self.power,
            total_degree=total_degree,
            max_base_degree=self.base.dimension,
        ))
        self._references[total_degree] = references
        self._reference_indices[total_degree] = {
            reference: index for index, reference in enumerate(references)
        }
        self._labels[total_degree] = labels
        self._label_indices[total_degree] = {
            label: index for index, label in enumerate(labels)
        }


def _reference_key(reference: SimplexReference):
    return reference.degree, reference.cell, reference.operator


def _restrict_reference_on_model(model, reference, positions):
    omitted = set(range(reference.dimension + 1)).difference(positions)
    for index in sorted(omitted, reverse=True):
        reference = model.face_reference(reference, index)
    return reference


def _validate_symmetric_power(power: int) -> None:
    if not isinstance(power, int) or isinstance(power, bool):
        raise TypeError("power must be an integer")
    if power < 1:
        raise ValueError("power must be positive")


def _symmetric_power_f_vector(
    base_f_vector: tuple[int, ...],
    power: int,
) -> tuple[int, ...]:
    dimension = len(base_f_vector) - 1
    answer = []
    for total_degree in range(power * dimension + 1):
        count = 0
        for common_count in range(total_degree + 1):
            available = sum(
                cell_count * comb(total_degree - common_count, degree)
                for degree, cell_count in enumerate(base_f_vector)
                if degree <= total_degree - common_count
            )
            count += (
                (-1) ** common_count
                * comb(total_degree, common_count)
                * comb(available + power - 1, power)
            )
        answer.append(count)
    return tuple(answer)
