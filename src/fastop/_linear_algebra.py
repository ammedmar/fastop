"""Sparse linear algebra over prime fields."""

from __future__ import annotations

from bisect import insort

Vector = dict[int, int]

try:
    from fastop._native import (
        column_image_and_kernel_basis as _native_column_image_and_kernel_basis,
    )
    from fastop._native import coordinate_basis_from_vectors as _native_coordinate_basis_from_vectors
except ImportError:  # pragma: no cover - depends on optional extension build
    _native_column_image_and_kernel_basis = None
    _native_coordinate_basis_from_vectors = None


def is_prime(value: int) -> bool:
    """Return whether ``value`` is prime."""
    if not isinstance(value, int) or isinstance(value, bool):
        return False
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, int(value**0.5) + 1))


def clean_vector(vector: Vector, p: int) -> Vector:
    """Return ``vector`` with coefficients reduced modulo ``p``."""
    return {
        index: coefficient % p
        for index, coefficient in vector.items()
        if coefficient % p
    }


def vector_add(left: Vector, right: Vector, p: int, scale: int = 1) -> Vector:
    """Return ``left + scale * right`` over ``F_p``."""
    answer = dict(left)
    return _vector_add_inplace(answer, right, p, scale)


def _vector_add_inplace(left: Vector, right: Vector, p: int, scale: int = 1) -> Vector:
    """Add ``scale * right`` into ``left`` over ``F_p``."""
    scale %= p
    if not scale:
        return left
    for index, coefficient in right.items():
        value = (left.get(index, 0) + scale * coefficient) % p
        if value:
            left[index] = value
        else:
            left.pop(index, None)
    return left


def vector_scale(vector: Vector, scalar: int, p: int) -> Vector:
    """Return ``scalar * vector`` over ``F_p``."""
    scalar %= p
    if not scalar:
        return {}
    return {index: (scalar * coefficient) % p for index, coefficient in vector.items()}


def leading_index(vector: Vector) -> int:
    """Return the largest index with nonzero coefficient."""
    return max(vector)


def rank(vectors, p: int) -> int:
    """Return the dimension of the span of ``vectors`` over ``F_p``."""
    basis: dict[int, Vector] = {}
    for vector in vectors:
        add_to_basis(basis, vector, p)
    return len(basis)


def add_to_basis(basis: dict[int, Vector], vector: Vector, p: int) -> bool:
    """Add ``vector`` to a reduced row basis if independent."""
    vector = reduce_vector(basis, vector, p)
    if not vector:
        return False
    pivot = leading_index(vector)
    inverse = pow(vector[pivot], -1, p)
    vector = vector_scale(vector, inverse, p)
    for other_pivot, other_vector in list(basis.items()):
        coefficient = other_vector.get(pivot, 0)
        if coefficient:
            basis[other_pivot] = vector_add(other_vector, vector, p, -coefficient)
    basis[pivot] = vector
    return True


def reduce_vector(basis: dict[int, Vector], vector: Vector, p: int) -> Vector:
    """Reduce ``vector`` by a reduced row basis."""
    vector = clean_vector(vector, p)
    for pivot in sorted(basis, reverse=True):
        coefficient = vector.get(pivot, 0)
        if coefficient:
            _vector_add_inplace(vector, basis[pivot], p, -coefficient)
    return vector


def rref(rows, p: int) -> dict[int, Vector]:
    """Return a reduced row basis indexed by pivot column."""
    basis: dict[int, Vector] = {}
    for row in rows:
        add_to_basis(basis, row, p)
    return basis


def column_kernel_basis(columns: list[Vector], p: int) -> list[Vector]:
    """Return kernel vectors by TDA-style column reduction."""
    return column_image_and_kernel_basis(columns, p)[1]


def column_image_and_kernel_basis(
    columns: list[Vector],
    p: int,
) -> tuple[list[Vector], list[Vector]]:
    """Return image and kernel bases by TDA-style column reduction."""
    if _native_column_image_and_kernel_basis is not None:
        return _native_column_image_and_kernel_basis(columns, p)
    return _column_image_and_kernel_basis_python(columns, p)


def _column_image_and_kernel_basis_python(
    columns: list[Vector],
    p: int,
) -> tuple[list[Vector], list[Vector]]:
    """Python fallback for TDA-style column reduction."""
    reduced_columns: list[Vector] = []
    transforms: list[Vector] = []
    pivot_to_reduced_column: dict[int, int] = {}
    cycles = []

    for column_index, column in enumerate(columns):
        reduced = clean_vector(column, p)
        transform = {column_index: 1}
        while reduced:
            pivot = leading_index(reduced)
            pivot_column_index = pivot_to_reduced_column.get(pivot)
            if pivot_column_index is None:
                break
            pivot_column = reduced_columns[pivot_column_index]
            coefficient = reduced[pivot] * pow(pivot_column[pivot], -1, p)
            _vector_add_inplace(reduced, pivot_column, p, -coefficient)
            _vector_add_inplace(transform, transforms[pivot_column_index], p, -coefficient)

        if reduced:
            pivot = leading_index(reduced)
            inverse = pow(reduced[pivot], -1, p)
            reduced = vector_scale(reduced, inverse, p)
            transform = vector_scale(transform, inverse, p)
            pivot_to_reduced_column[pivot] = len(reduced_columns)
            reduced_columns.append(reduced)
            transforms.append(transform)
        else:
            cycles.append(clean_vector(transform, p))

    return reduced_columns, cycles


def nullspace(columns: list[Vector], codomain_dimension: int, p: int) -> list[Vector]:
    """Return a basis for the kernel of the map with the given columns."""
    return _row_basis_and_nullspace(columns, codomain_dimension, p)[1]


def _row_basis_and_nullspace(
    columns: list[Vector],
    codomain_dimension: int,
    p: int,
) -> tuple[dict[int, Vector], list[Vector]]:
    """Return a row basis and kernel basis for a column-defined map."""
    rows_by_index: list[Vector] = [{} for _ in range(codomain_dimension)]
    for column_index, column in enumerate(columns):
        for row_index, coefficient in column.items():
            value = coefficient % p
            if value:
                rows_by_index[row_index][column_index] = value
    rows = [row for row in rows_by_index if row]

    row_basis = rref(rows, p)
    pivots = set(row_basis)
    cycles = []
    for free_column in range(len(columns)):
        if free_column in pivots:
            continue
        vector = {free_column: 1}
        for pivot, row in row_basis.items():
            coefficient = row.get(free_column, 0)
            if coefficient:
                vector[pivot] = (-coefficient) % p
        cycles.append(clean_vector(vector, p))
    return row_basis, cycles


class CoordinateBasis:
    """A pivot basis whose rows carry quotient coordinates."""

    def __init__(self, p: int) -> None:
        self.p = p
        self._rows: dict[int, tuple[Vector, Vector]] = {}
        self._pivots: list[int] = []
        self._pivot_set: set[int] = set()

    @classmethod
    def from_rows(
        cls,
        p: int,
        rows: dict[int, tuple[Vector, Vector]],
    ) -> "CoordinateBasis":
        """Create a coordinate basis from precomputed pivot rows."""
        basis = cls(p)
        basis._rows = rows
        basis._pivots = sorted(rows)
        basis._pivot_set = set(rows)
        return basis

    def add(self, vector: Vector, coordinate: Vector | None = None) -> bool:
        """Add ``vector`` with its quotient ``coordinate`` if independent."""
        if coordinate is None:
            coordinate = {}
        p = self.p
        vector = clean_vector(vector, p)
        coordinate = clean_vector(coordinate, p)
        self._reduce_with_coordinates(vector, coordinate, -1)
        if not vector:
            return False

        pivot = leading_index(vector)
        inverse = pow(vector[pivot], -1, p)
        vector = vector_scale(vector, inverse, p)
        coordinate = vector_scale(coordinate, inverse, p)
        self._rows[pivot] = (vector, coordinate)
        insort(self._pivots, pivot)
        self._pivot_set.add(pivot)
        return True

    def add_vector(self, vector: Vector) -> bool:
        """Add ``vector`` with zero quotient coordinate if independent."""
        p = self.p
        vector = clean_vector(vector, p)
        self._reduce_vector(vector)
        if not vector:
            return False

        pivot = leading_index(vector)
        inverse = pow(vector[pivot], -1, p)
        vector = vector_scale(vector, inverse, p)
        self._rows[pivot] = (vector, {})
        insort(self._pivots, pivot)
        self._pivot_set.add(pivot)
        return True

    def add_reduced_rows(self, rows: dict[int, Vector]) -> None:
        """Add already reduced rows with zero quotient coordinates."""
        for pivot, row in rows.items():
            self._rows[pivot] = (row, {})
        self._pivots = sorted(self._rows)
        self._pivot_set = set(self._rows)

    def coordinates(self, vector: Vector) -> Vector:
        """Return quotient coordinates of ``vector`` in this span."""
        p = self.p
        vector = clean_vector(vector, p)
        coordinate: Vector = {}
        self._reduce_with_coordinates(vector, coordinate, 1)
        if vector:
            raise ValueError("vector is not in the span")
        return coordinate

    def _reduce_with_coordinates(
        self,
        vector: Vector,
        coordinate: Vector,
        coordinate_sign: int,
    ) -> None:
        """Reduce ``vector`` and update quotient ``coordinate`` in place."""
        p = self.p
        candidate_pivots = set(vector).intersection(self._pivot_set)
        while candidate_pivots:
            pivot = max(candidate_pivots)
            candidate_pivots.remove(pivot)
            coefficient = vector.get(pivot, 0)
            if not coefficient:
                continue
            row, row_coordinate = self._rows[pivot]
            _vector_add_inplace(vector, row, p, -coefficient)
            _vector_add_inplace(
                coordinate,
                row_coordinate,
                p,
                coordinate_sign * coefficient,
            )
            candidate_pivots.update(set(row).intersection(vector, self._pivot_set))

    def _reduce_vector(self, vector: Vector) -> None:
        """Reduce ``vector`` in place without tracking quotient coordinates."""
        p = self.p
        candidate_pivots = set(vector).intersection(self._pivot_set)
        while candidate_pivots:
            pivot = max(candidate_pivots)
            candidate_pivots.remove(pivot)
            coefficient = vector.get(pivot, 0)
            if not coefficient:
                continue
            row, _ = self._rows[pivot]
            _vector_add_inplace(vector, row, p, -coefficient)
            candidate_pivots.update(set(row).intersection(vector, self._pivot_set))


def coordinate_basis_from_vectors(
    p: int,
    boundary_vectors: list[Vector],
    cycles: list[Vector],
) -> tuple[list[Vector], CoordinateBasis]:
    """Build cohomology representatives and their quotient projector."""
    if _native_coordinate_basis_from_vectors is not None:
        cocycle_basis, rows = _native_coordinate_basis_from_vectors(
            boundary_vectors,
            cycles,
            p,
        )
        return cocycle_basis, CoordinateBasis.from_rows(p, rows)
    return _coordinate_basis_from_vectors_python(p, boundary_vectors, cycles)


def _coordinate_basis_from_vectors_python(
    p: int,
    boundary_vectors: list[Vector],
    cycles: list[Vector],
) -> tuple[list[Vector], CoordinateBasis]:
    projector = CoordinateBasis(p)
    for boundary in boundary_vectors:
        projector.add_vector(boundary)

    cocycle_basis = []
    for cycle in cycles:
        basis_coordinate = {len(cocycle_basis): 1}
        if projector.add(cycle, basis_coordinate):
            cocycle_basis.append(cycle)
    return cocycle_basis, projector
