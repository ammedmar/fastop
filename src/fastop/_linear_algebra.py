"""Sparse linear algebra over prime fields, including bit helpers for F_2."""

from __future__ import annotations

from bisect import insort

Vector = dict[int, int]

try:
    from fastop._native import (
        column_image_and_kernel_basis as _native_column_image_and_kernel_basis,
    )
except ImportError:  # pragma: no cover - depends on optional extension build
    _native_column_image_and_kernel_basis = None


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


def iter_bits(vector: int):
    """Yield the positions of the nonzero entries of a bit vector."""
    while vector:
        bit = vector & -vector
        yield bit.bit_length() - 1
        vector ^= bit


def bit_rank(vectors) -> int:
    """Return the dimension of a span of bit vectors over F_2."""
    basis: dict[int, int] = {}
    for vector in vectors:
        bit_add_to_basis(basis, vector)
    return len(basis)


def bit_add_to_basis(basis: dict[int, int], vector: int) -> bool:
    """Add a bit vector to a reduced row basis if independent."""
    vector = bit_reduce_vector(basis, vector)
    if not vector:
        return False
    pivot = vector.bit_length() - 1
    for other_pivot, other_vector in list(basis.items()):
        if (other_vector >> pivot) & 1:
            basis[other_pivot] = other_vector ^ vector
    basis[pivot] = vector
    return True


def bit_reduce_vector(basis: dict[int, int], vector: int) -> int:
    """Reduce a bit vector by a reduced row basis."""
    for pivot in sorted(basis, reverse=True):
        if (vector >> pivot) & 1:
            vector ^= basis[pivot]
    return vector


def bit_rref(rows) -> dict[int, int]:
    """Return a reduced row basis of bit vectors indexed by pivot column."""
    basis: dict[int, int] = {}
    for row in rows:
        bit_add_to_basis(basis, row)
    return basis


def bit_nullspace(columns: list[int], codomain_dimension: int) -> list[int]:
    """Return a basis for the kernel of a bit-vector linear map."""
    rows = []
    for row_index in range(codomain_dimension):
        row = 0
        for column_index, column in enumerate(columns):
            if (column >> row_index) & 1:
                row |= 1 << column_index
        if row:
            rows.append(row)

    row_basis = bit_rref(rows)
    pivots = set(row_basis)
    domain_dimension = len(columns)
    cycles = []
    for free_column in range(domain_dimension):
        if free_column in pivots:
            continue
        vector = 1 << free_column
        for pivot, row in row_basis.items():
            if (row >> free_column) & 1:
                vector |= 1 << pivot
        cycles.append(vector)
    return cycles


class BitCoordinateBasis:
    """A reduced F_2 bit basis whose rows carry quotient coordinates."""

    def __init__(self) -> None:
        self._rows: dict[int, tuple[int, int]] = {}

    def add(self, vector: int, coordinate: int = 0) -> bool:
        """Add ``vector`` with its quotient ``coordinate`` if independent."""
        for pivot in sorted(self._rows, reverse=True):
            if (vector >> pivot) & 1:
                row, row_coordinate = self._rows[pivot]
                vector ^= row
                coordinate ^= row_coordinate
        if not vector:
            return False

        pivot = vector.bit_length() - 1
        for other_pivot, (other_row, other_coordinate) in list(self._rows.items()):
            if (other_row >> pivot) & 1:
                self._rows[other_pivot] = (
                    other_row ^ vector,
                    other_coordinate ^ coordinate,
                )
        self._rows[pivot] = (vector, coordinate)
        return True

    def coordinates(self, vector: int) -> int:
        """Return quotient coordinates of ``vector`` in this span."""
        coordinate = 0
        for pivot in sorted(self._rows, reverse=True):
            if (vector >> pivot) & 1:
                row, row_coordinate = self._rows[pivot]
                vector ^= row
                coordinate ^= row_coordinate
        if vector:
            raise ValueError("vector is not in the span")
        return coordinate
