"""Sparse linear algebra over prime fields, including bit helpers for F_2."""

from __future__ import annotations

Vector = dict[int, int]


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
    scale %= p
    if not scale:
        return answer
    for index, coefficient in right.items():
        value = (answer.get(index, 0) + scale * coefficient) % p
        if value:
            answer[index] = value
        else:
            answer.pop(index, None)
    return answer


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
            vector = vector_add(vector, basis[pivot], p, -coefficient)
    return vector


def rref(rows, p: int) -> dict[int, Vector]:
    """Return a reduced row basis indexed by pivot column."""
    basis: dict[int, Vector] = {}
    for row in rows:
        add_to_basis(basis, row, p)
    return basis


def nullspace(columns: list[Vector], codomain_dimension: int, p: int) -> list[Vector]:
    """Return a basis for the kernel of the map with the given columns."""
    rows = []
    for row_index in range(codomain_dimension):
        row = {
            column_index: column[row_index] % p
            for column_index, column in enumerate(columns)
            if column.get(row_index, 0) % p
        }
        if row:
            rows.append(row)

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
    return cycles


class CoordinateBasis:
    """A reduced basis whose rows carry quotient coordinates."""

    def __init__(self, p: int) -> None:
        self.p = p
        self._rows: dict[int, tuple[Vector, Vector]] = {}

    def add(self, vector: Vector, coordinate: Vector | None = None) -> bool:
        """Add ``vector`` with its quotient ``coordinate`` if independent."""
        if coordinate is None:
            coordinate = {}
        p = self.p
        vector = clean_vector(vector, p)
        coordinate = clean_vector(coordinate, p)
        for pivot in sorted(self._rows, reverse=True):
            coefficient = vector.get(pivot, 0)
            if coefficient:
                row, row_coordinate = self._rows[pivot]
                vector = vector_add(vector, row, p, -coefficient)
                coordinate = vector_add(coordinate, row_coordinate, p, -coefficient)
        if not vector:
            return False

        pivot = leading_index(vector)
        inverse = pow(vector[pivot], -1, p)
        vector = vector_scale(vector, inverse, p)
        coordinate = vector_scale(coordinate, inverse, p)
        for other_pivot, (other_row, other_coordinate) in list(self._rows.items()):
            coefficient = other_row.get(pivot, 0)
            if coefficient:
                self._rows[other_pivot] = (
                    vector_add(other_row, vector, p, -coefficient),
                    vector_add(other_coordinate, coordinate, p, -coefficient),
                )
        self._rows[pivot] = (vector, coordinate)
        return True

    def coordinates(self, vector: Vector) -> Vector:
        """Return quotient coordinates of ``vector`` in this span."""
        p = self.p
        vector = clean_vector(vector, p)
        coordinate: Vector = {}
        for pivot in sorted(self._rows, reverse=True):
            coefficient = vector.get(pivot, 0)
            if coefficient:
                row, row_coordinate = self._rows[pivot]
                vector = vector_add(vector, row, p, -coefficient)
                coordinate = vector_add(coordinate, row_coordinate, p, coefficient)
        if vector:
            raise ValueError("vector is not in the span")
        return coordinate


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

