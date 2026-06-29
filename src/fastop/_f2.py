"""Small linear algebra helpers over the field with two elements."""

from __future__ import annotations


def iter_bits(vector: int):
    """Yield the positions of the nonzero entries of ``vector``."""
    while vector:
        bit = vector & -vector
        yield bit.bit_length() - 1
        vector ^= bit


def rank(vectors) -> int:
    """Return the dimension of the span of the given bit vectors."""
    basis: dict[int, int] = {}
    for vector in vectors:
        add_to_basis(basis, vector)
    return len(basis)


def add_to_basis(basis: dict[int, int], vector: int) -> bool:
    """Add ``vector`` to a reduced row basis if it is independent."""
    vector = reduce_vector(basis, vector)
    if not vector:
        return False
    pivot = vector.bit_length() - 1
    for other_pivot, other_vector in list(basis.items()):
        if (other_vector >> pivot) & 1:
            basis[other_pivot] = other_vector ^ vector
    basis[pivot] = vector
    return True


def reduce_vector(basis: dict[int, int], vector: int) -> int:
    """Reduce ``vector`` by a reduced row basis."""
    for pivot in sorted(basis, reverse=True):
        if (vector >> pivot) & 1:
            vector ^= basis[pivot]
    return vector


def rref(rows) -> dict[int, int]:
    """Return a reduced row basis indexed by pivot column."""
    basis: dict[int, int] = {}
    for row in rows:
        add_to_basis(basis, row)
    return basis


def nullspace(columns: list[int], codomain_dimension: int) -> list[int]:
    """Return a basis for the kernel of the linear map with these columns."""
    rows = []
    for row_index in range(codomain_dimension):
        row = 0
        for column_index, column in enumerate(columns):
            if (column >> row_index) & 1:
                row |= 1 << column_index
        if row:
            rows.append(row)

    row_basis = rref(rows)
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


class CoordinateBasis:
    """A reduced basis whose rows carry quotient coordinates."""

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
