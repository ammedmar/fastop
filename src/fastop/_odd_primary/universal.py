"""Universal tensor data for odd-primary Steenrod operations."""

from __future__ import annotations

from dataclasses import dataclass

from fastop._odd_primary.indices import OperationIndex

TensorTerm = tuple[tuple[int, ...], ...]
OmissionPattern = tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class UniversalOperation:
    """Universal tensor formula for one odd-primary operation."""

    p: int
    r: int
    source_degree: int
    bockstein: bool
    target_degree: int
    missing_vertices_per_factor: int
    terms: dict[TensorTerm, int]

    @classmethod
    def from_terms(
        cls,
        index: OperationIndex,
        terms: dict[TensorTerm, int],
    ) -> "UniversalOperation":
        """Create universal data and reduce coefficients modulo ``p``."""
        reduced_terms = {
            tensor: coefficient % index.p
            for tensor, coefficient in terms.items()
            if coefficient % index.p
        }
        return cls(
            p=index.p,
            r=index.r,
            source_degree=index.source_degree,
            bockstein=index.bockstein,
            target_degree=index.target_degree,
            missing_vertices_per_factor=index.missing_vertices_per_factor,
            terms=reduced_terms,
        )

    def signature_table(self) -> "SignatureTable":
        """Return the equivalent omission-pattern coefficient table."""
        return SignatureTable.from_universal(self)


@dataclass(frozen=True)
class SignatureTable:
    """Universal coefficients indexed by omitted target-simplex positions."""

    p: int
    r: int
    source_degree: int
    bockstein: bool
    target_degree: int
    missing_vertices_per_factor: int
    coefficients: dict[OmissionPattern, int]

    @classmethod
    def from_universal(cls, universal: UniversalOperation) -> "SignatureTable":
        """Convert tensor factors to their omitted-position complements."""
        positions = set(range(universal.target_degree + 1))
        coefficients: dict[OmissionPattern, int] = {}
        for tensor, coefficient in universal.terms.items():
            if len(tensor) != universal.p:
                raise ValueError("universal tensor has the wrong number of factors")
            pattern = tuple(
                _omitted_positions(
                    factor,
                    positions,
                    universal.source_degree,
                    universal.missing_vertices_per_factor,
                )
                for factor in tensor
            )
            coefficients[pattern] = (
                coefficients.get(pattern, 0) + coefficient
            ) % universal.p
            if coefficients[pattern] == 0:
                del coefficients[pattern]
        return cls(
            p=universal.p,
            r=universal.r,
            source_degree=universal.source_degree,
            bockstein=universal.bockstein,
            target_degree=universal.target_degree,
            missing_vertices_per_factor=universal.missing_vertices_per_factor,
            coefficients=coefficients,
        )


def native_universal_operation(index: OperationIndex) -> UniversalOperation | None:
    """Return native universal data when the operation family is implemented."""
    if index.bockstein:
        return None
    if index.r <= 0 or index.source_degree != 2 * index.r:
        return None

    factor_length = index.source_degree + 1
    step = index.source_degree
    tensor = tuple(
        tuple(range(step * factor, step * factor + factor_length))
        for factor in range(index.p)
    )
    coefficient = (-1) ** index.r
    return UniversalOperation.from_terms(index, {tensor: coefficient})


def _omitted_positions(
    factor: tuple[int, ...],
    positions: set[int],
    source_degree: int,
    missing_vertices_per_factor: int,
) -> tuple[int, ...]:
    if len(factor) != source_degree + 1:
        raise ValueError("universal tensor factors must have source degree dimension")
    if tuple(sorted(factor)) != factor or len(set(factor)) != len(factor):
        raise ValueError("universal tensor factor positions must be strictly increasing")
    if not set(factor).issubset(positions):
        raise ValueError("universal tensor factor position is outside the target simplex")

    omitted = tuple(sorted(positions.difference(factor)))
    if len(omitted) != missing_vertices_per_factor:
        raise ValueError("universal tensor factor has the wrong codimension")
    return omitted
