"""Universal tensor data for odd-primary Steenrod operations."""

from __future__ import annotations

from dataclasses import dataclass

from fastop._odd_primary.indices import OperationIndex

TensorTerm = tuple[tuple[int, ...], ...]


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
