"""Universal tensor data for Steenrod operations."""

from __future__ import annotations

from dataclasses import dataclass

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
        *,
        p: int,
        r: int,
        source_degree: int,
        bockstein: bool,
        target_degree: int,
        missing_vertices_per_factor: int,
        terms: dict[TensorTerm, int],
    ) -> "UniversalOperation":
        """Create universal data and reduce coefficients modulo ``p``."""
        reduced_terms = {
            tensor: coefficient % p
            for tensor, coefficient in terms.items()
            if coefficient % p
        }
        return cls(
            p=p,
            r=r,
            source_degree=source_degree,
            bockstein=bockstein,
            target_degree=target_degree,
            missing_vertices_per_factor=missing_vertices_per_factor,
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


def native_universal_operation(
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    missing_vertices_per_factor: int,
) -> UniversalOperation | None:
    """Return native universal data when the operation family is implemented."""
    if bockstein:
        return None
    if r <= 0 or source_degree != 2 * r:
        return None

    factor_length = source_degree + 1
    step = source_degree
    tensor = tuple(
        tuple(range(step * factor, step * factor + factor_length))
        for factor in range(p)
    )
    coefficient = (-1) ** r
    return UniversalOperation.from_terms(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
        terms={tensor: coefficient},
    )


def universal_operation(
    *,
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    missing_vertices_per_factor: int,
    oddp_s: int | None = None,
    oddp_q: int | None = None,
) -> UniversalOperation:
    """Build universal tensor data, falling back to oddp when needed."""
    native = native_universal_operation(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
    )
    if native is not None:
        return native

    if oddp_s is None or oddp_q is None:
        raise NotImplementedError("this universal operation is not implemented natively")

    from fastop._oddp_bridge import universal_terms_oddp

    return UniversalOperation.from_terms(
        p=p,
        r=r,
        source_degree=source_degree,
        bockstein=bockstein,
        target_degree=target_degree,
        missing_vertices_per_factor=missing_vertices_per_factor,
        terms=universal_terms_oddp(
            p=p,
            bockstein=bockstein,
            oddp_s=oddp_s,
            oddp_q=oddp_q,
        ),
    )


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
