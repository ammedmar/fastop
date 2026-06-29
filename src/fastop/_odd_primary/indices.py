"""Index conventions for odd-primary Steenrod operations."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class OperationIndex:
    """Translate ``fastop`` operation indices to the current oddp convention."""

    p: int
    r: int
    source_degree: int
    bockstein: bool = False

    @property
    def oddp_s(self) -> int:
        """Return the operation index used by oddp."""
        return -self.r

    @property
    def oddp_q(self) -> int:
        """Return the cochain degree used by oddp."""
        return -self.source_degree

    @property
    def target_degree(self) -> int:
        """Return the target cohomological degree."""
        return self.source_degree + self.missing_vertices_per_factor

    @property
    def missing_vertices_per_factor(self) -> int:
        """Return how many target vertices each source simplex omits."""
        return 2 * self.r * (self.p - 1) + int(self.bockstein)

