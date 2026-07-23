# Archived: a prime-seven computation

The first reduced power at the prime seven raises the degree-two generator
of complex projective space to its seventh power:

`P¹: H²(CP⁷; F₇) → H¹⁴(CP⁷; F₇)`.

Fastop constructs this space as `Sym⁷(S²)` and computes a rank-one operation.
This is a genuine cochain-level calculation on the simplicial model, not a
result inserted from the known cohomology ring.

## Why degree-lazy construction matters

The model has 13,478,264 nondegenerate cells in total, but the computation
does not need all degrees simultaneously. Symmetric powers now construct
cell labels by degree on demand, and cohomology likewise obtains cell sets
only when a requested degree requires them. The largest materialized layer
for this calculation has 810,810 cells.

On the development arm64 macOS host with Python 3.14.3, the top cohomology
calculation took 72.4 seconds, evaluation of `P¹` took 24.9 seconds, and peak
memory was 1.76 GiB. This is practical as an opt-in benchmark, but deliberately
too large for routine CI.

Reproduce it with:

```bash
python benchmarks/prime_seven.py
```

The established alternatives grow much faster. `Sym⁷(T²)` has 4,464,114,342
cells, while the join-quotient model of `L¹⁵(7)` has 366,127,232 cells. They
should be screened by cell counts rather than constructed with the current
representations.
