# Archived: symmetric products of surfaces

For a closed orientable surface Σ_g, the package constructs a compact
one-vertex simplicial-set model and then forms its symmetric powers directly.
The cubic family

`Sym³(Σ_g) = (Σ_g)³ / S₃`

consists of smooth six-manifolds. It is a useful controlled experiment: the
dimension and target operation stay fixed while the topology and cell count
grow with genus.

## Ground truth

Macdonald's presentation gives the Betti numbers of Sym³(Σ_g) as

```text
(1,
 2g,
 C(2g, 2) + 1,
 C(2g, 3) + 2g,
 C(2g, 2) + 1,
 2g,
 1).
```

A modern complete treatment of the integral presentation is available in
[Gugnin, *On Integral Cohomology Ring of Symmetric Products*](https://arxiv.org/abs/1502.01862).
At the prime 3, the degree-two divisor class supplies the expected nonzero
P¹: H² → H⁶. The computation does not insert that class by hand: it constructs
cohomology from the face maps and evaluates the cochain-level operation on a
basis.

## Computed family

| Genus | Nondegenerate cells | Betti numbers | P¹ rank | Status |
| ---: | ---: | --- | ---: | --- |
| 0 | 84 | `(1, 0, 1, 0, 1, 0, 1)` | 1 | CP³ ground truth |
| 1 | 1,638 | `(1, 2, 2, 2, 2, 2, 1)` | 1 | full regression |
| 2 | 41,478 | `(1, 4, 7, 8, 7, 4, 1)` | 1 | full regression |
| 3 | 189,626 | `(1, 6, 16, 26, 16, 6, 1)` expected | 1 | extended operation run |

On the development arm64 macOS host with Python 3.14.3, genus two built in
3.33 seconds, computed full cohomology in 4.48 seconds, and computed the
operation from a fresh cohomology object in 0.41 seconds. Genus three built in
15.38 seconds and computed the focused operation in 17.74 seconds; full
middle-degree cohomology was skipped.

Reproduce these runs from the repository root with:

```bash
python benchmarks/symmetric_surfaces.py
python benchmarks/symmetric_surfaces.py 3
```

## Construction and scaling

For positive genus, the base surface is the standard `4g`-gon presentation
triangulated by a fan. It has one vertex, `6g - 3` edges, and `4g - 2`
triangles. A symmetric-power simplex is represented immediately by an
unordered tuple of factor simplices. This produces the quotient without
first allocating the ordered cartesian power.

Genus three is the first case in the unstable range for the classical
cohomology presentation, and it is also where full middle-degree linear
algebra becomes the dominant cost. That makes it the right extended
performance target.

## The p = 5 continuation

The parallel family at the next prime is Sym⁵(Σ_g). The ground-truth case
Sym⁵(S²) = CP⁵ and the positive-genus case Sym⁵(T²) both compute rank-one
P¹: H² → H¹⁰. Compact unordered labels make the 1,797,894-cell torus model
practical, but the predicted 414,092,094 cells at genus two establish a clear
boundary for this presentation. See [the prime-five report](prime-five.md)
for timings and the next search strategy.
