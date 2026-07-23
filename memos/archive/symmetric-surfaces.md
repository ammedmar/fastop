# Profile: symmetric products of surfaces

Profiled on 2026-07-23 with Python 3.14.3 and the optional C extension on an
arm64 macOS 26.5.2 development host. Timings are representative single cold
runs, not regression thresholds.

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

## Updated cold profile

The symmetric-power constructor is now lazy. Its sub-millisecond time only
creates a wrapper; cell labels, face restrictions, boundary matrices, and
cohomology bases are materialized when degrees are requested. The useful
real-world measurement is therefore the public focused workflow:

```python
X = spaces.orientable_surface(genus).symmetric_power(3)
H = X.cohomology(p=3)
H.operation_rank(2, 1)
```

`operation_rank` first needs the source and target cohomology. The benchmark
times those prerequisites separately, then times the operation after both
degrees are available.

| Genus | Source H² | Target H⁶ | First P¹ evaluation | Focused cold total | Full cold cohomology |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.0005 s | 0.0013 s | 0.0006 s | 0.0024 s | 0.0052 s |
| 1 | 0.0202 s | 0.0118 s | 0.0119 s | 0.0440 s | 0.1013 s |
| 2 | 0.7923 s | 0.2975 s | 0.2839 s | 1.3737 s | 6.7736 s |
| 3 | 18.3443 s | 1.4779 s | 8.2482 s | 28.0705 s | skipped |

At genus two, source and target cohomology account for 79.4% of the focused
cold query; the operation accounts for 20.7%. At genus three the split is
70.6% versus 29.4%, with source H² alone taking 65.4% of the total. Full
genus-two cohomology reached a cumulative process peak of about 356 MiB. The
focused genus-three run reached about 335 MiB; these high-water figures are
not directly additive or per-object measurements.

Cached operation timings at genera zero through two were respectively
0.0005, 0.0114, and 0.2840 seconds. They are nearly identical to first
evaluation timings. Caching universal formula data therefore has little
effect on this family: traversal and restriction of target cells dominate
the operation.

## Bottlenecks

A deterministic `cProfile` run of the combined genus-two focused and full
workflows attributes the degree-data hot path as follows:

| Work | Cumulative profiled time | Share of profiled degree-data time |
| --- | ---: | ---: |
| Boundary-column construction through normalized face maps | 8.02 s | 61.7% |
| Sparse image/kernel reduction in the C extension | 3.58 s | 27.5% |
| Cohomology quotient-basis construction in the C extension | 0.93 s | 7.2% |

The profiler slows Python-heavy paths, so these are diagnostic shares rather
than replacements for the wall-clock table. It records 277,254 top-level
symmetric-power face calls, expanding to about 1.19 million recursive base
face-reference calls. This is the primary cohomology bottleneck.

Within P¹ evaluation, almost all profiled time is in 10,246 target-cell
restriction calls and their normalization; final projection to cohomology
takes only about 0.004 seconds. The optimization order suggested by this
workload is:

1. Batch or cache normalized faces while constructing boundary columns.
2. Reduce repeated symmetric-cell restriction and normalization during
   cochain evaluation, preferably across all factors of a target cell.
3. Continue improving sparse reduction only after face generation, since the
   accelerated reducer is the second rather than first bottleneck.

The steep genus-two to genus-three increase in source H² (about 23× for a
4.6× cell-count increase) also shows that cell count alone is not a sufficient
cost predictor. Sparse-elimination fill-in and the distribution of cells by
degree matter.

Reproduce these runs from the repository root with:

```bash
python benchmarks/symmetric_surfaces.py
python benchmarks/symmetric_surfaces.py 3 --hot-repeats 0
python -m cProfile -s cumulative benchmarks/symmetric_surfaces.py 2 --hot-repeats 0
```

## Construction and scaling

For positive genus, the base surface is the standard `4g`-gon presentation
triangulated by a fan. It has one vertex, `6g - 3` edges, and `4g - 2`
triangles. A symmetric-power simplex is represented immediately by an
unordered tuple of factor simplices. This produces the quotient without
first allocating the ordered cartesian power.

Genus three is the first case in the unstable range for the classical
cohomology presentation. It is also the right extended performance target:
the focused query still completes, while full middle-degree cohomology is
large enough that it should remain opt-in.

## The p = 5 continuation

The parallel family at the next prime is Sym⁵(Σ_g). The ground-truth case
Sym⁵(S²) = CP⁵ and the positive-genus case Sym⁵(T²) both compute rank-one
P¹: H² → H¹⁰. Compact unordered labels make the 1,797,894-cell torus model
practical, but the predicted 414,092,094 cells at genus two establish a clear
boundary for this presentation. See [the prime-five report](prime-five.md)
for timings and the next search strategy.
