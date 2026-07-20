# Roadmap: compact odd-primary examples

The immediate showcase is complete: `fastop` now computes odd-primary
operations on abstract simplicial complexes, face-map Delta-complexes, and
finite simplicial sets, including compact quotients and symmetric products.

## 1. Recorded simplicial examples — complete

- Simplicial suspension is a public construction.
- On ΣᵏCP³, P¹ remains nonzero after suspension; in particular,
  P¹: H⁴(Σ²CP³) → H⁸(Σ²CP³) has rank 1, outside the unstable cup-power
  degree.

## 2. Face-map and quotient input — complete

- `DeltaComplex` stores dense face-index tables compatible with Sage's
  `DeltaComplex.cells()` representation.
- Chains and operation evaluation use local face restriction, so loops and
  repeated faces are supported without introducing global vertex labels.
- Strict finite group actions can be quotiented when they commute with every
  face map.
- The catalog model of L⁷(3) has 800 cells, mod-3 Betti number 1 in every
  degree from 0 through 7, and detects both rank-one operations
  βP⁰: H¹ → H² and P¹: H² → H⁶.

## 3. Normalized simplicial sets — complete

- A simplex is represented by an underlying nondegenerate cell and its full
  order-preserving surjection; degeneracies are retained until a complete
  face restriction has been evaluated.
- Finite Sage simplicial sets can be copied through an optional adapter.
- Cartesian products, factor permutations, strict quotients, and symmetric
  powers are public constructions.
- The 84-cell model Sym³(S²) = CP³ has the expected mod-3 cohomology and a
  rank-one P¹: H² → H⁶, validating the construction against the catalog
  triangulation.

## 4. Unexpected six-manifold example — complete

- Sym³(T²) = (T²)³/S₃ is constructed from the minimal simplicial torus.
- Its model has 1,638 nondegenerate simplices and mod-3 Betti numbers
  `(1, 2, 2, 2, 2, 2, 1)`.
- The computation finds rank-one P¹: H² → H⁶. This gives a compact,
  naturally occurring six-manifold example beyond a projective-space
  catalog model.
- Reproducible construction, cohomology, and operation timings live in
  `benchmarks/showcase.py`.

## 5. Symmetric cubes across surface genus — complete

- Compact one-vertex models of the closed orientable surface Σ_g have
  f-vector `(1, 6g - 3, 4g - 2)` for positive genus.
- `symmetric_product_of_curve(genus, power=3)` exposes the complete family;
  the older surface-named constructor remains as an alias.
- Symmetric powers are built directly from unordered simplex tuples rather
  than by constructing the ordered product and then quotienting it.
- Sym³(Σ₂) has 41,478 cells, measured mod-3 Betti numbers
  `(1, 4, 7, 8, 7, 4, 1)`, and rank-one P¹: H² → H⁶.
- The first unstable case Sym³(Σ₃) has 189,626 cells and also computes a
  rank-one P¹: H² → H⁶. Full middle-degree cohomology is deliberately left
  out of the default regression suite because it is much more expensive.

## 6. Prime-five baseline — complete

- Sym⁵(S²) = CP⁵ has 22,010 cells, the expected even-degree mod-5
  cohomology, and rank-one P¹: H² → H¹⁰.
- Sym⁵(T²) has 1,797,894 cells. The focused computation measures
  `b₂ = 2`, `b₁₀ = 1`, and rank-one P¹: H² → H¹⁰.
- L¹¹(5) gives an independent face-map validation with 354,312 cells. It
  detects both rank-one βP⁰: H¹ → H² and rank-one P¹: H² → H¹⁰.
- Symmetric powers now store compact unordered labels and compute faces on
  demand. Cohomology boundary data is also built only in requested degrees.
- A public `symmetric_power_f_vector(power)` preflight predicts model size
  without constructing the model.

## Next prime-five search

Sym⁵(Σ₂) would contain 414,092,094 cells in the present model, so increasing
surface genus is no longer the efficient next move. The next search should
screen compact ten-dimensional face-map models from other geometric or
combinatorial families, using cohomology rings or cell counts before any
cochain operation is attempted.

Generalized Bott manifolds remain a possible later family, but are postponed
while the prime-five search focuses on examples less directly predicted by a
projective-space-type degree-two generator.

## Package stabilization

The example-search phase now feeds a small public API:

- `DeltaComplex` is the compact face-map input for semi-simplicial models.
- `FiniteGroupAction` packages graded permutation generators, can build them
  from cell maps, and exposes order and freeness checks.
- `orientable_surface`, `symmetric_product_of_curve`, and `lens_space` are the
  user-facing constructors for the two showcase families.
- expensive Sym⁵(T²) and L¹¹(5) checks live in the opt-in `large` test tier;
  smaller ground truths stay in the routine suite.
- the native extension is an optional accelerator, so a compiler failure does
  not prevent installation of the pure-Python package.
- GitHub Actions verifies the routine suite on Python 3.10 through 3.14,
  builds the generated documentation, and validates source and wheel artifacts.

The remaining work before a first alpha release is release administration:
choose the alpha version, build and inspect distribution artifacts, publish to
a test index, and exercise the installation on the supported Python versions.
