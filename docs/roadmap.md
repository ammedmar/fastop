# Roadmap: compact odd-primary examples

The immediate showcase is complete: `fastop` now computes odd-primary
operations on abstract simplicial complexes, face-map Delta-complexes, and
finite simplicial sets, including compact quotients and symmetric products.

## 1. Recorded simplicial examples — complete

- The matching complex M₇ = Cl(KG(7,2)) has mod-3 Betti numbers
  `(1, 1, 21)` in degrees 0 through 2, and
  βP⁰: H¹(M₇) → H²(M₇) has rank 1.
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
- `symmetric_product_of_surface(genus, power=3)` exposes the complete family.
- Symmetric powers are built directly from unordered simplex tuples rather
  than by constructing the ordered product and then quotienting it.
- Sym³(Σ₂) has 41,478 cells, measured mod-3 Betti numbers
  `(1, 4, 7, 8, 7, 4, 1)`, and rank-one P¹: H² → H⁶.
- The first unstable case Sym³(Σ₃) has 189,626 cells and also computes a
  rank-one P¹: H² → H⁶. Full middle-degree cohomology is deliberately left
  out of the default regression suite because it is much more expensive.

## Next prime: p = 5

The next search will keep the surface family and increase the symmetric
power. Sym⁵(S²) = CP⁵ provides ground truth for
P¹: H² → H¹⁰, where the unstable identity gives P¹(x) = x⁵. Higher-genus
Sym⁵(Σ_g) will then be the search family. Before constructing those models,
the direct unordered algorithm will need a multiplicity-based enumerator so
it does not iterate over every five-element tuple.

Generalized Bott threefolds remain a possible later family, but are postponed
until the symmetric-product program at primes 3 and 5 has been explored.
