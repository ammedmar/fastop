# Sage cohomology representatives and Steenrod-square warm-up

This memo records how Sage produces explicit cohomology representatives and
how its existing mod-2 Steenrod-square implementation uses them. Sage is a
reference implementation here; `fastop` should not depend on Sage at runtime.

Primary Sage entry points studied:

- `sage.homology.homology_vector_space_with_basis`
- `sage.homology.algebraic_topological_model`
- `sage.homology.chains`
- `sage.topology.cell_complex`
- `sage.topology.simplicial_complex`
- `sage.topology.delta_complex`
- `sage.topology.simplicial_set`

Upstream documentation:
<https://doc.sagemath.org/html/en/reference/homology/sage/homology/homology_vector_space_with_basis.html>

## Executive summary

Sage stores a cohomology class as a coordinate vector in a finite graded vector
space with basis. A basis class is indexed by `(degree, index)` and printed as
`h^{degree,index}`. The explicit cocycle representative is not the stored
identity of the class; it is recovered from the algebraic topological model by
calling `to_cycle()`.

For cohomology, `to_cycle()` returns a Sage `Cochains.Element`: a sparse linear
combination of dual cells written as terms like `\chi_(0, 1, 2)`. Iterating over
that object gives `(cell, coefficient)` pairs. For a simplicial complex, `cell`
is a simplex object whose vertices can be converted to a tuple. For a
simplicial set, `cell` is a nondegenerate simplex object. For a Delta complex,
`cell` is a tagged cell `(index, face_tuple)` because different cells may have
the same face data.

The key conversion target for `fastop` is therefore:

```text
CohomologyClass {
  degree: int,
  coordinates: dict[(degree, basis_index), F_p],
  representative: SparseCochain,
  source_basis: ChainBasis,
  model_metadata: ...
}

SparseCochain {
  degree: int,
  p: int,
  terms: dict[CellKey, int mod p]
}
```

For pure simplicial complexes, `CellKey` can be `tuple[int, ...]`. For
semi-simplicial or Delta-style inputs, `CellKey` should preserve a stable cell
id rather than only a vertex tuple.

## Sage's homology/cohomology object model

The generic user-facing methods are:

- `X.homology_with_basis(base_ring)`
- `X.homology_with_basis(base_ring, cohomology=True)`
- `X.cohomology_ring(base_ring)`

For field coefficients, Sage constructs a `HomologyVectorSpaceWithBasis` or a
`CohomologyRing`. Over `GF(2)`, Sage dispatches to mod-2 subclasses that add a
Steenrod algebra module structure and the method `Sq(i)`.

The constructor calls:

```text
phi, M = cell_complex.algebraic_topological_model(base_ring)
```

For cohomology, Sage dualizes the contraction:

```text
phi = phi.dual()
```

The resulting object stores:

- `_contraction`: the chain contraction `phi`;
- `_complex`: the original cell complex;
- `_cohomology`: whether this is cohomology;
- `_graded_indices`: degree-by-degree ranges of basis indices.

Basis elements are monomials in a `CombinatorialFreeModule`, indexed by
`(degree, index)`. A linear combination of basis elements is still a class in
the abstract cohomology vector space, not itself a cocycle.

## How representatives are produced

The representative method is:

```text
class HomologyVectorSpaceWithBasis.Element:
    def to_cycle(self):
        if not self.is_homogeneous():
            raise ValueError
        return sum(c * parent._to_cycle_on_basis(i) for i, c in self)
```

For one basis element `(degree, index)`, Sage does:

```text
vec = contraction.iota().in_degree(degree).column(index)
chains = complex.n_chains(degree, base_ring, cochains=_cohomology)
return chains.from_vector(vec)
```

In homology this is a chain representative. In cohomology this is a cocycle
representative in the dual free module of cochains. The cocycle is sparse:
iteration gives only nonzero terms.

The representative is chosen by the chain contraction. It is deterministic for
a fixed Sage version, fixed cell ordering, and fixed linear algebra behavior,
but it is not a canonical mathematical representative. Changing cell order,
model, or reduction choices can change the cocycle while preserving its
cohomology class.

## What counts as a representative

For `fastop`, the useful representative is the actual cochain returned by
`to_cycle()`, not only the abstract class coordinates. It supports:

- iteration as `(cell, coefficient)` pairs;
- `to_vector()` for coordinates in Sage's cochain basis;
- `is_cocycle()` to verify the coboundary is zero;
- `eval(chain)` for pairing with a chain;
- conversion to a dictionary by iterating over terms.

For a simplicial complex over `GF(p)`, a direct conversion is:

```text
{tuple(cell): int(coefficient) % p for cell, coefficient in class.to_cycle()}
```

This dictionary is close to the current `oddp` style input, provided that the
simplex orientation convention agrees with the operation code. For `p = 2`,
orientation signs disappear. For odd primes, `fastop` must keep a precise
orientation convention.

For a Delta complex, do not convert to `dict[tuple[int, ...], int]` unless the
tuple includes the cell id. Sage's own basis keys look like `(idx, face_tuple)`.
Flattening to `face_tuple` can identify distinct cells incorrectly.

For a simplicial set, Sage's normalized chain basis is the set of
nondegenerate simplices. A representative should preserve simplex identity and
its face/degeneracy structure, not only a tuple of vertices.

## Model-by-model notes

### SimplicialComplex

This is the cleanest starting model for `fastop`.

- `n_chains(n, base_ring, cochains=True)` uses sorted `n`-simplices as the
  cochain basis.
- `face(simplex, i)` deletes the `i`th vertex.
- `alexander_whitney(simplex, dim_left)` returns the standard front/back split
  `(1, left, right)`.
- Sage's `Sq(i)` implementation is available for cohomology over `GF(2)`.

For this model, explicit cocycles can naturally be converted to dictionaries
keyed by oriented vertex tuples.

### DeltaComplex

Sage can compute cohomology with basis and explicit cocycle representatives for
Delta complexes. The chain module deliberately tags cells by index:

```text
(idx, face_tuple)
```

This avoids ambiguity when several cells have identical face tuples. Sage also
has an Alexander-Whitney split for Delta-complex cells, using iterated top and
bottom faces.

However, Sage's `Sq(i)` explicitly raises `NotImplementedError` for Delta
complexes. A future `fastop` implementation can still support Delta or
semi-simplicial models if it has its own diagonal/higher diagonal formulas and
a stable indexed-cell basis.

### SimplicialSet

Sage uses normalized chains: `n_chains(n)` is the free module on
nondegenerate `n`-simplices.

The face maps may produce degenerate simplices. Sage's Alexander-Whitney method
returns coefficient `0` if either factor is degenerate and coefficient `1`
otherwise, making it easy to ignore degenerate terms in cup-product style
computations. The existing `Sq(i)` implementation supports simplicial sets and
checks nondegeneracy while applying face maps.

This model is important because Sage's own examples show it can be much faster
than large triangulated simplicial complexes, especially for spaces like real
projective spaces.

## Basis choices and reproducibility

Sage's representatives come from an algebraic topological model: a chain
complex `C`, a smaller chain complex `M` with zero differential, chain maps
`pi: C -> M` and `iota: M -> C`, and a chain homotopy `phi`. Since `M` has zero
differential, its generators are the chosen homology basis; `iota` lifts those
generators to cycles. After dualizing, `iota` lifts cohomology basis elements
to cocycles.

The model-construction algorithm scans cells in degree order and uses the
complex's sorted cell ordering. When it needs a pivot, it takes a nonzero
available pivot from the current vector. This is deterministic if all upstream
orderings are deterministic.

For `fastop`, reproducibility should be made explicit:

- store the ordered chain basis in each degree;
- define one deterministic pivot convention;
- record the coefficient field;
- record the geometric model and any normalization convention;
- keep a stable distinction between class coordinates and cochain
  representatives.

## Sage's existing mod-2 Steenrod-square algorithm

Sage implements `Sq(i)` on `CohomologyRing_mod2.Element`.

Important limitations:

- coefficients must be `GF(2)`;
- the main implementation supports `SimplicialComplex` and
  `SimplicialSet`;
- cubical complexes are converted to simplicial complexes;
- Delta complexes are not supported.

For a homogeneous class `x` of degree `j`, Sage handles:

- `Sq^0(x) = x`;
- `Sq^j(x) = x cup x`;
- `Sq^i(x) = 0` if the target cohomology group is zero or `i > j`;
- otherwise, a Gonzalez-Diaz--Real face-map formula.

The nontrivial case works like this:

1. Lift the class to a cocycle: `cycle = x.to_cycle()`.
2. Let `m = j + i` and `n = j - i`.
3. Build the finite list of face-index sequences appearing in the GDR formula.
4. For each homology basis class in degree `m`, lift it to a cycle.
5. For each cell in that cycle and each face-index sequence, form two
   `j`-simplices by applying face maps.
6. Evaluate the original cocycle on both resulting `j`-simplices.
7. Sum the products in `GF(2)`.
8. Use the homology/cohomology dual bases to convert those evaluations into
   coordinates of the output cohomology class.

This is a useful warm-up architecture for `fastop`: compute an operation on
representatives, then project the resulting cocycle/class back to the chosen
cohomology basis by evaluation against homology representatives.

## Proposed `fastop` data model

`fastop` should expose class-level operations while preserving representative
data:

```text
CellBasis {
  degree: int,
  cells: tuple[CellKey, ...],
  orientation: OrientationConvention,
}

Cochain {
  degree: int,
  p: int,
  basis: CellBasis,
  terms: dict[CellKey, int],
}

CohomologyBasis {
  p: int,
  classes_by_degree: dict[int, tuple[CohomologyClass, ...]],
  homology_cycles_by_degree: dict[int, tuple[Chain, ...]],
}

CohomologyClass {
  degree: int,
  p: int,
  coordinates: dict[int, int],
  representative: Cochain,
}
```

`CellKey` should have model-specific variants:

- simplicial complex: oriented vertex tuple;
- semi-simplicial/Delta: stable cell id plus ordered face data;
- simplicial set: stable nondegenerate simplex id plus face/degeneracy
  operators.

The old `dict[tuple[int, ...], int]` format should remain an adapter for
simplicial complexes, not the core representation for all models.

## Proposed path to Steenrod-operation computation

1. Start with finite simplicial complexes over `GF(2)`.
2. Build ordered chain and cochain bases in every degree.
3. Compute a deterministic contraction or equivalent homology reduction.
4. Store cohomology classes with explicit cocycle representatives.
5. Implement a Sage-compatible `Sq^i` warm-up using the same high-level flow:
   representative -> face formula -> evaluation against homology basis ->
   output class coordinates.
6. Compare against Sage on small examples: spheres, tori, projective spaces,
   suspensions, and Klein bottle variants.
7. Generalize the core representation from vertex tuples to stable cell keys.
8. Add Delta/semi-simplicial support once the cell-key and face-map interface is
   solid.
9. Connect odd-prime operations to the existing `oddp` algorithms through a
   conversion layer from `CohomologyClass.representative` to the operation
   input format.

The immediate implementation target should not be the full cohomology ring.
The minimum useful target is: show the user a basis of cohomology classes,
allow them to select a class, expose its explicit representative, and compute
`Sq^i` on that selected class.

## Reference value for `fastop`

Sage is most useful as a reference for:

- the distinction between abstract class coordinates and cocycle
  representatives;
- the use of chain contractions to choose representatives;
- deterministic but noncanonical basis selection;
- conversion from operation-level evaluations back to cohomology coordinates;
- special handling needed for simplicial sets and Delta complexes.

Sage is less useful as a direct runtime model for `fastop` because its objects
are broad, generic, and optimized for the Sage ecosystem rather than a small
operation-focused library.
