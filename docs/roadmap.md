# Roadmap: compact odd-primary examples

The immediate goal is a self-contained `fastop` package with convincing
odd-primary computations on both familiar and naturally occurring spaces.
Work is organized into independently testable slices.

## 1. Recorded simplicial examples

- Include the matching complex (M_7=\operatorname{Cl}(KG(7,2))).  Over
  \(\mathbf F_3\), it has a nonzero Bockstein
  \(\beta P^0:H^1(M_7)\to H^2(M_7)\).
- Provide simplicial suspension as a public construction.
- Use \(\Sigma^k\mathbf{CP}^3\) to benchmark stable \(P^1\) computations in
  degrees where the operation is no longer a top cup power.

Acceptance criteria: the catalog models have tested Betti numbers and
operation ranks, and the examples are runnable from the public API.

## 2. Face-map and quotient input

- Add a finite `DeltaComplex` represented by dense face-index tables
  `faces[d][simplex][i] = d_i(simplex)`.
- Match Sage's `DeltaComplex.cells()` representation closely enough for a
  dependency-free adapter.
- Refactor chains and the all-target operation evaluator around local face
  restriction rather than global vertex tuples.
- Add strict finite-group quotients whose cell actions commute with every
  face map.

Acceptance criteria: abstract simplicial complexes agree with their
face-map conversion, repeated boundary faces cancel correctly, and a compact
quotient model of \(L^7(3)\) detects both \(\beta P^0\) and \(P^1\).

## 3. Normalized simplicial sets

- Represent a simplex by an underlying nondegenerate cell and a degeneracy
  word.
- Apply a complete composite of face maps before deciding whether a factor
  vanishes in normalized cochains.
- Support finite Sage simplicial sets through an optional adapter.

Acceptance criteria: a small simplicial-set model of
\(\operatorname{Sym}^3(S^2)=\mathbf{CP}^3\) agrees with the catalog
triangulation.

## 4. Unexpected six-manifold example

- Construct
  \(\operatorname{Sym}^3(T^2)=(T^2)^3/S_3\) from a minimal torus simplicial
  set.
- Verify the expected mod-3 cohomology and the nonzero operation
  \(P^1(\eta)=\eta^3\), where \(\eta\) is the divisor class.
- Record timings and model sizes.

Generalized Bott threefolds
\(\mathbf P(\mathcal O\oplus\mathcal O(k))\to\mathbf{CP}^2\) form the next
structured family.  Their cohomology rings can prefilter parameters before a
model is constructed: the natural degree-two class has nonzero cube mod 3
when (3\nmid k).

