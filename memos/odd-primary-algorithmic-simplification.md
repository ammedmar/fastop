# Odd-primary algorithmic simplification proposal

> **Implementation status (July 2026).** The tuple/dictionary universal
> builder described here is now implemented in
> `src/fastop/_odd_primary_formula.py`. `oddp` is test-only, and the catalog is
> a cache rather than a coverage boundary. The remaining phases below are kept
> as the design record that led to the implementation.

This note studies the current `oddp` implementation as an algorithm rather
than as a transcription of the proof. The goal is to identify the smaller
combinatorial kernel that `fastop` should eventually own.

No code is proposed here. This is the blueprint for the next implementation
round.

## Executive summary

The `oddp` implementation is organized around the maps in the paper:

```text
minimal resolution
  -> periodic resolution
  -> Milnor resolution
  -> tensor chains on simplices
  -> evaluation on a simplicial cochain
```

That structure is valuable for proving the formula, but it is not the shape we
want for computation. Once the prime, operation index, source degree, and
Bockstein flag are fixed, all the resolution maps produce universal
coefficients. The input complex only enters at the last step, where we evaluate
those universal coefficients on faces of target simplices.

So the computational core can be much simpler:

```text
operation index
  -> cached universal omission-pattern coefficients
  -> sparse evaluation on cochain support
  -> cohomology projection
```

The key simplification is to describe everything in terms of omitted vertex
positions inside a target simplex. For a target simplex `T` and a source
simplex `A`, the relevant datum is:

```text
positions of vertices of T that are not in A
```

The paper-facing maps decide the coefficient attached to an ordered tuple of
such omission sets. The evaluator only needs that coefficient.

## What `oddp` is really computing

The public operation in `oddp` takes `(p, s, q, bockstein)`. In `fastop`
language this corresponds to:

```python
p = prime
r = reduced power index
d = source cohomological degree
s = -r
q = -d
target_degree = d + 2*r*(p - 1) + int(bockstein)
omission_size = 2*r*(p - 1) + int(bockstein)
```

The universal part depends only on `(p, r, d, bockstein)`. It does not depend
on the simplicial complex or on the cochain. In `oddp` it is built by:

- `phi_dual`, which chooses the periodic-resolution representative and fixes
  the odd-prime coefficient convention;
- `psi_aux` or `psi_dual`, which expands that representative into Milnor
  shapes;
- `abc` or `abc_dual`, which translates Milnor shapes into tensor factors on
  simplices;
- `alex`, which takes Alexander duals and contributes orientation signs.

The direct evaluator then does exactly this:

```python
for target_simplex in target_faces:
    for tensor, universal_coefficient in universal_terms:
        source_0 = tuple(target_simplex[j] for j in tensor[0])
        ...
        source_(p-1) = tuple(target_simplex[j] for j in tensor[p-1])
        add universal_coefficient * x(source_0) * ... * x(source_(p-1))
```

That is already much simpler than the proof-level construction. We have this
shape in `fastop` now through `UniversalOperation` and the native all-target
evaluator.

For Steenrod operations on a homogeneous cochain, the tensor factors that are
evaluated all have the same dimension: each factor selects a source face of
dimension `d`. In other words, every tensor factor has length `d + 1`. Native
code should treat this as an invariant of the universal operation data and may
validate it at the boundary between universal construction and evaluation.

The support-based evaluator in `oddp` reveals the next simplification. It
starts from source simplices in the cochain support, takes their union as a
possible target simplex, and computes the omitted positions of each source
simplex inside that target. It then looks up the universal coefficient for
that omission pattern.

That is the algorithm we should make native.

## Direct combinatorial model

Define an internal object like this conceptually:

```python
OmissionPattern = tuple[tuple[int, ...], ...]
```

It has `p` factors. Each factor records the positions omitted by one source
face inside a target simplex.

For a homogeneous source cochain of degree `d`, each source simplex has
`d + 1` vertices. A target simplex for `P^r` or `beta P^r` has:

```text
d + 2*r*(p - 1) + int(bockstein) + 1
```

vertices. Therefore each factor omits:

```text
2*r*(p - 1) + int(bockstein)
```

positions.

The native universal data should eventually be a table:

```python
SignatureTable:
    key: OmissionPattern, possibly after a fixed normalization
    value: coefficient in F_p
```

Then evaluation is independent of the resolution machinery:

```python
for ordered source tuple (A_0, ..., A_(p-1)), allowing repeats:
    T = sorted_union(A_0, ..., A_(p-1))

    if T is not a target simplex:
        continue

    D_i = omitted_positions(T, A_i)

    if any(len(D_i) != omission_size):
        continue

    coefficient = signature_table[D_0, ..., D_(p-1)]
    output[T] += coefficient * x[A_0] * ... * x[A_(p-1)]
```

There are two equivalent ways to get the initial table:

1. Build universal tensor terms, as we do now, and convert each tensor factor
   to its complement omission set.
2. Build omission-pattern coefficients directly, bypassing tensor terms.

The first route is safer for the next step because it can use `oddp` as an
oracle. The second route is the desired final shape.

## Why this is simpler

The maps in `oddp` are universal. Recomputing their object hierarchy during
cochain evaluation mixes two separate problems:

- constructing the operation formula;
- evaluating the operation formula on a complex.

In `fastop`, these should be separate.

Universal construction should run once per operation index and be cached.
Simplicial evaluation should run many times on different cochains and should
use only tuples, dictionaries, modular arithmetic, and face lookup tables.

That separation also makes the eventual C++ path clearer. The hot loop is not
the resolution construction. The hot loop is:

- enumerate candidate source tuples;
- form the union target;
- compute omitted positions;
- look up a coefficient;
- multiply cochain values.

That kernel is small enough to port later without dragging the proof-facing
objects into C++.

## Evaluation strategies

### `all_targets`

This is the simplest evaluator and is already native in `fastop`.

It loops over every target simplex and every universal tensor term. Its cost is
roughly:

```text
number of target faces * number of universal terms
```

It is excellent as a correctness baseline because it directly mirrors the
universal tensor formula. It will become too expensive for larger complexes.

### `source_focused`

This should become the main general odd-primary evaluator.

It loops over source simplices in the support of the input cochain and only
creates target candidates that can actually contribute. Its naive cost is:

```text
number of support faces ^ p
```

This is the faithful tensor-evaluation loop: each tensor factor evaluates the
same cochain independently, so repeated source faces must be allowed. In
practice, the target-simplex membership test and omitted-position size test
discard many candidates.

For `p = 3`, `fastop` now has a specialized source-side evaluator named
`source_mod_3`. It works directly with omission patterns. For each universal
pattern, it chooses the pair of tensor factors covering the largest number of
target positions, reconstructs the partial target simplex from those two
source faces, then looks up the third source face by the positions forced by
that partial target. In ordinary covered patterns, this changes the basic
search from support cubed to support squared plus an indexed lookup.

If an exotic pattern has a target position omitted by all three factors, the
source-side reconstruction cannot see that vertex. In that case `source_mod_3`
falls back to the target-omission evaluator for that pattern, preserving
correctness while keeping the fast path for the formulas whose tensor factors
cover the target simplex.

Idea: keep `all_targets` as the default for p=3 until the selector has a
clear sparse-support regime. Native timings show that CP3 and the mod-3 Moore
space are still better served by `all_targets`: their target sets are small
relative to the input cochain support, so the source-pair search performs more
work than the target scan. The specialized `source_mod_3` kernel becomes the
right tool when support is tiny and the target degree has many simplices; in a
synthetic sparse case with 280840 target faces and support size 3, native
`source_mod_3` was orders of magnitude faster than native `all_targets`.

### `target_omissions`

There is another useful variant that `oddp` does not isolate as cleanly:

```python
for target_simplex in target_faces:
    source_faces = source faces of target_simplex that appear in x.support
    enumerate p-tuples from source_faces
```

This can be better when the complex has relatively few target simplices but
the cochain support is broad. It also makes omitted positions cheap because
they are positions in the current target simplex by construction.

In `fastop`, this is the `target_omissions` internal evaluator. It traverses
target simplices, but uses the omission-pattern signature table rather than
the raw tensor terms. This gives a useful comparison point between the raw
`all_targets` tensor baseline and the `source_focused` evaluator.

The native evaluator can eventually choose between:

- `all_targets`;
- `target_omissions`;
- `source_focused`;
- `source_mod_2`;
- `source_mod_3`.

The public API should not expose these names. They are internal strategies.
At present, `source_mod_2` is the extracted version of the original
cohomology-level square support rule. The next mod-2 step is to compare it
against the `fast_sq` implementation and replace the kernel if needed.

## Simplicial lookup simplification

Once the universal formula is fixed, evaluation can use a basic property of
simplicial complexes: vertices determine simplices. We do not need to carry
extra simplex objects through the hot loop. A source face is just the sorted
tuple of vertices selected from a target simplex, and membership or coefficient
lookup is a dictionary/set lookup on that tuple.

This is the same kind of simplification used in the `fast_sq` implementation:
make the universal operation a list of vertex-position patterns, then let the
complex-specific part be fast tuple lookup.

Useful indexes:

- `faces_by_degree[d]`: the set or ordered list of `d`-faces;
- `face_to_index[d]`: tuple of vertices to basis index;
- `cochain_support`: tuple of `(face, coefficient)` pairs with nonzero
  coefficient;
- `target_face_set`: membership table for possible output faces;
- optionally, for `target_omissions`, all `d`-faces contained in a
  target face that also occur in the cochain support.

With these indexes, the evaluator only manipulates vertex tuples and position
tuples. It does not need the resolution objects, and it does not need an
additional simplex identity beyond the ordered vertex tuple.

## Universal data construction

The final native universal builder should not copy the full `oddp` class
hierarchy. It should compute the same coefficients with small functions.

A practical migration sequence:

1. Keep `oddp.chain_operations` as the oracle for universal tensor terms.
2. Convert universal tensor terms into omission-pattern coefficients.
3. Write native evaluators against omission-pattern coefficients.
4. Replace universal construction family by family.

The first native family is already clear:

```text
P^r on a degree 2r class, without Bockstein
```

For this top reduced-power case, the universal tensor has one term:

```python
(
    (0, ..., 2r),
    (2r, ..., 4r),
    ...
)
```

with coefficient `(-1)^r mod p`. This explains why `P^1` on the degree-2
class of `CP^3` is already self-contained over `F_3`.

The next families need derivation and oracle checks:

- `beta P^0`, because it controls the Moore-space example;
- low-degree non-top reduced powers for `p = 3`;
- the same low-degree cases for `p = 5`;
- then a general direct construction of the signature table.

## Sign handling

Signs are the main place to be conservative.

In the support evaluator, `oddp` gets the final coefficient from three pieces:

```text
signature coefficient from psi_aux(phi_dual(...))
coefficient from abc applied to the omission tuple
Alexander-dual sign
```

The simplified model should combine these into one number in the signature
table. The evaluator should not recompute proof-level signs. It should only
look up:

```python
signature_table[omission_pattern]
```

During migration, we can compare table entries against `oddp` for small
operation indices. Once the table is correct, the evaluator is just arithmetic.

This also gives us a good debugging story: if an end-to-end operation fails,
we can separately test whether the universal table is wrong or the simplicial
enumeration is wrong.

## Proposed implementation sequence

### Phase 1: extraction and equivalence

Add an internal conversion:

```python
UniversalOperation -> SignatureTable
```

where each tensor factor is replaced by its complement inside the target
simplex positions.

Then add `source_focused`, using the signature table but still obtaining
the table from the current universal operation source. This gives a new
evaluation path without changing universal construction.

Validation:

- compare `source_focused` output with `all_targets` output on raw cochains;
- compare both with the `oddp` direct evaluator while the bridge remains;
- keep CP3 and Moore examples as end-to-end checks.

### Phase 2: evaluator choice

Make odd-primary operations use an `auto` evaluator.  The conservative default
is `all_targets`; switch to `source_focused` only when the input cochain support
is small enough. Keep `all_targets` as the fallback and as a test oracle.

Useful heuristic inputs:

- number of target faces;
- number of source support faces;
- prime `p`;
- number of universal signature entries.

This phase is purely computational. It should not change mathematical output.

### Phase 3: native universal families (completed)

The general native builder now replaces `oddp.chain_operations` in production.

Suggested order:

1. top reduced powers, already started;
2. `beta P^0`;
3. p=3 low-degree reduced powers needed by examples;
4. p=5 low-degree checks;
5. general signature construction.

Each family should be accepted only after comparing universal data and
cohomology operations against the oracle.

### Phase 4: remove runtime dependency (completed)

`oddp` is now test-only. The package does not import it at runtime.

At that point we can decide whether to keep `oddp` as a development oracle in
optional tests.

### Phase 5: C++ kernel

Only after the Python tuple kernel stabilizes should we move hot loops to C++.

The C++ boundary should be narrow:

```text
faces by dimension
cochain support
signature table
prime
```

It should return the raw target cochain. Cohomology projection, basis objects,
and the public `operation(...)` API should remain Python-level.

## Validation matrix

The tests should cover three levels.

Universal data:

- compare tensor terms from native builders with `oddp.chain_operations`;
- compare omission-pattern tables after tensor-to-signature conversion;
- include primes 3 and 5 when feasible.

Raw cochain evaluation:

- compare `all_targets`, `target_omissions`, and `source_focused`;
- compare native output with `oddp.cochain_operation`;
- use cochains with multiple support terms, not only basis generators.

Cohomology operations:

- `CP^2` over several primes;
- `CP^3` over `F_3`, especially `P^1: H^2 -> H^6`;
- Moore space over `F_3`, especially `beta P^0: H^1 -> H^2`;
- a p=5 example before claiming general odd-primary support.

## Open questions

The main mathematical question is how much of `psi_aux(phi_dual(...))` can be
collapsed into a closed formula for the signature table. The top reduced-power
case suggests there are simple families, but the general case may still need a
small dynamic-programming version of the resolution maps.

The main engineering question is which evaluator wins most often:
`source_focused`, `target_omissions`, or a prime-specific source-side
specialization. We should implement them in Python first, benchmark on the
example spaces, and only then choose the C++ kernel shape.

The main API question is already settled: users should keep calling
`class.operation(r, bockstein=False)`. The prime is determined by the
cohomology object, and evaluator choices remain internal.
