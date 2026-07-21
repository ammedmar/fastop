# Odd-primary Steenrod operation internalization

This memo records how the useful formula-building path from `oddp` was
internalized into `fastop`. The work is complete: installed `fastop` does not
import `oddp` at runtime. It keeps the mathematics and oracle comparisons while
replacing the resolution object hierarchy with tuple/dictionary kernels that
match `fastop`'s cohomology classes and sparse vectors.

## Terminology

The public language in `fastop` should be mathematical:

- `operation(r)` over `F_2` means `Sq^r`.
- `operation(r)` over odd `F_p` means the reduced power `P^r`.
- `operation(r, bockstein=True)` over odd `F_p` means the Bockstein followed
  by `P^r`.

The names `support` and `prime-three` should not become preferred `fastop`
terminology. They are names of two evaluation strategies in `oddp`, not
mathematical operations. They may remain as compatibility aliases only.

Better internal names:

- `all_targets`: evaluate the universal tensor formula on every target simplex.
  This corresponds to `oddp`'s `direct` algorithm.
- `target_omissions`: evaluate every target simplex using omission-pattern
  signatures instead of raw tensor terms.
- `source_focused`: build candidate target simplices from source faces in the
  input cochain. This corresponds to `oddp`'s `support` algorithm.
- `source_mod_2`: the mod-2 source-focused square evaluator, matching the
  `fast_sq` direction.
- `source_mod_3`: a p=3 specialization of `source_focused`. This corresponds
  to `oddp`'s `prime-three` algorithm.

The first two are target-side evaluators. The last three are source-side
evaluators.  The public odd-primary default should use an `auto` policy:
prefer `all_targets` as the conservative target-side baseline, and choose
`source_focused` only when the input support is small enough to make source
tuple enumeration cheaper.  `source_mod_2`/`source_mod_3` remain
prime-specific source-side optimizations.
The current `source_mod_2` implementation is an extraction of the existing
mod-2 square support rule; it still needs a direct comparison with the
`fast_sq` implementation before we claim it is the final mod-2 kernel.

## What `oddp` currently does

The public operation entry point is:

```python
Steenrod.cochain_operation(complex, cochain, p, s, q, bockstein=False)
```

Its conventions differ from the `fastop` API:

- `q` is the negative cohomological degree.
- `s` is the negative reduced-power index.

For a `fastop` class in degree `d`, the bridge is:

```python
q = -d
s = -r
target_degree = d + 2*r*(p - 1) + int(bockstein)
```

The core pipeline in `oddp` is:

```text
operation parameters
  -> minimal resolution generator
  -> phi_dual
  -> psi_aux or psi_dual
  -> abc or abc_dual
  -> universal tensor/signature data
  -> simplicial evaluation on the cochain
```

The last step is the part that touches the input complex. Everything before it
is universal data depending only on `(p, r, source_degree, bockstein)`.

## What to keep

Keep the formulas and tests, not the class hierarchy.

Useful pieces:

- index conversion in `powers._operation_parameters`;
- `phi_dual`, especially its odd-prime coefficient convention;
- `psi_aux`, which gives the compact "signature" data used by sparse
  evaluation;
- `abc` and `abc_dual`, which translate between the tensor shape and the
  Milnor-resolution shape;
- the CP3 and Moore-space regression examples.

Pieces to avoid copying directly:

- the public `Steenrod` wrapper;
- the broad `ParametricCounter` arithmetic model;
- the resolution object hierarchy as public or semi-public API;
- the algorithm names `support` and `prime-three`, except as compatibility
  aliases.

`fastop` already has sparse vectors over `F_p`; native code should use plain
dicts with explicit modular arithmetic.

## Implemented internal shape

A native implementation should live behind the existing `operation(...)` API.
The likely private modules are:

```text
fastop/
  _universal.py
  _odd_primary_formula.py
  _cochain_evaluation.py
  _linear_algebra.py
tests/
  _oddp_oracle.py
```

The cohomology object should own convention conversion. It has a
`convention` flag:

```python
complex.cohomology(p=3, convention=1)   # cohomological, default
complex.cohomology(p=3, convention=-1)  # homological index signs
```

The production path uses the public `fastop` convention directly:

```python
operation_degree = convention * r
target_degree = source_degree + 2*operation_degree*(p - 1) + int(bockstein)
missing_vertices_per_factor = 2*operation_degree*(p - 1) + int(bockstein)
```

`universal.py` should produce cached universal data:

```python
UniversalOperation(
    p: int,
    r: int,
    source_degree: int,
    bockstein: bool,
    target_degree: int,
    signatures: dict[tuple[tuple[int, ...], ...], int],
)
```

The cache key should be:

```python
(p, r, source_degree, bockstein)
```

`_cochain_evaluation.py` should contain the cochain-level evaluators. The
source-focused evaluator takes:

```python
faces_by_degree
cochain_support
universal_operation
```

and returns a sparse cochain in the target degree. It should not know about
cohomology projection.

`tests/_oddp_oracle.py` wraps `oddp` only for optional comparison tests.
Production modules contain no `oddp` adapter and never import it.

## Evaluation strategies

### `all_targets`

For every target simplex in the target degree, try every tensor summand in the
universal operation. This mirrors `oddp`'s `direct` algorithm.

This is useful as a simple reference implementation but should not be the main
path for large complexes.

### `source_focused`

Choose `p` source simplices from the input cochain support. Their union is a
possible target simplex. If that union is actually present in the complex and
the omitted vertex positions have the required sizes, evaluate the universal
signature on that configuration.

This mirrors `oddp`'s `support` algorithm.

### `source_mod_3`

When `p = 3`, the source tuple has exactly three factors. `oddp` uses a
candidate list to reduce triple enumeration. This is a good optimization idea,
but the name `prime-three` should not survive.

In `fastop`, this should be an internal branch of the source-side evaluator:

```python
if p == 3:
    use_source_mod_3(...)
else:
    use_source_focused(...)
```

The caller should still ask for `operation(r)`, not for an algorithm named
after the prime.

## Test oracle

`oddp` remains a valuable independent oracle. Optional integration tests
compare:

- raw cochains from native all-target evaluation against `oddp` direct output;
- raw cochains from native sparse evaluation against `oddp` support output;
- projected cohomology classes through the current `operation(...)` API.

Regression spaces:

- CP2 over `F_3`: `P^1` on the degree-2 class is zero by dimension.
- CP3 over `F_3`: `P^1: H^2 -> H^6` has rank 1.
- Moore space over `F_3`: Bockstein `P^0: H^1 -> H^2` has rank 1.
- the symmetric-product and lens-space examples at `p = 5`.

## Completed implementation slices

The implementation intentionally did not port all of `oddp`.

Completed slices:

1. Isolate the `oddp` bridge as a test-only oracle.
2. Implement and validate cataloged formulas and the top reduced-power family.
3. Implement all-target, omission-table, and source-focused evaluators.
4. Absorb the general `phi_dual -> psi_dual -> abc_dual -> alex` formula path
   as private tuple/dictionary routines.
5. Compare cataloged and uncataloged formulas term-by-term at primes 3 and 5.
6. Retain the optional C extension for the stable evaluation kernels.

The adapted formula code carries the upstream MIT notice in
`THIRD_PARTY_NOTICES.md`.
