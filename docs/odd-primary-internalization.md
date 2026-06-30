# Odd-primary Steenrod operation internalization plan

This memo records how to internalize the useful parts of `oddp` into `fastop`.
It is not an implementation plan that preserves `oddp`'s module structure.
The goal is to keep the mathematics and tests, while replacing the current
verbose object model with smaller kernels that match `fastop`'s cohomology
classes and sparse vectors.

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
evaluators, with `source_focused` as the general odd-primary strategy and
`source_mod_2`/`source_mod_3` as prime-specific optimizations.
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

## Proposed internal shape

A native implementation should live behind the existing `operation(...)` API.
The likely private modules are:

```text
fastop/
  _universal.py
  _cochain_evaluation.py
  _linear_algebra.py
  _oddp_bridge.py
```

The cohomology object should own convention conversion. It has a
`convention` flag:

```python
complex.cohomology(p=3, convention=1)   # cohomological, default
complex.cohomology(p=3, convention=-1)  # homological index signs
```

For the current oddp bridge, the cohomology object converts from the public
operation index to the oddp convention:

```python
operation_degree = convention * r
oddp_s = -operation_degree
oddp_q = -source_degree
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

`_oddp_bridge.py` can temporarily wrap `oddp` for comparison tests while native
kernels are being written.

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

Until native kernels replace the bridge, `oddp` should remain the oracle. The
initial native tests should compare:

- raw cochains from native all-target evaluation against `oddp` direct output;
- raw cochains from native sparse evaluation against `oddp` support output;
- projected cohomology classes through the current `operation(...)` API.

Regression spaces:

- CP2 over `F_3`: `P^1` on the degree-2 class is zero by dimension.
- CP3 over `F_3`: `P^1: H^2 -> H^6` has rank 1.
- Moore space over `F_3`: Bockstein `P^0: H^1 -> H^2` has rank 1.
- A future `p = 5` example before claiming confidence beyond p=3.

## Suggested next implementation slice

Do not port all of `oddp`.

First native slice:

1. Add `_oddp_bridge.py` as the current `oddp` bridge.
2. Move the existing cohomology-level bridge out of `cohomology.py` into that
   bridge module.
3. Keep behavior unchanged.

Second native slice:

1. Implement universal signature extraction for a very small range of cases.
2. Compare it against `oddp`'s universal tensor chains.
3. Use CP3 and Moore only as end-to-end tests.

Third native slice:

1. Implement all-target evaluation.
2. Compare raw output with `oddp` direct output.
3. Then implement `source_focused` and compare raw output with the
   `all_targets` evaluator.

Only after these slices should we consider C++ kernels.
