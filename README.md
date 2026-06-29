# fastop

`fastop` is a new package for computing cohomology operations on spaces presented simplicially.

The goal is to provide user-facing methods that start from a simplicial or
semi-simplicial object, compute cohomology over `F_p`, choose cohomology
classes, and evaluate Steenrod operations on those classes.

## Examples

```python
from fastop import spaces

H = spaces.real_projective_plane().cohomology(p=2)
x = H.basis(1)[0]
assert x.operation(1) == H.basis(2)[0]

CP2 = spaces.complex_projective_plane().cohomology(p=2)
u = CP2.basis(2)[0]
assert u.operation(2) == CP2.basis(4)[0]
```

Odd-primary operations are exposed on the same cohomology classes. These
currently use the neighboring `oddp` package as a cochain-level engine, then
project the resulting cocycle back to cohomology.

```python
H = spaces.complex_projective_space(3).cohomology(p=3)
u = H.basis(2)[0]
assert u.operation(1, algorithm="prime-three") == H.basis(6)[0]

M = spaces.moore_space(3).cohomology(p=3)
a = M.basis(1)[0]
assert a.operation(0, bockstein=True, algorithm="prime-three") == M.basis(2)[0]
```

## Development

```bash
.venv/bin/python -m pytest
```
