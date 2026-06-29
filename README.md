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
assert x.sq(1) == H.basis(2)[0]

CP2 = spaces.complex_projective_plane().cohomology(p=2)
u = CP2.basis(2)[0]
assert u.sq(2) == CP2.basis(4)[0]
```

## Development

```bash
.venv/bin/python -m pytest
```
