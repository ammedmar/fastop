# fastop

`fastop` is a new package for computing cohomology operations on finite simplicial data.

The goal is to provide the user-facing layer around fast operation engines such as `oddp`: start from a simplicial complex, compute cohomology over `F_p`, choose cohomology classes, and evaluate Steenrod operations on those classes.

This repository is intentionally separate from `oddp` while the API and data model settle.

## Development

```bash
python -m pip install -e ".[test]"
python -m pytest
```
