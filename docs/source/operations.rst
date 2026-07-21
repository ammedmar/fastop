Operations
==========

Every finite model exposes ``cohomology(p=...)``. The resulting object can
return Betti numbers, basis classes, cocycle representatives, operation
matrices, and operation ranks.

At the prime two, ``operation(k)`` computes ``Sq^k``. At an odd prime,
``operation(r)`` computes ``P^r`` and
``operation(r, bockstein=True)`` computes ``beta P^r``.

.. code-block:: python

   from fastop import spaces

   H = spaces.lens_space(7, 3).cohomology(p=3)
   assert H.operation_rank(1, 0, bockstein=True) == 1
   assert H.operation_rank(2, 1) == 1

To act on a particular cohomology class, retrieve a basis class by degree:

.. code-block:: python

   H = spaces.complex_projective_space(3).cohomology(p=3)
   u = H.basis(2)[0]
   assert u.operation(1) == H.basis(6)[0]

``operation_rank(degree, r)`` is normally the quickest first question: it
tests the entire map from ``H^degree``. ``operation_matrix`` returns its
columns in the selected basis. At odd primes, set ``bockstein=True`` for
``beta P^r``; at the prime two, the same interface computes ``Sq^r``.

The installed package contains both a small cache of cataloged universal
formulas and a native builder for general formulas. Passing
``formula_source="catalog"`` requires a cached formula, while
``formula_source="computed"`` bypasses the cache and exercises the native
builder. The sibling ``oddp`` project is only a development-time oracle in
optional parity tests; it is not a runtime dependency.
