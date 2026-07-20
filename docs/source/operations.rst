Operations
==========

Every finite model exposes ``cohomology(p=...)``. The resulting object can
return Betti numbers, basis classes, operation matrices, and operation ranks.

At the prime two, ``operation(k)`` computes ``Sq^k``. At an odd prime,
``operation(r)`` computes ``P^r`` and
``operation(r, bockstein=True)`` computes ``beta P^r``.

.. code-block:: python

   from fastop import spaces

   H = spaces.lens_space(7, 3).cohomology(p=3)
   assert H.operation_rank(1, 0, bockstein=True) == 1
   assert H.operation_rank(2, 1) == 1

The installed package contains its cataloged universal formulas and the native
top reduced-power formula used in the prime-five examples. The sibling
``oddp`` project remains a development fallback for uncataloged formulas; it
is not a runtime dependency.
