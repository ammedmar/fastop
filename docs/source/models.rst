Models and constructions
========================

``fastop`` uses three interchangeable finite models:

* ``SimplicialComplex`` for ordinary triangulations;
* ``DeltaComplex`` for dense face maps, loops, repeated faces, and compact
  quotient models; and
* ``SimplicialSet`` for normalized nondegenerate cells, products, and
  symmetric powers.

The ``spaces`` module provides the two showcase families:

.. code-block:: python

   from fastop import spaces

   curve_power = spaces.symmetric_product_of_curve(genus=1, power=3)
   lens = spaces.lens_space(dimension=11, order=5)

For a symmetric power, preflight its size before construction:

.. code-block:: python

   surface = spaces.orientable_surface(genus=2)
   cells = surface.symmetric_power_f_vector(power=5)

Finite actions are strict degreewise cell permutations. The public action
object creates them directly from cell maps and a quotient validates both face
compatibility and, when requested, freeness:

.. code-block:: python

   from fastop import DeltaComplex, FiniteGroupAction

   polygon = DeltaComplex([
       [() for _ in range(5)],
       [((i + 1) % 5, i) for i in range(5)],
   ])
   rotation = FiniteGroupAction.from_cell_maps(
       polygon, lambda degree, cell: (cell + 1) % 5
   )
   circle = polygon.quotient(rotation, require_free=True)
