Choose an input model
=====================

``fastop`` has three finite input models. They share ``dimension``, ``cells``
and ``cohomology``, but they deliberately represent different kinds of input.

=======================  ==============================================  ==========================================
Model                    Start from                                      Use it when
=======================  ==============================================  ==========================================
``SimplicialComplex``    facets as tuples of vertex labels                you have an ordinary triangulation
``DeltaComplex``         dense indexed face maps                          identifications create loops or repeated faces
``SimplicialSet``        nondegenerate simplices and their face maps      degeneracies, products, or symmetric powers matter
=======================  ==============================================  ==========================================

Most user-supplied triangulations should begin with ``SimplicialComplex``.
Use the other two only when their extra structure is actually present.

Abstract simplicial complexes
-----------------------------

Pass an iterable of facets. Vertex labels should be integers; the constructor
sorts each facet and supplies all of its faces automatically. The listed
simplices need not already be maximal.

.. code-block:: python

   from fastop import SimplicialComplex

   # The boundary of a triangle: a simplicial circle.
   circle = SimplicialComplex([
       (0, 1),
       (1, 2),
       (0, 2),
   ])

   assert circle.dimension == 1
   assert circle.vertices == (0, 1, 2)
   assert circle.faces(0) == frozenset({(0,), (1,), (2,)})
   assert circle.faces(1) == frozenset({(0, 1), (0, 2), (1, 2)})
   assert circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}

``faces(dimension)`` and ``cells(dimension)`` return the simplices in that
dimension. ``as_delta_complex()`` retains the same realization using only
local face maps; ``suspension(times)`` returns an iterated simplicial
suspension.

Delta-complexes and semi-simplicial input
-----------------------------------------

For quotient models, global vertex labels can be misleading or unavailable.
``DeltaComplex`` instead stores cells by degree and position. For a
``d``-cell with index ``sigma``, ``face_maps[d][sigma][i]`` is the index of
its ``i``-th face among the ``(d - 1)``-cells. Vertices have empty face tuples.

.. code-block:: python

   from fastop import DeltaComplex

   # One vertex and one loop edge. Both faces of the edge are vertex 0.
   circle = DeltaComplex([
       [()],
       [(0, 0)],
   ])

   assert circle.f_vector() == (1, 1)
   assert circle.face(1, 0, 0) == 0
   assert circle.face(1, 0, 1) == 0
   assert circle.cohomology(p=3).betti_numbers() == {0: 1, 1: 1}

Repeated faces are intentional: they express the loop without subdivision.
The constructor checks the semi-simplicial identities by default. Pass
``check=False`` only for trusted, already-validated tables.

Finite simplicial sets
----------------------

``SimplicialSet`` is for models where degeneracies must be retained. Its
direct format uses ``SimplexReference`` objects and is therefore more
specialized than the two formats above. In ordinary use, begin from a
Delta-complex or a built-in small model:

.. code-block:: python

   from fastop import SimplicialSet

   sphere = SimplicialSet.sphere(2)
   projective_space_model = sphere.symmetric_power(3)

   assert sphere.f_vector() == (1, 0, 1)
   assert projective_space_model.cohomology(p=3).operation_rank(2, 1) == 1

``from_delta_complex`` freely adds degeneracies to a Delta-complex model.
Pass normalized face references directly to ``SimplicialSet(...)`` only for
specialized constructions.
``cartesian_product`` and ``symmetric_power`` build normalized simplicial-set
models. Before a large symmetric power, call ``symmetric_power_f_vector`` to
predict its size without constructing it.
