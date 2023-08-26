Initialize simulation
=====================

The *InitializeSimulation* class deals with creating the simulation box,
populating it with a chosen number of atoms with initial velocity. The
*InitializeSimulation* class also takes care of converting the inputs units 
into non-dimensionalized.

Initialize box
--------------

The *initialize_box* functions takes the inputs values for :math:`Lx` (and optionally 
:math:`Ly` and :math:`Lz`) to create a cuboid box of dimensions :math:`Lx \times Ly \times Lz`,
if :math:`Ly` and :math:`Lz` are provided, or :math:`Lx \times Lx \times Lx`, 
if :math:`Ly` and :math:`Lz` are not provided. For system with *dimensions* lower than 3,
2D or even 1D simulation boxes are created.

.. code-block:: python

    def initialize_box(self):
        """Define box boundaries based on Lx, Ly, and Lz.

        If Ly or Lz or both are None, then Lx is used instead"""
        box_boundaries = np.zeros((self.dimensions, 2))
        for dim, L in zip(range(self.dimensions), [self.Lx, self.Ly, self.Lz]):
            if L is not None:
                box_boundaries[dim] = -L/2, L/2
            else:
                box_boundaries[dim] = -self.Lx/2, self.Lx/2
        box_boundaries = self.nondimensionalise_units(box_boundaries, "distance")
        self.box_boundaries = box_boundaries

Note that choosing a dimension different fom 3 may cause problem in other parts of the code, to be fixed.