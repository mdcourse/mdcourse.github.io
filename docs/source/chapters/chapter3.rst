.. _chapter3-label:

Initialize the simulation
=========================

Here, the *InitializeSimulation* class is created. This class is used to
prepare the simulation box and populate the box randomly with atoms.

Start coding
------------

Let us improve the previously created *InitializeSimulation* class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    class InitializeSimulation(Prepare):
        def __init__(self,
                    box_dimensions=[10, 10, 10],  # List - Angstroms
                    seed=None,  # Int
                    initial_positions=None,  # Array - Angstroms
                    thermo_period=None,
                    dumping_period=None,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)
            self.box_dimensions = box_dimensions
            self.dimensions = len(box_dimensions)
            self.seed = seed
            self.initial_positions = initial_positions
            self.thermo_period = thermo_period
            self.dumping_period = dumping_period

.. label:: end_InitializeSimulation_class

Several parameters are provided to the *InitializeSimulation* class:

- The first parameter is the box dimensions, which is a list with a length
  corresponding to the dimension of the system. Each element of the list
  corresponds to a dimension of the box in Ångström in the x, y, and z
  directions, respectively.
- Optionally, a seed can be provided as an integer. If the seed is given
  by the user, the random number generator will always produce the same
  displacements.
- Optionally, initial positions for the atoms can be provided as an array
  of length corresponding to the number of atoms. If *initial_positions* 
  is left equal to *None*, positions will be randomly assigned to the
  atoms (see below).
- Optionally, a thermo period and a dumping_period can be provided to
  control the outputs from the simulation (it will be implemented
  in :ref:`chapter5-label`).

Finally, the *dimensions* of the system are calculated as the length of
*box_dimensions*.


Provide a seed
--------------

A *seed* can be provided to launch the *exact* same simulations (i.e. same atom positions,
same atom velocities) several times in a row, which is useful during debugging. Add
the following *if* condition to the *__init__()* method:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.initial_positions = initial_positions
        if self.seed is not None:
            np.random.seed(self.seed)

.. label:: end_InitializeSimulation_class

If a *seed* is provided, it is passed to the *random.seed()* function of *NumPy*.
If *seed* is left to its default value of *None*, the simulation will be randomized.

Nondimensionalize units
-----------------------

Just like we did in the pervious chapter, let us nondimensionalize the provided
parameters, here the box_dimensions.

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def nondimensionalize_units_1(self):
        """Use LJ prefactors to convert units into non-dimensional."""
        # Normalize box dimensions
        norm_box_dimensions = []
        for L in self.box_dimensions:
            norm_box_dimensions.append(L/self.reference_distance)
        self.box_dimensions = norm_box_dimensions
        if self.initial_positions is not None:
            self.initial_positions = self.initial_positions/self.reference_distance

.. label:: end_InitializeSimulation_class

Define the box
--------------

Let us define a box from the *box_dimensions* list. Add the following method
to the InitializeSimulation class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def define_box(self):
        box_boundaries = np.zeros((self.dimensions, 2))
        for dim, L in zip(range(self.dimensions), self.box_dimensions):
            box_boundaries[dim] = -L/2, L/2
        self.box_boundaries = box_boundaries
        box_size = np.diff(self.box_boundaries).reshape(3)
        box_geometry = np.array([90, 90, 90])
        self.box_size = np.array(box_size.tolist()+box_geometry.tolist())

.. label:: end_InitializeSimulation_class

The *box_boundaries* are calculated from the *box_dimensions*. It corresponds to
the lowest and highest coordinate in all directions. By symmetry, the box is centered
in 0 for all axes. A *box_size* is also defined. It follows the MDAnalysis
conventions: Lx, Ly, Lz, 90, 90, 90, where the last three numbers are angles in
degrees. Values different from *90* for the angles would define a triclinic
(non-orthogonal) boxe, which is not currently supported by the current code.

Populate the box
----------------

Here, the atoms are placed within the simulation box. If initial
positions were not provided (i.e. *initial_positions = None*), atoms
are placed randomly within the box. If initial positions were provided
as an array named *initial_positions*, they are used instead. Note that in that
case, the array number be of size 'number of atoms' x ''number of dimensions.

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def populate_box(self):
        if self.initial_positions is None:
            atoms_positions = np.zeros((self.total_number_atoms,
                                        self.dimensions))
            for dim in np.arange(self.dimensions):
                diff_box = np.diff(self.box_boundaries[dim])
                random_pos = np.random.random(self.total_number_atoms)
                atoms_positions[:, dim] = random_pos*diff_box-diff_box/2
            self.atoms_positions = atoms_positions
        else:
            self.atoms_positions = self.initial_positions

.. label:: end_InitializeSimulation_class

In case initial positions were not provided by the user, and array of size
total_number_atoms x dimensions is created, random positions are defined
using the random function of NumPy.

Here, the newly added atoms are added randomly within the box, without taking care
of avoiding any overlap with existing atoms.

Finally, let us call the methods from the *__init__* class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.initial_positions = initial_positions
        if self.seed is not None:
            np.random.seed(self.seed)
        self.nondimensionalize_units_1()
        self.define_box()
        self.populate_box()

.. label:: end_InitializeSimulation_class

Test the code
-------------

Let us test the *InitializeSimulation* class to make sure that it does what
is expected, i.e. that it does create a simulation box of desired size and
attribute positions to the atoms.

.. label:: start_test_InitializeSimulation_class

.. code-block:: python

    import numpy as np
    from InitializeSimulation import InitializeSimulation

    init = InitializeSimulation(number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 1], # g/mol
        box_dimensions=[20, 20, 20], # A
        seed=48031,
        )

    for d in range(3):
        assert np.round(init.box_size[d],3) == np.round(20/3,3)
    assert np.shape(init.atoms_positions) == (init.total_number_atoms, 3)
    for d in range(3):
        assert init.atoms_positions[0][d] >= init.box_boundaries[0][0]
        assert init.atoms_positions[0][d] <= init.box_boundaries[0][1]

.. label:: end_test_InitializeSimulation_class

The first assert statement ensures that the box size has the expected size along
all 3 dimensions of space, in normalised units.
The second assert ensures that the size of the created *atoms_positions* array
is consistent with the anumber of atoms. The third assert ensures that the created
atoms are located within the box boundaries along all 3 dimensions of space.