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

Set a random seed
-----------------

Providing a *seed* allows for reproducible simulations. When a seed is
specified, the simulation will produce the exact same results each time it
is run, including identical atom positions and velocities. This can be
particularly useful for debugging purposes.

Add the following *if* condition to the *__init__()* method:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.initial_positions = initial_positions
        if self.seed is not None:
            np.random.seed(self.seed)

.. label:: end_InitializeSimulation_class

If a *seed* is provided, it is used to initialize the *random.seed()* function
from *NumPy*. If *seed* is set to its default value of *None*, the simulation
will proceed with randomization.

Nondimensionalize units
-----------------------

Just like we did in :ref:`chapter2-label`, let us nondimensionalize the provided
parameters, here the *box_dimensions* as well as the *initial_positions*:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def nondimensionalize_units_1(self):
        """Use LJ prefactors to convert units into non-dimensional."""
        # Normalize box dimensions
        box_dimensions = []
        for L in self.box_dimensions:
            box_dimensions.append(L/self.reference_distance)
        self.box_dimensions = box_dimensions # errase the previously defined box_dimensions
        # Normalize the box dimensions
        if self.initial_positions is not None:
            self.initial_positions = self.initial_positions/self.reference_distance

.. label:: end_InitializeSimulation_class

Let us call the *nondimensionalize_units_1* method from the *__init__* class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        if self.seed is not None:
            np.random.seed(self.seed)
        self.nondimensionalize_units_1()

.. label:: end_InitializeSimulation_class

Define the box
--------------

Let us define the simulation box using the values from *box_dimensions*. Add the following
method to the *InitializeSimulation* class:

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

The *box_boundaries* are calculated from the *box_dimensions*. They
represent the lowest and highest coordinates in all directions. By symmetry,
the box is centered at 0 along all axes. A *box_size* is also defined,
following the MDAnalysis conventions: Lx, Ly, Lz, 90, 90, 90, where the
last three numbers are angles in degrees. Values different from *90* for
the angles would define a triclinic (non-orthogonal) box, which is not
currently supported by the existing code.

Let us call the *define_box* method from the *__init__* class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.nondimensionalize_units_1()
        self.define_box()

.. label:: end_InitializeSimulation_class

Populate the box
----------------

Here, the atoms are placed within the simulation box. If initial
positions were not provided (i.e., *initial_positions = None*), atoms
are placed randomly within the box. If *initial_positions* was provided
as an array, the provided positions are used instead. Note that, in this
case, the array must be of size 'number of atoms' times 'number of dimensions'.

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

In case *initial_positions is None*, and array is first created. Then, random
positions that are constrained within the box boundaries are defined using the
random function of NumPy. Note that, here, the newly added atoms are added
randomly within the box, without taking care of avoiding overlaps with
existing atoms. Overlaps will be dealt with using energy minimization,
see :ref:`chapter4-label`.

Let us call the *populate_box* method from the *__init__* class:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.define_box()
        self.populate_box()

.. label:: end_InitializeSimulation_class
        
Test the code
-------------

Let us test the *InitializeSimulation* class to make sure that it does what
is expected, i.e. that it does create atoms that are located within the box
boundaries along all 3 dimensions of space:

.. label:: start_test_3a_class

.. code-block:: python

    import numpy as np
    from InitializeSimulation import InitializeSimulation

    init = InitializeSimulation(number_atoms=[2, 3],
        epsilon=[0.2, 0.4], # kcal/mol
        sigma=[3, 4], # A
        atom_mass=[10, 20], # g/mol
        box_dimensions=[20, 20, 20], # A
        )

    def test_placement(box_boundaries, atoms_positions):
        """Ensure that atoms are placed within the box"""
        for atom_position in atoms_positions:
            for x, boundary in zip(atom_position, box_boundaries):
                assert (x >= boundary[0]) & (x <= boundary[1]), f"Test failed: Atoms outside of the box"
        print("Test passed")
    test_placement(init.box_boundaries, init.atoms_positions)

.. label:: end_test_3a_class
