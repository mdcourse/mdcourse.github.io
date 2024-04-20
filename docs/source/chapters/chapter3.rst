Initialize the simulation
=========================

.. container:: justify

    Here the *InitializeSimulation* class is created. This class is used to
    prepare the simulation box and populate the box with atoms.

Simulation box
--------------

.. container:: justify

    In all simulations performed here, a certain number of atoms are first 
    randomly placed in a cuboid box.

Start coding
------------

.. container:: justify

    Let us start by importing NumPy, as well as the previously created *Prepare*
    class:

.. code-block:: python

    import numpy as np
    from Prepare import Prepare

.. container:: justify

    Then, let us create the *InitializeSimulation* class that is inheriting
    the *Prepare* class.

.. code-block:: python

    class InitializeSimulation(Prepare):
        def __init__(self,
                    box_dimensions=[10, 10, 10],  # List - Angstroms
                    seed=None,  # Int
                    initial_positions=None,  # Array - Angstroms
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)
            self.box_dimensions = box_dimensions
            self.dimensions = len(box_dimensions)
            self.seed = seed
            self.initial_positions = initial_positions

.. container:: justify

    Three parameters are provided to the *InitializeSimulation* class, namely
    the box dimensions which is a list with a length corresponding to *dimensions* (default is 3)
    with parameters in units of Angstrom, a seed, and some initial positions for the atoms that can be provided as an array
    of three columns. If *initial_positions* is
    left equal to *None*, positions with be randomly attributed to the atoms.

Seed
----

.. container:: justify

    The *seed* can be provided to launch the same simulations several
    times in a row, which is useful during debugging. Add the following *if*
    condition to the *__init__()* method.

.. code-block:: python

    (...)
    self.initial_positions = initial_positions
    if self.seed is not None:
        np.random.seed(self.seed)

.. container:: justify

    If a seed is provided, it is passed to the *random.seed()* function of *NumPy*.
    If the seed is left to None, the simulation will be randomized.

Define the box
--------------

.. container:: justify

    Let us define a box from the dimensions provided as a list named *box_dimensions*.

.. code-block:: python

    def define_box(self):
        box_boundaries = np.zeros((self.dimensions, 2))
        for dim, L in zip(range(self.dimensions), self.box_dimensions):
            box_boundaries[dim] = -L/2, L/2
        self.box_boundaries = box_boundaries
        box_size = np.diff(self.box_boundaries).reshape(3)
        box_geometry = np.array([90, 90, 90])
        self.box_size = np.array(box_size.tolist()+box_geometry.tolist())

.. container:: justify

    By symmetry, the box is centered in 0 for all axes. A *box_size* is also
    defined. It follows the  MDAnalysis conventions; Lx, Ly, Lz, 90, 90, 90,
    where the last three numbers are angles in degrees that are usually used to
    define triclinic (non-orthogonal) boxes, which is not a possibility of the current code.

Populate the box
----------------

.. container:: justify

    Here, a number of atoms are placed within the simulation box. If initial
    positions were not provided (i.e. *initial_positions = None*), atoms
    are placed randomly within the box. If initial positions were provided
    as an array, they are used instead. Note that in that case, the array
    number be of size 'number of atoms' x ''number of dimensions.

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

Final code
----------

.. container:: justify

    After following these steps, this is what the final code should
    look like. For clarity, some comments and descriptions were added for each
    method.

.. label:: start_InitializeSimulation_class

.. code-block:: python

    import numpy as np
    from Prepare import Prepare


    class InitializeSimulation(Prepare):
        def __init__(self,
                    box_dimensions=[10, 10, 10],  # List - Angstroms
                    seed=None,  # Int
                    initial_positions=None,  # Array - Angstroms
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)
            self.box_dimensions = box_dimensions
            self.dimensions = len(box_dimensions)
            self.seed = seed
            self.initial_positions = initial_positions
            if self.seed is not None:
                np.random.seed(self.seed)
            self.nondimensionalize_units_1()
            self.define_box()
            self.populate_box()

        def nondimensionalize_units_1(self):
            """Use LJ prefactors to convert units into non-dimensional."""
            # Normalize box dimensions
            box_dimensions = []
            for L in self.box_dimensions:
                box_dimensions.append(L/self.reference_distance)
            self.box_dimensions = box_dimensions

        def define_box(self):
            """Define box boundaries based on the box dimensions."""
            box_boundaries = np.zeros((self.dimensions, 2))
            for dim, L in zip(range(self.dimensions), self.box_dimensions):
                box_boundaries[dim] = -L/2, L/2
            self.box_boundaries = box_boundaries
            # Also define the box size following MDAnalysis conventions
            box_size = np.diff(self.box_boundaries).reshape(3)
            box_geometry = np.array([90, 90, 90])
            self.box_size = np.array(box_size.tolist()+box_geometry.tolist())

        def populate_box(self):
            """Place atoms at random positions within the box."""
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

Test the code
-------------

.. container:: justify

    Let us test the *InitializeSimulation* class to make sure that it does what
    is expected.

.. label:: start_test_InitializeSimulation_class

.. code-block:: python

    from InitializeSimulation import InitializeSimulation

    self = InitializeSimulation(number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 1], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    print("Atom positions:")
    print(self.atoms_positions)

.. container:: justify

    Which should return:

.. code-block:: python

    Atom positions:
    [[-1.15270975  1.25033545  0.39460297]
    [ 2.10225087 -2.12285757 -2.43760443]
    [ 0.86169508 -0.77310475 -0.74742818]
    [ 0.81255861  2.26285536  1.76611306]
    [-0.31367217 -1.55867269 -2.71347742]]

.. label:: end_test_InitializeSimulation_class