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
    the box dimensions which is a list of size three with parameters in units of Angstrom, 
    a seed, and some initial positions for the atoms that can be provided as an array
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

Final code
----------

.. container:: justify

    After following these steps, this is what the final code should
    look like. For clarity, some comments and descriptions were added for each
    method.

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
