Initialize simulation
=====================

The *InitializeSimulation* class deals with creating the simulation box,
populating it with a chosen number of atoms with initial velocity. The
*InitializeSimulation* class also takes care of converting the inputs units 
into non-dimensionalized.

The *__init__* function
-----------------------

Open a new blank *.py* file, call it *main.py*, and copy the following
files in it: 

.. code-block:: python

    import numpy as np

    import warnings
    warnings.filterwarnings('ignore')

The *numpy* functions will be used here, and some annoying warning
are deactivated (will be useful when evaluating Metropolis criteria).
Then, copy:

.. code-block:: python

    class InitializeSimulation:
        def __init__(self,
                     number_atoms,
                     Lx,
                     dimensions=3,
                     Ly = None,
                     Lz = None,
                     epsilon=0.1,
                     sigma=1,
                     atom_mass=1,
                     seed=None,
                     desired_temperature=300,
                     desired_pressure=1)

        self.number_atoms = number_atoms
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = dimensions
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass
        self.seed = seed
        self.desired_temperature = desired_temperature
        self.desired_pressure = desired_pressure

Important variables are being defined within the *__init__* function of the InitializeSimulation class. 
*number_atoms* is the desired initial number of atoms, *Lx* the box size along the *x* direction (or 
along all 3 directions if *Ly* and *Lz* are not provided, which is the default), *dimensions* is 
the number of dimensions of the system, *epsilon* the reference Lennard-Jones energy parameter in 
Kcal/mol, *sigma* the reference Lennard-Jones distance parameter in Angstrom, and *atom_mass* the 
reference mass in g/mol. The *desired_temperature* and *desired_pressure* parameters, respectively in
Kelvin and in atmosphere, are optional parameters that will be used in some case (as will be seen later
on), and *seed* offer the possibility to fix the seed of the randomly generated number, allowing for
repeating the same simulation if necessary -- this is usually useful during debugging. 

Add the following lines to the *__init__* function:

.. code-block:: python

        if self.seed is not None:
            np.random.seed(self.seed)

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

Populate box
------------

The *populate_box* function add a number corresponding to *number_atoms* of atoms in the 
box at random positions. 

.. code-block:: python

    def populate_box(self):
        """Place atoms at random positions within the box."""
        atoms_positions = np.zeros((self.number_atoms, self.dimensions))
        if self.provided_positions is not None:
            atoms_positions = self.provided_positions/self.reference_distance
        else:
            for dim in np.arange(self.dimensions):
                atoms_positions[:, dim] = np.random.random(self.number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2    
        self.atoms_positions = atoms_positions

More advanced version of the *populate_box* function could for instance ensure that no overlap exists 
between the atoms. For simplicity, we stick here with simple positioning of the atoms. Another more
advanced alternative is to place the atoms on given lattice, like simple cubic lattice, which is 
useful for studying solids. Another common option is to restrain the placement of the atoms within a
certain sub-region of the system.

Give initial velocity to the atoms
----------------------------------

Providing the atoms with an initial velocity is useful to reach the target temperature faster.

.. code-block:: python

    def give_velocity(self):
        """Give velocity to atoms so that the initial temperature is the desired one."""
        atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
        if self.provided_velocities is not None:
            atoms_velocities = self.provided_velocities/self.reference_distance*self.reference_time
        else:
            for dim in np.arange(self.dimensions):  
                atoms_velocities[:, dim] = np.random.normal(size=self.number_atoms)
        self.atoms_velocities = atoms_velocities
        self.calculate_temperature()
        scale = np.sqrt(1+((self.desired_temperature/self.temperature)-1))
        self.atoms_velocities *= scale

Commonly, one can make sure that no overall translational nor rotational momentum is given to the
atoms. 

