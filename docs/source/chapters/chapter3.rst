Create the box
==============

.. container:: justify

    The simulations box is defined and populated with atoms. 

Define the box dimensions and atom number
-----------------------------------------

.. container:: justify

    First, let us modify the *__init__* method of the *InitializeSimulation* class
    to allow for 3 new parameters to be defined: :math:`L_x, L_y`, and :math:`L_z`,
    i.e. the box dimensions along *x*, *y*, and *z*, respectively.
    We also define the number of atoms, *number_atoms*, which should be
    provided by users.

.. code-block:: python

    def __init__(self,
                 number_atoms
                 Lx,
                 epsilon=[0.1],
                 sigma=[1],
                 atom_mass=[1],
                 Ly=None,
                 Lz=None,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs) 

        self.number_atoms = number_atoms
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = 3
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units()
        self.define_box()
        self.identify_atom_properties()
        self.populate_box()

.. container:: justify

    A variable named *dimensions* was also created, and its value is enforced to be 3.
    In the future, the code could be adapted to create 2D or 1D simulations. 

    The 3 new methods called *define_box*, *identify_atom_properties*, and *populate_box*
    will be defined below.

    The *nondimensionalize_units* units method must be modified to allow for
    the normalization of *Lx*, *Ly*, and *Lz*:

.. code-block:: python

    def nondimensionalize_units(self):
        self.Lx /= self.reference_distance
        if self.Ly is not None:
            self.Ly /= self.reference_distance
        if self.Lz is not None:
            self.Lz /= self.reference_distance
        epsilon, sigma, atom_mass = [], [], []
        for e0, s0, m0 in zip(self.epsilon, self.sigma, self.atom_mass):
            epsilon.append(e0/self.reference_energy)
            sigma.append(s0/self.reference_distance)
            atom_mass.append(m0/self.reference_mass)
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

.. container:: justify

    If *Ly* or *Lz* has a value of *None*, it is not normalized.

Identify the atom properties
----------------------------

.. container:: justify

    Since it is anticipated that the code allows for more than one types
    of atoms, let us store the information for each atoms in arrays.
    Add the following method to the *InitializeSimulation* class.

.. code-block:: python

    def identify_atom_properties(self):
        self.total_number_atoms = np.sum(self.number_atoms)
        atoms_sigma = []
        atoms_epsilon = []
        atoms_mass = []
        atoms_type = []
        for sigma, epsilon, mass, number_atoms, type in zip(self.sigma, self.epsilon,
                                                            self.atom_mass, self.number_atoms,
                                                            np.arange(len(self.number_atoms))+1):
            atoms_sigma += [sigma] * number_atoms
            atoms_epsilon += [epsilon] * number_atoms
            atoms_mass += [mass] * number_atoms
            atoms_type += [type] * number_atoms
        self.atoms_sigma = np.array(atoms_sigma)
        self.atoms_epsilon = np.array(atoms_epsilon)
        self.atoms_mass = np.array(atoms_mass)
        self.atoms_type = np.array(atoms_type)

.. container:: justify

    Here, *atoms_sigma*, *atoms_epsilon*, *atoms_mass*, and *atoms_type* 
    are arrays of lengths equal to the total number of atoms. They will be used,
    among others, for calculating the forces between the atoms.

Define the box
--------------

.. container:: justify

    The values for *Lx*, *Ly*, and *Lz* that are provided by the
    users are used to define the box boundaries and box size.
    Add the following method to the *InitializeSimulation* class.

.. code-block:: python

     def define_box(self):
        box_boundaries = np.zeros((self.dimensions, 2))
        for dim, L in zip(range(self.dimensions),
                          [self.Lx, self.Ly, self.Lz]):
            if L is not None:
                box_boundaries[dim] = -L/2, L/2
            else:
                box_boundaries[dim] = -self.Lx/2, self.Lx/2
        self.box_boundaries = box_boundaries
        self.box_size = np.diff(box_boundaries).reshape(3)

.. container:: justify

    If *Ly* or *Lz* are not provided, their values is *None*, and
    *Lx* is used along the *y* or *z* direction.

Populate the box
----------------

.. container:: justify

    The previously defined box is populated with a number of atoms.

.. code-block:: python

    def populate_box(self):
        atoms_positions = np.zeros((self.total_number_atoms, self.dimensions))
        for dim in np.arange(self.dimensions):
            atoms_positions[:, dim] = np.random.random(self.total_number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2    
        self.atoms_positions = atoms_positions

Test the code
-------------

.. container:: justify

    Let us call the *MolecularDynamics* class and 
    create a cubic box with 20 atoms, and print the box size and
    the position of the first atom:

.. code-block:: python

    from MolecularDynamics import MolecularDynamics

    md = MolecularDynamics(number_atoms=[20],
                        Lx=30,
                        sigma=[3],
                        epsilon=[0.1],
                        atom_mass=[1],
                        data_folder = "md-output/")
    md.run()
    print("normalized box size:", md.box_size)
    print("first atom position:", md.atoms_positions[0])

.. container:: justify

    This is what I see:

.. code-block:: python

    normalized box size: [10. 10. 10.]
    first atom position: [ 0.36546066 -1.65788067  3.49811286]