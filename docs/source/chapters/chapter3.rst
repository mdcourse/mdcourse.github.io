Create the box
==============

.. container:: justify

    The simulations box is defined and populated with atoms. 

Define the box dimensions
-------------------------

.. container:: justify

    First, let us modify the *__init__* method of the *InitializeSimulation* class
    to allow for 3 new parameters to be defined: :math:`L_x, L_y`, and :math:`L_z`.

.. code-block:: python

    def __init__(self,
                 Lx,
                 epsilon=[0.1],
                 sigma=[1],
                 atom_mass=[1],
                 Ly=None,
                 Lz=None,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs) 

        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = 3
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units()

.. container:: justify

    A variable named *dimensions* was also created, and its value is enforced to be 3.
    In the future, the code could be adapted to create 2D or 1D simulations.  

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

Define the box
--------------

.. container:: justify

    Here, the values for *Lx*, *Ly*, and *Lz* that are provided by the
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