LJ parameters and units
=======================

.. container:: justify

    The normalization of the basic parameters is defined here. 

Unit systems
------------

.. container:: justify

    Two separate unit systems are used. The first unit system is called
    *real*, and followed the convention from LAMMPS. All inputs will be 
    provided to the code in *real* units, for which:

    - mass in grams/mole
    - distance in Ångstroms
    - time in femtoseconds
    - energy in kcal/mol
    - velocity in Ångstroms/femtosecond
    - force in (kcal/mol)/Angstrom
    - temperature in Kelvin
    - pressure in atmospheres
    - density in g/cm^dim

.. container:: justify

    The *real* unit system is conventional in molecular simulations. However,
    it is not practical to perform the calculations with such complex units,
    as it would involve complex prefactors. Instead, the LJ (for Lennard-Jones)
    unit system will be used. With the LJ unit systems, all quantities are
    unitless. All masses, distances, and energies are specified as multiples 
    of :math:`m`, :math:`\sigma`, and :math:`\epsilon`, which are the mass and LJ
    parameters of the atoms. Other quantities are specified from these 3 quantities:

    - mass in :math:`m`
    - distance in :math:`\sigma`
    - time in :math:`\sqrt{m \sigma^2 / \epsilon}`
    - energy in :math:`\epsilon`
    - velocity in :math:`\sqrt{m \sigma^4 / \epsilon}`
    - force in :math:`\epsilon/\sigma`
    - temperature in :math:`\epsilon/k_\text{B}`
    - pressure in :math:`\epsilon/\sigma^3`
    - density in :math:`m/\sigma^3`

    where :math:`k_\text{B}` is the Boltzmann constant. 

Lennard-Jones parameters
------------------------

.. container:: justify

    Let us modify the *InitializeSimulation* class
    and allow for the :math:`m`, :math:`\sigma`, and :math:`\epsilon`
    parameters to be specified:

.. code-block:: python

    from scipy import constants as cst
    import numpy as np

    class InitializeSimulation:
        def __init__(self,
                    epsilon=[0.1],
                    sigma=[1],
                    atom_mass=[1],
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs) 

            self.epsilon = epsilon
            self.sigma = sigma
            self.atom_mass = atom_mass

            self.calculate_LJunits_prefactors()
            self.nondimensionalize_units()

.. container:: justify

    First, the NumPy library is imported, together with the constant module of SciPy.

    Then, :math:`m`, :math:`\sigma`, and :math:`\epsilon` are provided
    to *InitializeSimulation* as lists, respectively called *atom_mass*, *sigma*, and *epsilon*.
    Each list is made of a single element with default value :math:`1~\text{g/mol}`,
    :math:`1~\text{Å}`, and :math:`1~\text{kcal/mol}`.

    The two methods *calculate_LJunits_prefactors* and *nondimensionalize_units()*,
    which will be defined below, are called and will be used
    to convert the *real* units into *LJ* units.

Calculate the LJ prefactors
---------------------------

.. container:: justify

    Here, prefactors needed for the normalization are calculated
    from :math:`m`, :math:`\sigma`, and :math:`\epsilon`. Add the
    following method to the *InitializeSimulation* class:

.. code-block:: python

    def calculate_LJunits_prefactors(self):
        self.reference_distance = self.sigma[0] # Angstrom
        self.reference_energy = self.epsilon[0] # Kcal/mol
        self.reference_mass = self.atom_mass[0] # g/mol

.. container:: justify

    Here, *reference_distance*, *reference_energy*, and *reference_mass*
    are calculated and saved into *self*.
    
    If more than one atom type is used, and the lists *sigma*, *epsilon*, and *atom_mass*
    contain more than one element, the first values will be used for the normalization.

Non-dimensionalize the units
----------------------------

.. container:: justify

    Here, the initial values for :math:`m`, :math:`\sigma`, and :math:`\epsilon`
    are normalized. Add the following method to the *InitializeSimulation* class.

.. code-block:: python

    def nondimensionalize_units(self):
        epsilon, sigma, atom_mass = [], [], []
        for e0, s0, m0 in zip(self.epsilon, self.sigma, self.atom_mass):
            epsilon.append(e0/self.reference_energy)
            sigma.append(s0/self.reference_distance)
            atom_mass.append(m0/self.reference_mass)
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

Test the code
-------------

.. container:: justify

    Let us make sure that the normalization is properly made
    by our script, by calling the *MolecularDynamics* class:

.. code-block:: python

    from MolecularDynamics import MolecularDynamics

    md = MolecularDynamics(sigma=[3],
                        epsilon=[0.1],
                        atom_mass=[1],
                        data_folder = "md-output/")
    md.run()
    print("normalized epsilon value:", md.epsilon[0])
    print("normalized sigma value:", md.sigma[0])
    print("normalized mass value:", md.atom_mass[0])

.. container:: justify

    If it works, the code should return:

.. code-block:: python

    normalized epsilon value: 1.0
    normalized sigma value: 1.0
    normalized mass value: 1.0