Prepare the simulation
======================

.. container:: justify

    For simplicity, all the parameters that are specified as inputs by the user
    will be normalized. This makes all the calculations much simpler. The
    normalization is done within the *Prepare* class, which is defined here.

Unit systems
------------

.. container:: justify

    Two separate unit systems are used. The first unit system is called
    *real*. It follows the convention from the |lammps-unit-systems|.
    All input parameters are to be provided to the code in *real*
    units, for which:

    - masses are in grams/mole,
    - distances are in Ångstrom,
    - the time is in femtoseconds,
    - energies are in kcal/mol,
    - velocities are in Ångstrom/femtosecond,
    - forces are in (kcal/mol)/Ångstrom,
    - the temperature is in Kelvin,
    - the pressure is in atmospheres,
    - and the density is in g/cm^dim.

.. |lammps-unit-systems| raw:: html

   <a href="https://docs.lammps.org/units.html" target="_blank">LAMMPS unit systems</a>

.. container:: justify

    The *real* unit system is conventional in molecular simulations. However,
    it is not practical to perform calculations with such a complex unit system,
    as it would involve complicated prefactors. Instead, the LJ (for Lennard-Jones)
    unit system will be used. With the LJ unit systems, all quantities are
    unitless. All masses, distances, and energies are specified as multiples 
    of :math:`m`, :math:`\sigma`, and :math:`\epsilon`, which are the mass and LJ
    parameters of the atoms. Other quantities are specified from these 3 parameters,
    such as:

    - the time is in :math:`\sqrt{m \sigma^2 / \epsilon}`,
    - energies are in :math:`\epsilon`,
    - velocities are in :math:`\sqrt{m \sigma^4 / \epsilon}`,
    - forces are in :math:`\epsilon/\sigma`,
    - the temperature is in :math:`\epsilon/k_\text{B}`,
    - the pressure is in :math:`\epsilon/\sigma^3`,
    - and the density is in :math:`m/\sigma^3`.

    where :math:`k_\text{B}` is the Boltzmann constant. 

Start coding
------------

.. container:: justify

    Let us create a class named *Prepare*. To make the
    unit conversion easier, let us also import *numpy*, as
    well as the *constants* library from *numpy*.

    In a Python file named *Prepare.py*, copy the following lines:

.. code-block:: python

    import numpy as np
    from scipy import constants as cst

    import warnings
    warnings.filterwarnings('ignore')


    class Prepare:

.. container:: justify

    Four parameters are given to the Prepare class,
    the atom masses :math:`m`, the LJ parameters
    :math:`\sigma` and :math:`\epsilon`, and the
    number of atoms. These quantities must be provided as 
    lists, which will be useful later when we want to mix
    atoms of different types within the same simulation box.

.. container:: justify

    Within the *Prepare* class, add the following *__init__*
    method:  

.. code-block:: python

    def __init__(self,
                number_atoms=[10],  # List
                epsilon=[0.1],  # List - Kcal/mol
                sigma=[1],  # List - Angstrom
                atom_mass=[1],  # List - g/mol
                *args,
                **kwargs):
        self.number_atoms = number_atoms
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass
        super().__init__(*args, **kwargs)

.. container:: justify

    All four lists, *number_atoms*, *epsilon*, *sigma*, and *atom_mass* are
    given default values of :math:`10`,
    :math:`0.1~\text{[Kcal/mol]}`,
    :math:`1~\text{[Å]}`,
    and :math:`0.1~\text{[g/mol]}`, respectively. All four parameters are passed
    as *self*, which will allow for other methods to access them. Here, *args* and
    *kwargs* are used to accept an arbitrary number of positional
    and keyword arguments, respectively.

Calculate LJ units prefactors
-----------------------------

.. container:: justify

    Let us create a method called *calculate_LJunits_prefactors* that will be
    used to calculate the prefactors necessary to convert units from the *real*
    unit system to the *LJ* unit system.

.. container:: justify

    Within the *Prepare* class, copy the following method:

.. code-block:: python

    def calculate_LJunits_prefactors(self):
        self.reference_distance = self.sigma[0]  # Angstrom
        self.reference_energy = self.epsilon[0]  # Kcal/mol
        self.reference_mass = self.atom_mass[0]  # g/mol
        mass_kg = self.atom_mass[0]/cst.kilo/cst.Avogadro  # kg
        epsilon_J = self.epsilon[0]*cst.calorie*cst.kilo/cst.Avogadro  # J
        sigma_m = self.sigma[0]*cst.angstrom  # m
        time_s = np.sqrt(mass_kg*sigma_m**2/epsilon_J)  # s
        self.reference_time = time_s / cst.femto  # fs
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo  # kCal/mol/K
        self.reference_temperature = self.epsilon[0]/kB  # K
        pressure_pa = epsilon_J/sigma_m**3  # Pa
        self.reference_pressure = pressure_pa/cst.atm  # atm

.. container:: justify

    This method defines the *reference_distance* as the first element in the
    *sigma* list, i.e. :math:`\sigma_{11}`. Therefore atoms of type one will
    always be used for the normalization. Similarly, the first element
    in the *epsilon* list (:math:`\epsilon_{11}`) is used as a *reference_energy*, 
    and the first element in the *atom_mass* list (:math:`m_1`) is used as *reference_mass*.
    Then, the *reference_time* in femtosecond is calculated as :math:`\sqrt{m_1 \sigma_{11}^2 / \epsilon_{11}}`,
    and the *reference_pressure* is atmospheres is calculated as :math:`\epsilon_{11}/\sigma_{11}^3`.

.. container:: justify

    Finally, let us call the *calculate_LJunits_prefactors()*
    by adding the following line to the *__init__* method:

.. code-block:: python

    self.calculate_LJunits_prefactors()

.. container:: justify

    Every time the *Prepare* class will be initialized, all five reference values
    will be calculated and passed as *self*. 

Nondimensionalize units
-----------------------

.. container:: justify

    Let us take advantage of the calculated reference values and normalize the 
    three inputs of the *Prepare* class that have a physical dimension, i.e.
    *epsilon*, *sigma*, and *atom_mass*.

.. container:: justify

    Create a new method called 

.. code-block:: python

   def nondimensionalize_units_0(self):
        r"""Use LJ prefactors to convert units into non-dimensional."""
        # Normalize LJ properties
        epsilon, sigma, atom_mass = [], [], []
        for e0, s0, m0 in zip(self.epsilon, self.sigma, self.atom_mass):
            epsilon.append(e0/self.reference_energy)
            sigma.append(s0/self.reference_distance)
            atom_mass.append(m0/self.reference_mass)
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

.. code-block:: python


    self.nondimensionalize_units_0()
    self.calculate_cross_coefficients()

.. code-block:: python

    def nondimensionalize_units_0(self):
        r"""Use LJ prefactors to convert units into non-dimensional."""
        # Normalize LJ properties
        epsilon, sigma, atom_mass = [], [], []
        for e0, s0, m0 in zip(self.epsilon, self.sigma, self.atom_mass):
            epsilon.append(e0/self.reference_energy)
            sigma.append(s0/self.reference_distance)
            atom_mass.append(m0/self.reference_mass)
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

    def identify_atom_properties(self):
        r"""Create initial atom array from input parameters"""
        self.total_number_atoms = np.sum(self.number_atoms)
        atoms_sigma = []
        atoms_epsilon = []
        atoms_mass = []
        atoms_type = []
        for parts in zip(self.sigma,
                        self.epsilon,
                        self.atom_mass,
                        self.number_atoms,
                        np.arange(len(self.number_atoms))+1):
            sigma, epsilon, mass, number_atoms, type = parts
            atoms_sigma += [sigma] * number_atoms
            atoms_epsilon += [epsilon] * number_atoms
            atoms_mass += [mass] * number_atoms
            atoms_type += [type] * number_atoms
        self.atoms_sigma = np.array(atoms_sigma)
        self.atoms_epsilon = np.array(atoms_epsilon)
        self.atoms_mass = np.array(atoms_mass)
        self.atoms_type = np.array(atoms_type)

    def calculate_cross_coefficients(self):
        r"""The LJ cross coefficients are calculated and returned as arrays"""
        self.identify_atom_properties()
        epsilon_ij = []
        for i in range(self.total_number_atoms):
            epsilon_i = self.atoms_epsilon[i]
            for j in range(i + 1, self.total_number_atoms):
                epsilon_j = self.atoms_epsilon[j]
                epsilon_ij.append((epsilon_i+epsilon_j)/2)
        self.array_epsilon_ij = np.array(epsilon_ij)
        sigma_ij = []
        for i in range(self.total_number_atoms):
            sigma_i = self.atoms_sigma[i]
            for j in range(i + 1, self.total_number_atoms):
                sigma_j = self.atoms_sigma[j]
                sigma_ij.append((sigma_i+sigma_j)/2)
        self.array_sigma_ij = np.array(sigma_ij)


