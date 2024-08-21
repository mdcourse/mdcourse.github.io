Prepare the simulation
======================

For simplicity, all the parameters that are specified as inputs by the user
will be non-dimensionalized. This makes the calculations, such as force evaluation,
much simpler. Here, the non-dimensionalization is done within the *Prepare* class

Note: Although this is a necessary step that will make our life much easier later,
this is by far the less exciting part of the code.

Unit systems
------------

In this code, two unit systems are used: *real* and *LJ*, where *LJ* stands
for Lennard-Jones. The *real* unit system is used for both inputs and outputs,
and the *LJ* unit system is used to all the internal calculations.

The *real* unit system follows the convention from the |lammps-unit-systems|,
and:

- masses are in grams/mole,
- distances are in Ångstrom,
- the time is in femtoseconds,
- energies are in kcal/mol,
- velocities are in Ångstrom/femtosecond,
- forces are in (kcal/mol)/Ångstrom,
- the temperature is in Kelvin,
- the pressure is in atmospheres,
- the density is in g/cm^dim.

.. |lammps-unit-systems| raw:: html

   <a href="https://docs.lammps.org/units.html" target="_blank">LAMMPS unit systems</a>

The *real* unit system is conventional in molecular simulations. However,
it would not be practical to perform calculations with such a complex unit system,
as it would involve complicated prefactors. Instead, the LJ unit system will be
used for all calculations. With the LJ unit systems, all quantities are
unitless: all masses, distances, and energies are specified as multiples 
of :math:`m`, :math:`\sigma`, and :math:`\epsilon`, which are the mass and LJ
parameters of the atoms. Other quantities are specified from these 3 parameters:

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

Let us fill the previously created class named *Prepare*. To make the
unit conversion easier, let us import *numpy*, as
well as the *constants* library from *numpy*.

In the file named *Prepare.py*, copy the following lines:

.. label:: start_Prepare_class

.. code-block:: python

    import numpy as np
    from scipy import constants as cst

.. label:: end_Prepare_class

Here, the |NumPy| library is imported, together with the constants module of |SciPy|.

.. |NumPy| raw:: html

   <a href="https://numpy.org/" target="_blank">NumPy</a>

.. |SciPy| raw:: html

   <a href="https://scipy.org/" target="_blank">SciPy</a>

Four parameters are given to the *Prepare* class:

- the atom masses :math:`m`,
- the LJ parameters :math:`\sigma` and :math:`\epsilon`,
- and the number of atoms.

All these quantities must be provided as  lists, which will be useful later when we want to mix
atoms of different types within the same simulation box.

Modify the *Prepare* class as follows:  

.. label:: start_Prepare_class

.. code-block:: python

    class Prepare:
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

.. label:: end_Prepare_class

Here the four lists, *number_atoms*, *epsilon*, *sigma*, and *atom_mass* are
given default values of :math:`10`, :math:`0.1~\text{[Kcal/mol]}`, :math:`1~\text{[Å]}`,
and :math:`0.1~\text{[g/mol]}`, respectively.

All four parameters are passed as *self*, which will allow for other methods to
access them. Here, *args* and *kwargs* are used to accept an arbitrary number of
positional and keyword arguments, respectively.

Calculate LJ units prefactors
-----------------------------

Within the *Prepare* class, let us create a method called *calculate_LJunits_prefactors*
that will be used to calculate the prefactors necessary to convert units from the *real*
unit system to the *LJ* unit system.

Within the *Prepare* class, copy the following method:

.. label:: start_Prepare_class

.. code-block:: python

    def calculate_LJunits_prefactors(self):
        # Distance, energy, and mass
        self.reference_distance = self.sigma[0]  # Angstrom
        self.reference_energy = self.epsilon[0]  # Kcal/mol
        self.reference_mass = self.atom_mass[0]  # g/mol
        # Time
        mass_kg = self.atom_mass[0]/cst.kilo/cst.Avogadro  # kg
        epsilon_J = self.epsilon[0]*cst.calorie*cst.kilo/cst.Avogadro  # J
        sigma_m = self.sigma[0]*cst.angstrom  # m
        time_s = np.sqrt(mass_kg*sigma_m**2/epsilon_J)  # s
        self.reference_time = time_s / cst.femto  # fs
        # Temperature
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo  # kCal/mol/K
        self.reference_temperature = self.epsilon[0]/kB  # K
        # Pressure
        pressure_pa = epsilon_J/sigma_m**3  # Pa
        self.reference_pressure = pressure_pa/cst.atm  # atm

.. label:: end_Prepare_class

This method defines the *reference_distance* as the first element in the
*sigma* list, i.e. :math:`\sigma_{11}`. Therefore atoms of type one will
always be used for the normalization. Similarly, the first element
in the *epsilon* list (:math:`\epsilon_{11}`) is used as a *reference_energy*, 
and the first element in the *atom_mass* list (:math:`m_1`) is used as *reference_mass*.
Then, the *reference_time* in femtosecond is calculated as :math:`\sqrt{m_1 \sigma_{11}^2 / \epsilon_{11}}`,
the reference temperature in Kelvin as :math:`\epsilon_{11} / k_\text{B}`,
and the *reference_pressure* in atmospheres is calculated as :math:`\epsilon_{11}/\sigma_{11}^3`.

Finally, let us make sure that the *calculate_LJunits_prefactors()* method is
called systematically by adding the following line to the *__init__()* method:

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        super().__init__(*args, **kwargs)
        self.calculate_LJunits_prefactors()

.. label:: end_Prepare_class

Every time the *Prepare* class will be initialized, all reference values
will be calculated and passed as *self*. 

Nondimensionalize units
-----------------------

Let us take advantage of the calculated reference values and normalize the 
three inputs of the *Prepare* class that have a physical dimension, i.e.
*epsilon*, *sigma*, and *atom_mass*.

Create a new method called *nondimensionalize_units_0* within the *Prepare*
class. The index *0* is used to differentiate this method from the other methods
that will be used to nondimensionalize units in future classes. 

.. label:: start_Prepare_class

.. code-block:: python

    def nondimensionalize_units_0(self):
        # Normalize LJ properties
        epsilon, sigma, atom_mass = [], [], []
        for e0, s0, m0 in zip(self.epsilon, self.sigma, self.atom_mass):
            epsilon.append(e0/self.reference_energy)
            sigma.append(s0/self.reference_distance)
            atom_mass.append(m0/self.reference_mass)
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass

.. label:: end_Prepare_class

Here, we anticipate that *epsilon*, *sigma*, and *atom_mass* may contain
more than one element in the future, and normalize each element with the
corresponding reference value. The *zip()* function allows us to loop over
all three lists at once.  

Let us also call the *nondimensionalize_units_0* from the *__init__()* method:

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units_0()

.. label:: end_Prepare_class

Identify atom properties
------------------------

Anticipating the future use of multiple atom types, where each type will be
associated with its own :math:`\sigma`, :math:`\epsilon` and  :math:`m`,
let us create arrays containing the properties of each atom in the simulation. 
For instance, in the case of a simulation with two atoms of type 1 and three
atoms of type 2, the corresponding *atoms_sigma* will be:

.. math::

    \text{atoms_sigma} = [\sigma_{11}, \sigma_{11}, \sigma_{22}, \sigma_{22}, \sigma_{22}]

where :math:`\sigma_{11}` and :math:`\sigma_{22}` are the sigma values for 
atoms of type 1 and 2 respectively. The *atoms_sigma* array will allow
for future calculation of force.

Create a new method called *identify_atom_properties*, and place it
within the *Prepare* class:

.. label:: start_Prepare_class

.. code-block:: python

    def identify_atom_properties(self):
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
    
.. label:: end_Prepare_class
    
Let us call the *identify_atom_properties* from the *__init__()* method:

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units_0()
        self.identify_atom_properties()

.. label:: end_Prepare_class

Calculate cross coefficients
----------------------------

Let us calculate all cross coefficients that can be used to calculate the interactions
between atom :math:`i` and atom :math:`j`. From the example described previously,
where:

.. math::

    \text{atoms_sigma} = [\sigma_{11}, \sigma_{11}, \sigma_{22}, \sigma_{22}, \sigma_{22}]

one expects all direct and cross coefficients to be:

.. math::

    \text{array_sigma_ij} = \begin{bmatrix}
            \sigma_{11} \text{(0-1)} & \sigma_{12} \text{(0-2)} & \sigma_{12} \text{(0-3)} & \sigma_{12} \text{(0-4)} \\
                                     & \sigma_{12} \text{(1-2)} & \sigma_{12} \text{(1-3)} & \sigma_{12} \text{(1-4)} \\
                                     &                          & \sigma_{22} \text{(2-3)} & \sigma_{22} \text{(2-4)} \\
                                     &                          &                          & \sigma_{22} \text{(3-4)} \\
        \end{bmatrix}


where it is assumed that the matrix is symmetric, so the coefficients in the bottom
left are not specified. The first value in the top left corner of the matrix,
:math:`\sigma_{11} \text{(0-1)}`, indicates that the :math:`\sigma` value for the interaction
between the first (0) and second atom (1) is :math:`\sigma_{11}`, as both atoms 0
and 1 are of type 1. A similar matrix can be written for epsilon_sigma_ij.

The values of the cross coefficients :math:`\sigma_{12}` and :math:`\epsilon_{12}`
are assumed to follow the arithmetic mean :

.. math::

    \sigma_{12} = (\sigma_{11}+\sigma_{22})/2 \\
    \epsilon_{12} = (\epsilon_{11}+\epsilon_{22})/2

Create the following method called *calculate_cross_coefficients* within the 
*Prepare* class:

.. label:: start_Prepare_class

.. code-block:: python

    def calculate_cross_coefficients(self):
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

.. label:: end_Prepare_class

After calling for the *identify_atom_properties()* method, double loops
are performed over all direct coefficients, and the cross coefficients
are stored within *array_sigma_ij* and *array_epsilon_ij*.

Finally, let us call the *calculate_cross_coefficients* method from the
*__init__()* method.

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units_0()
        self.identify_atom_properties()
        self.calculate_cross_coefficients()

.. label:: end_Prepare_class

Test the code
-------------

Let us test the *Prepare* class to make sure that it does what is expected.
This initiates a system with 2 atoms of type 1, and 3 atoms of type 2:

.. label:: start_test_Prepare_class

.. code-block:: python

    import numpy as np
    from Prepare import Prepare

    prep = Prepare(number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 9], # g/mol
        )
    assert prep.reference_energy == 0.1
    assert prep.reference_distance == 3
    assert prep.reference_mass == 1
    assert prep.total_number_atoms == 5
    assert np.int32(prep.total_number_atoms*(prep.total_number_atoms-1)/2) == len(self.array_epsilon_ij)

.. label:: end_test_Prepare_class

This test assert that the reference energy, distance, and mass were properly defined,
that the number of atoms indeed corresponds to the sum of the number of atoms of each
type, and that the length of the array_epsilon_ij array is consistent.
