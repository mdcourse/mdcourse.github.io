.. _chapter2-label:

Setting Up the Simulation
==========================

To streamline the simulation process, all user-specified parameters will be
non-dimensionalized. By removing units from the calculations, we simplify
complex operations like force evaluation, making the code more efficient and
easier to manage. This non-dimensionalization is handled within the *Prepare*
class.

While this step may not be the most thrilling aspect of the simulation, it is
essential groundwork that will significantly ease our work as we progress.

Unit systems
------------

In this code, two unit systems are employed: *real* and *LJ*, with *LJ* standing
for Lennard-Jones. The *real* unit system is used for both inputs and outputs,
while the *LJ* unit system is used for all internal calculations.

The *real* unit system follows the conventions outlined in the |lammps-unit-systems|:

- Masses are in grams per mole,
- Distances are in Ångströms,
- Time is in femtoseconds,
- Energies are in kcal/mol,
- Velocities are in Ångströms per femtosecond,
- Forces are in (kcal/mol)/Ångström,
- Temperature is in Kelvin,
- Pressure is in atmospheres,
- Density is in g/cm\ :sup:`3` (in 3D).

.. |lammps-unit-systems| raw:: html

   <a href="https://docs.lammps.org/units.html" target="_blank">LAMMPS unit systems</a>

The *real* unit system is conventional in molecular simulations. However, using
such a complex unit system for calculations would involve cumbersome prefactors.
To simplify, the *LJ* unit system is used for all calculations. In the *LJ* unit
system, all quantities are dimensionless. Masses, distances, and energies are
expressed as multiples of :math:`m`, :math:`\sigma`, and :math:`\epsilon`,
which represent the mass and LJ parameters of the atoms. Other quantities are
derived from these three parameters:

- Time is in :math:`tau = \sqrt{m \sigma^2 / \epsilon}`,
- Energies are in :math:`\epsilon`,
- Velocities are in :math:`\sigma / \tau`,
- Forces are in :math:`\epsilon / \sigma`,
- Temperature is in :math:`\epsilon / k_\text{B}`,
- Pressure is in :math:`\epsilon / \sigma^3`,
- Density is in :math:`m / \sigma^3`,

where :math:`k_\text{B}` is the Boltzmann constant.

Start coding
------------

Let's fill in the previously created class named Prepare. To facilitate unit
conversion, we will import |NumPy| and the constants module from |SciPy|.

.. |NumPy| raw:: html

   <a href="https://numpy.org/" target="_blank">NumPy</a>

.. |SciPy| raw:: html

   <a href="https://scipy.org/" target="_blank">SciPy</a>

In the file named *Prepare.py*, add the following lines:

.. label:: start_Prepare_class

.. code-block:: python

    import numpy as np
    from scipy import constants as cst

.. label:: end_Prepare_class

Four parameters are provided to the *Prepare* class:

- the atom masses :math:`m`,
- the LJ parameters :math:`\sigma` and :math:`\epsilon`,
- and the number of atoms.

All these quantities must be supplied as lists. This will be useful later when
we want to mix atoms of different types within the same simulation box.

Modify the *Prepare* class as follows:  

.. label:: start_Prepare_class

.. code-block:: python

    class Prepare:
        def __init__(self,
                    number_atoms=[10],  # List - no unit
                    epsilon=[0.1],  # List - Kcal/mol
                    sigma=[3],  # List - Angstrom
                    atom_mass=[10],  # List - g/mol
                    *args,
                    **kwargs):
            self.number_atoms = number_atoms
            self.epsilon = epsilon
            self.sigma = sigma
            self.atom_mass = atom_mass
            super().__init__(*args, **kwargs)

.. label:: end_Prepare_class

Here, the four lists *number_atoms* :math:`N`, *epsilon* :math:`\epsilon`,
*sigma* :math:`\sigma`, and *atom_mass* :math:`m` are given default values of
:math:`10`, :math:`0.1~\text{[Kcal/mol]}`, :math:`3~\text{[Å]}`, and
:math:`10~\text{[g/mol]}`, respectively.

All four parameters are assigned to *self*, allowing other methods to access
them. The *args* and *kwargs* parameters are used to accept an arbitrary number
of positional and keyword arguments, respectively.

Calculate LJ units prefactors
-----------------------------

Within the *Prepare* class, let us create a method called *calculate_LJunits_prefactors*
that will be used to calculate the prefactors necessary to convert units from the *real*
unit system to the *LJ* unit system:

.. label:: start_Prepare_class

.. code-block:: python

    def calculate_LJunits_prefactors(self):
        """Calculate the Lennard-Jones units prefacors."""
        # Define the reference distance, energy, and mass
        self.reference_distance = self.sigma[0]  # Angstrom
        self.reference_energy = self.epsilon[0]  # kcal/mol
        self.reference_mass = self.atom_mass[0]  # g/mol

        # Calculate the prefactor for the time
        mass_kg = self.atom_mass[0]/cst.kilo/cst.Avogadro  # kg
        epsilon_J = self.epsilon[0]*cst.calorie*cst.kilo/cst.Avogadro  # J
        sigma_m = self.sigma[0]*cst.angstrom  # m
        time_s = np.sqrt(mass_kg*sigma_m**2/epsilon_J)  # s
        self.reference_time = time_s / cst.femto  # fs

        # Calculate the prefactor for the temperature
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo  # kcal/mol/K
        self.reference_temperature = self.epsilon[0]/kB  # K

        # Calculate the prefactor for the pressure
        pressure_pa = epsilon_J/sigma_m**3  # Pa
        self.reference_pressure = pressure_pa/cst.atm  # atm

.. label:: end_Prepare_class

This method defines the reference distance as the first element in the
*sigma* list, i.e., :math:`\sigma_{11}`. Therefore, atoms of type one will
always be used for the normalization. Similarly, the first element
in the *epsilon* list (:math:`\epsilon_{11}`) is used as the reference energy,
and the first element in the *atom_mass* list (:math:`m_1`) is used as the
reference mass. Then, the reference_time in femtoseconds is calculated
as :math:`\sqrt{m_1 \sigma_{11}^2 / \epsilon_{11}}`, the reference temperature
in Kelvin as :math:`\epsilon_{11} / k_\text{B}`, and the reference_pressure
in atmospheres is calculated as :math:`\epsilon_{11}/\sigma_{11}^3`.

Finally, let us ensure that the *calculate_LJunits_prefactors* method is
called systematically by adding the following line to the *__init__()* method:

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        super().__init__(*args, **kwargs)
        self.calculate_LJunits_prefactors()

.. label:: end_Prepare_class

Every time the *Prepare* class is initialized, all reference values will
be calculated and stored as attributes of *self*.

Nondimensionalize units
-----------------------

Let us take advantage of the calculated reference values and normalize the
three inputs of the *Prepare* class that have physical dimensions, i.e.,
*epsilon*, *sigma*, and *atom_mass*.

Create a new method called *nondimensionalize_units_0* within the *Prepare*
class:

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

The index *0* is used to differentiate this method from other methods that
will be used to nondimensionalize units in future classes. We anticipate that
*epsilon*, *sigma*, and *atom_mass* may contain more than one element, so
each element is normalized with the corresponding reference value. The
*zip()* function allows us to loop over all three lists at once.

Let us also call the *nondimensionalize_units_0* from the *__init__()* method
of the *Prepare* class:

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
associated with its own :math:`\sigma`, :math:`\epsilon`, and :math:`m`, let
us create arrays containing the properties of each atom in the simulation. For
instance, in the case of a simulation with two atoms of type 1 and three atoms
of type 2, the corresponding *atoms_sigma* array will be:

.. math::

    \text{atoms_sigma} = [\sigma_{11}, \sigma_{11}, \sigma_{22}, \sigma_{22}, \sigma_{22}]

where :math:`\sigma_{11}` and :math:`\sigma_{22}` are the sigma values for
atoms of type 1 and 2, respectively. The *atoms_sigma* array will allow for
future calculations of force.

Create a new method called *identify_atom_properties*, and place it
within the *Prepare* class:

.. label:: start_Prepare_class

.. code-block:: python

    def identify_atom_properties(self):
        """Identify the atom properties for each atom."""
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
        self.nondimensionalize_units_0()
        self.identify_atom_properties()

.. label:: end_Prepare_class

..
    Calculate cross coefficients
    ----------------------------

    Let us calculate all cross coefficients that are required to calculate the interactions
    between atom :math:`i` and atom :math:`j`. From the example described previously,
    where:

    .. math::

        \text{atoms_sigma} = [\sigma_{11}, \sigma_{11}, \sigma_{22}, \sigma_{22}, \sigma_{22}]

    one expects all direct and cross coefficients to be:

    .. math::

        \text{array_sigma_ij} = \begin{bmatrix}
                \sigma_{11} & \sigma_{11} & \sigma_{12} & \sigma_{12} & \sigma_{12} \\
                \sigma_{11} & \sigma_{11} & \sigma_{12} & \sigma_{12} & \sigma_{12} \\
                \sigma_{12} & \sigma_{12} & \sigma_{22} & \sigma_{22} & \sigma_{22} \\
                \sigma_{12} & \sigma_{12} & \sigma_{22} & \sigma_{22} & \sigma_{22} \\
                \sigma_{12} & \sigma_{12} & \sigma_{22} & \sigma_{22} & \sigma_{22} \\
            \end{bmatrix}


    The matrix is symmetric, so the coefficients in the bottom left corner are 
    the same as the coefficient in the top right corner. The first value in the top left corner of the matrix,
    :math:`\sigma_{11}`, indicates that the :math:`\sigma` value for the interaction
    between the atom 1 and itself is is :math:`\sigma_{11}`. A similar matrix can
    be written for epsilon_sigma_ij.

    Here, the values of the cross coefficients :math:`\sigma_{12}` and :math:`\epsilon_{12}`
    are assumed to follow the arithmetic mean :

    .. math::

        \sigma_{12} = (\sigma_{11}+\sigma_{22})/2 \\
        \epsilon_{12} = (\epsilon_{11}+\epsilon_{22})/2

    Create the following method called *calculate_cross_coefficients* within the 
    *Prepare* class:

    . . label:: start_Prepare_class

    .. code-block:: python

        def calculate_cross_coefficients(self):
            """Calculate all the cross-coefficients for the LJ interations."""
            self.identify_atom_properties() # TOFIX: this was left because of GCMC. Remove? Move?
            matrix_epsilon_ij = []
            matrix_sigma_ij = []
            for i in range(self.total_number_atoms):
                matrix_epsilon_ij.append([])
                matrix_sigma_ij.append([])
                for j in range(self.total_number_atoms):
                    matrix_epsilon_ij[-1].append((self.atoms_epsilon[i]+self.atoms_epsilon[j])/2)
                    matrix_sigma_ij[-1].append((self.atoms_sigma[i]+self.atoms_sigma[j])/2)
            self.matrix_sigma_ij = matrix_sigma_ij
            self.matrix_epsilon_ij = matrix_epsilon_ij

    . . label:: end_Prepare_class

    After calling for the *identify_atom_properties()* method, a double loop
    is performed over all direct coefficients, and the cross coefficients
    are stored within *matrix_sigma_ij* and *matrix_epsilon_ij*, two matrices.

    Finally, let us call the *calculate_cross_coefficients* method from the
    *__init__()* method.

    . . label:: start_Prepare_class

    .. code-block:: python

        def __init__(self,
            (...)
            self.identify_atom_properties()
            self.calculate_cross_coefficients()

    . . label:: end_Prepare_class

Test the code
-------------

Let us test the *Prepare* class to make sure that it does what is expected.
Just like in the previous example, let us initiates a system with 2 atoms of
type 1, and 3 atoms of type 2:

.. label:: start_test_2a_class

.. code-block:: python

    import numpy as np
    from Prepare import Prepare

    prep = Prepare(number_atoms=[2, 3],
        epsilon=[0.2, 0.4], # kcal/mol
        sigma=[3, 4], # A
        atom_mass=[10, 20], # g/mol
        )

    def test_array(result, expected):
        """Test function comparing *result* and *expected*"""
        assert np.array_equal(result, expected), f"Test failed: {result} != {expected}"
        print("Test passed")

    # Make sure the *atoms_epsilon* gives the expected values
    # of 1 1 2 2 2 (in LJ units)
    test_array(prep.atoms_epsilon, np.array([1., 1., 2., 2., 2.]))

.. label:: end_test_2a_class

This test assert that the generated *atoms_epsilon* array is consistent with
its expected value (see the previous paragraphs).
