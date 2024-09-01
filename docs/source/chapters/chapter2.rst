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
- Density is in g/cm\ :sup:`3`.

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

Let's fill in the previously created class named *Prepare*. To facilitate unit
conversion, we will import NumPy and the constants module from |SciPy|.

.. |SciPy| raw:: html

   <a href="https://scipy.org/" target="_blank">SciPy</a>

In the file named *Prepare.py*, add the following lines:

.. label:: start_Prepare_class

.. code-block:: python

    import numpy as np
    from scipy import constants as cst

.. label:: end_Prepare_class

Four atom parameters are provided to the *Prepare* class:

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
                    ureg, # Pint unit registry
                    number_atoms, # List - no unit
                    epsilon, # List - Kcal/mol
                    sigma, # List - Angstrom
                    atom_mass,  # List - g/mol
                    *args,
                    **kwargs):
            self.ureg = ureg
            self.number_atoms = number_atoms
            self.epsilon = epsilon
            self.sigma = sigma
            self.atom_mass = atom_mass
            super().__init__(*args, **kwargs)

.. label:: end_Prepare_class

Here, *number_atoms* :math:`N`, *epsilon* :math:`\epsilon`,
*sigma* :math:`\sigma`, and *atom_mass* :math:`m` must be provided as lists
where the elements have no units, kcal/mol, angstrom, and g/mol units,
respectively. The units will be enforced with the |Pint| unit registry, *ureg*,
which must also be provided as a parameter (more details later).

.. |Pint| raw:: html

   <a href="https://pint.readthedocs.io" target="_blank">Pint</a>

All the parameters are assigned to *self*, allowing other methods to access
them. The *args* and *kwargs* parameters are used to accept an arbitrary number
of positional and keyword arguments, respectively.

Calculate LJ units prefactors
-----------------------------

Within the *Prepare* class, create a method called *calculate_LJunits_prefactors*
that will be used to calculate the prefactors necessary to convert units from the
*real* unit system to the *LJ* unit system:

.. label:: start_Prepare_class

.. code-block:: python

    def calculate_LJunits_prefactors(self):
        """Calculate the Lennard-Jones units prefactors."""
        # First define Boltzmann and Avogadro constants
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo  # kcal/mol/K
        kB *= self.ureg.kcal/self.ureg.mol/self.ureg.kelvin
        Na = cst.Avogadro/self.ureg.mol
        # Define the reference distance, energy, and mass
        self.ref_length = self.sigma[0]  # Angstrom
        self.ref_energy = self.epsilon[0]  # kcal/mol
        self.ref_mass = self.atom_mass[0]  # g/mol
        # Optional: assert that units were correctly provided by users
        assert self.ref_length.units == self.ureg.angstrom, \
            f"Error: Provided sigma has wrong units, should be angstrom"
        assert self.ref_energy.units == self.ureg.kcal/self.ureg.mol, \
            f"Error: Provided epsilon has wrong units, should be kcal/mol"
        assert self.ref_mass.units == self.ureg.g/self.ureg.mol, \
            f"Error: Provided mass has wrong units, should be g/mol"
        # Calculate the prefactor for the time (in femtosecond)
        self.ref_time = np.sqrt(self.ref_mass \
            *self.ref_length**2/self.ref_energy).to(self.ureg.femtosecond)
        # Calculate the prefactor for the temperature (in Kelvin)
        self.ref_temperature = self.ref_energy/kB  # Kelvin
        # Calculate the prefactor for the pressure (in Atmosphere)
        self.ref_pressure = (self.ref_energy \
            /self.ref_length**3/Na).to(self.ureg.atmosphere)
        # Group all the reference quantities into a list for practicality
        self.ref_quantities = [self.ref_length, self.ref_energy,
            self.ref_mass, self.ref_time, self.ref_pressure, self.ref_temperature]
        self.ref_units = [ref.units for ref in self.ref_quantities]

.. label:: end_Prepare_class

This method defines the reference length as the first element in the
*sigma* list, i.e., :math:`\sigma_{11}`. Therefore, atoms of type 1 will
always be used for normalization. Similarly, the first element in the
*epsilon* list (:math:`\epsilon_{11}`) is used as the reference energy, and
the first element in the *atom_mass* list (:math:`m_1`) is used as the
reference mass.

The reference time in femtoseconds is then calculated as
:math:`\sqrt{m_1 \sigma_{11}^2 / \epsilon_{11}}`, the reference
temperature in Kelvin as :math:`\epsilon_{11} / k_\text{B}`, and the
reference pressure in atmospheres as :math:`\epsilon_{11} / \sigma_{11}^3`.

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

Create a new method called *nondimensionalize_units* within the *Prepare*
class:

.. label:: start_Prepare_class

.. code-block:: python

    def nondimensionalize_units(self, quantities_to_normalise):
        for name in quantities_to_normalise:
            quantity = getattr(self, name)  # Get the attribute by name
            if isinstance(quantity, list):
                for i, element in enumerate(quantity):
                    assert element.units in self.ref_units, \
                        f"Error: Units not part of the reference units"
                    ref_value = self.ref_quantities[self.ref_units.index(element.units)]
                    quantity[i] = element/ref_value
                    assert quantity[i].units == self.ureg.dimensionless, \
                        f"Error: Quantities are not properly nondimensionalized"
                    quantity[i] = quantity[i].magnitude # get rid of ureg
                setattr(self, name, quantity)
            elif len(np.shape(quantity)) > 0: # for position array
                assert element.units in self.ref_units, \
                    f"Error: Units not part of the reference units"
                ref_value = self.ref_quantities[self.ref_units.index(element.units)]
                quantity = quantity/ref_value
                assert quantity.units == self.ureg.dimensionless, \
                    f"Error: Quantities are not properly nondimensionalized"
                quantity = quantity.magnitude # get rid of ureg
                setattr(self, name, quantity)
            else:
                if quantity is not None:
                    assert np.shape(quantity) == (), \
                        f"Error: The quantity is a list or an array"
                    assert quantity.units in self.ref_units, \
                        f"Error: Units not part of the reference units"
                    ref_value = self.ref_quantities[self.ref_units.index(quantity.units)]
                    quantity = quantity/ref_value
                    assert quantity.units == self.ureg.dimensionless, \
                        f"Error: Quantities are not properly nondimensionalized"
                    quantity = quantity.magnitude # get rid of ureg
                    setattr(self, name, quantity)

.. label:: end_Prepare_class

When a *quantities_to_normalise* list containing parameter names is provided
to the *nondimensionalize_units* method, a loop is performed over all the
quantities. The value of each quantity is extracted using *getattr*. The units
of the quantity of interest are then detected and normalized by the appropriate
reference quantities defined by *calculate_LJunits_prefactors*.

Let us also call the *nondimensionalize_units* from the *__init__()* method
of the *Prepare* class:

.. label:: start_Prepare_class

.. code-block:: python

    def __init__(self,
        (...)
        self.calculate_LJunits_prefactors()
        self.nondimensionalize_units(["epsilon", "sigma", "atom_mass"])

.. label:: end_Prepare_class

Here, the *epsilon*, *sigma*, and *atom_mass* parameters will be
nondimensionalized.

Identify Atom Properties
------------------------

Anticipating the future use of multiple atom types, where each type will be
associated with its own :math:`\sigma`, :math:`\epsilon`, and :math:`m`, let
us create arrays containing the properties of each atom in the simulation. For
instance, in a simulation with two atoms of type 1 and three atoms of type 2,
the corresponding *atoms_sigma* array will be:

.. math::

    \text{atoms_sigma} = [\sigma_{11}, \sigma_{11}, \sigma_{22}, \sigma_{22}, \sigma_{22}]

where :math:`\sigma_{11}` and :math:`\sigma_{22}` are the sigma values for
atoms of type 1 and 2, respectively. The *atoms_sigma* array will facilitate
future calculations of forces.

Create a new method called *identify_atom_properties*, and place it
within the *Prepare* class:

.. label:: start_Prepare_class

.. code-block:: python

    def identify_atom_properties(self):
        """Identify the properties for each atom."""
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
        self.nondimensionalize_units(["epsilon", "sigma", "atom_mass"])
        self.identify_atom_properties()

.. label:: end_Prepare_class

Test the code
-------------

Let's test the *Prepare* class to make sure that it does what is expected.
Here, a system containing 2 atoms of type 1, and 3 atoms of type 2 is
prepared. LJs parameters and masses for each groups are also defined, and given
physical units thanks to the *UnitRegistry* of Pint.

.. label:: start_test_2a_class

.. code-block:: python

    import numpy as np
    from Prepare import Prepare
    from pint import UnitRegistry
    ureg = UnitRegistry()

    # Define atom number of each group
    nmb_1, nmb_2= [2, 3]
    # Define LJ parameters (sigma)
    sig_1, sig_2 = [3, 4]*ureg.angstrom
    # Define LJ parameters (epsilon)
    eps_1, eps_2 = [0.2, 0.4]*ureg.kcal/ureg.mol
    # Define atom mass
    mss_1, mss_2 = [10, 20]*ureg.gram/ureg.mol

    # Initialize the prepare object
    prep = Prepare(
        ureg = ureg,
        number_atoms=[nmb_1, nmb_2],
        epsilon=[eps_1, eps_2], # kcal/mol
        sigma=[sig_1, sig_2], # A
        atom_mass=[mss_1, mss_2], # g/mol
    )

    # Test function using pytest
    def test_atoms_epsilon():
        expected = np.array([1., 1., 2., 2., 2.])
        result = prep.atoms_epsilon
        assert np.array_equal(result, expected), \
        f"Test failed: {result} != {expected}"
        print("Test passed")

    # In the script is launched with Python, call Pytest
    if __name__ == "__main__":
        import pytest
        pytest.main(["-s", __file__])

.. label:: end_test_2a_class

This test asserts that the generated *atoms_epsilon* array is consistent
with its expected value (see the previous paragraphs). Note that the expected
values are in LJ units, i.e., they have been nondimensionalized by :math:`\epsilon_1`.
