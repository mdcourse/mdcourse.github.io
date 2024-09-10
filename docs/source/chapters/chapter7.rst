.. _chapter7-label:

Pressure Measurement
====================

To extract the equation of state in our simulation, we need to measure the
pressure of the system, :math:`p`. The pressure in a molecular simulation can
be calculated from the interactions between particles. The pressure can be
measured as the sum of the ideal contribution,
:math:`p_\text{ideal} = N_\text{DOF} k_\text{B} T / V d`, 
which comes from the ideal gas law, and a Virial term that
accounts for the pressure contribution from the forces between particles,
:math:`p_\text{non_ideal} = \left< \sum_i r_i \cdot F_i \right> / V d`. The final
expression reads :cite:`frenkel2023understanding`:

.. math::

    p = \dfrac{1}{V d} \left[ N_\text{DOF} k_\text{B} T +  \left< \sum_i r_i \cdot F_i \right> \right]

:math:`N_\text{DOF}` is the number of degrees-of-freedom, which can be calculated
from the number of particles, :math:`N`, and the dimension of the system,
:math:`d = 3`, as :math:`N_\text{DOF} = d N - d` :cite:`frenkel2023understanding`.

The calculation of :math:`p_\text{ideal}` is straightforward. For Monte Carlo
simulation, as atoms do not have a temperature, the *imposed* temperature will
be used as ":math:`T`". The calculation of :math:`p_\text{non_ideal}` requires the
measurement of all the forces and distances between atoms. The calculation of
the forces was already implemented in a previous chapter (see *compute_force*
the :ref:`chapter4-label` chapter), but a new function that returns all the
vector directions between atom pairs will need to be written here.

Implement the Virial equation
-----------------------------

Let us add the following method to the *Measurements* class:

.. label:: start_Measurements_class

.. code-block:: python

    def calculate_pressure(self):
        """Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smit,
        Understanding molecular simulation: from algorithms to applications, 2002)"""
        # Compute the ideal contribution
        Nat = np.sum(self.number_atoms) # total number of atoms
        dimension = 3 # 3D is the only possibility here
        Ndof = dimension*Nat-dimension    
        volume = np.prod(self.box_size[:3]) # box volume
        try:
            self.calculate_temperature() # this is for later on, when velocities are computed
            temperature = self.temperature
        except:
            temperature = self.desired_temperature # for MC, simply use the desired temperature
        p_ideal = Ndof*temperature/(volume*dimension)
        # Compute the non-ideal contribution
        distances_forces = np.sum(self.compute_force(return_vector = False) \
                                  *self.evaluate_rij_matrix())
        p_nonideal = distances_forces/(volume*dimension)
        # Final pressure
        self.pressure = p_ideal+p_nonideal

.. label:: end_Measurements_class

To evaluate all the vectors between all the particles, let us also add the
*evaluate_rij_matrix* method to the *Utilities* class:

.. label:: start_Utilities_class

.. code-block:: python

    def evaluate_rij_matrix(self):
        """Matrix of vectors between all particles."""
        Nat = np.sum(self.number_atoms)
        Box = self.box_size[:3]
        rij_matrix = np.zeros((Nat, Nat,3))
        pos_j = self.atoms_positions
        for Ni in range(Nat-1):
            pos_i = self.atoms_positions[Ni]
            rij_xyz = (np.remainder(pos_i - pos_j + Box/2.0, Box) - Box/2.0)
            rij_matrix[Ni] = rij_xyz
        return rij_matrix

.. label:: end_Utilities_class

Test the code
-------------

Let us test the outputted pressure. 

.. label:: start_test_7a_class

.. code-block:: python

    from MonteCarlo import MonteCarlo
    from pint import UnitRegistry
    ureg = UnitRegistry()
    import os

    # Define atom number of each group
    nmb_1= 50
    # Define LJ parameters (sigma)
    sig_1 = 3*ureg.angstrom
    # Define LJ parameters (epsilon)
    eps_1 = 0.1*ureg.kcal/ureg.mol
    # Define atom mass
    mss_1 = 10*ureg.gram/ureg.mol
    # Define box size
    L = 20*ureg.angstrom
    # Define a cut off
    rc = 2.5*sig_1
    # Pick the desired temperature
    T = 300*ureg.kelvin
    # choose the displace_mc
    displace_mc = sig_1/4

    # Initialize the prepare object
    mc = MonteCarlo(
        ureg = ureg,
        maximum_steps=100,
        thermo_period=10,
        dumping_period=10,
        number_atoms=[nmb_1],
        epsilon=[eps_1], # kcal/mol
        sigma=[sig_1], # A
        atom_mass=[mss_1], # g/mol
        box_dimensions=[L, L, L], # A
        cut_off=rc,
        thermo_outputs="Epot-press",
        desired_temperature=T, # K
        neighbor=20,
        displace_mc = displace_mc,
    )

    # Run the Monte Carlo simulation
    mc.run()

    # Test function using pytest
    def test_output_files():
        assert os.path.exists("Outputs/dump.mc.lammpstrj"), \
        "Test failed: dump file was not created"
        assert os.path.exists("Outputs/simulation.log"), \
        "Test failed: log file was not created"
        print("Test passed")

    # If the script is run directly, execute the tests
    if __name__ == "__main__":
        import pytest
        # Run pytest programmatically
        pytest.main(["-s", __file__])

.. label:: end_test_7a_class

The pressure should be returned alongside the potential energy within *simulation.log*:

.. code-block:: bw

    step Epot press
    0 134248.72 4608379.27
    10 124905.76 4287242.75
    20 124858.91 4285584.40
    (...)