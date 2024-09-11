.. _chapter9-label:

Monte Carlo swap
================

A Monte Carlo swap move consists of randomly attempting to exchange the positions of two
particles. Just like in the case of Monte Carlo displace move, the swap is accepted
based on energy criteria. For systems where the dynamics are slow, such as in glasses, a Monte
Carlo swap move can considerably speed up equilibration.

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_swap(self):
        if self.swap_type[0] is not None:
            self.update_neighbor_lists()
            self.update_cross_coefficients()
            if hasattr(self, 'Epot') is False:
                self.Epot = self.compute_potential()
            initial_Epot = self.Epot
            initial_positions = copy.deepcopy(self.atoms_positions)
            # Pick an atom of type one randomly
            atom_id_1 = np.random.randint(self.number_atoms[self.swap_type[0]])
            # Pick an atom of type two randomly
            atom_id_2 = np.random.randint(self.number_atoms[self.swap_type[1]])

            shift_1 = 0
            for N in self.number_atoms[:self.swap_type[0]]:
                shift_1 += N
            shift_2 = 0
            for N in self.number_atoms[:self.swap_type[1]]:
                shift_2 += N
            # attempt to swap the position of the atoms
            position1 = copy.deepcopy(self.atoms_positions[shift_1+atom_id_1])
            position2 = copy.deepcopy(self.atoms_positions[shift_2+atom_id_2])
            self.atoms_positions[shift_2+atom_id_2] = position1
            self.atoms_positions[shift_1+atom_id_1] = position2
            # force the recalculation of neighbor list
            initial_atoms_sigma = self.atoms_sigma
            initial_atoms_epsilon = self.atoms_epsilon
            initial_atoms_mass = self.atoms_mass
            initial_atoms_type = self.atoms_type
            initial_sigma_ij_list = self.sigma_ij_list
            initial_epsilon_ij_list = self.epsilon_ij_list
            initial_neighbor_lists = self.neighbor_lists
            self.update_neighbor_lists(force_update=True)
            self.identify_atom_properties()
            self.update_cross_coefficients(force_update=True)
            # Measure the potential energy of the new configuration
            trial_Epot = self.compute_potential()
            # Evaluate whether the new configuration should be kept or not
            beta =  1/self.desired_temperature
            delta_E = trial_Epot-initial_Epot
            random_number = np.random.random() # random number between 0 and 1
            acceptation_probability = np.min([1, np.exp(-beta*delta_E)])
            if random_number <= acceptation_probability: # Accept new position
                self.Epot = trial_Epot
                self.successful_swap += 1
            else: # Reject new position
                self.atoms_positions = initial_positions # Revert to initial positions
                self.failed_swap += 1
                self.atoms_sigma = initial_atoms_sigma
                self.atoms_epsilon = initial_atoms_epsilon
                self.atoms_mass = initial_atoms_mass
                self.atoms_type = initial_atoms_type
                self.sigma_ij_list = initial_sigma_ij_list
                self.epsilon_ij_list = initial_epsilon_ij_list
                self.neighbor_lists = initial_neighbor_lists
                
.. label:: end_MonteCarlo_class

Let us initialise swap counter:

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Outputs):
        def __init__(self,
            (...)
            self.failed_move = 0
            self.successful_swap = 0
            self.failed_swap = 0

.. label:: end_MonteCarlo_class

Complete the *__init__* method as follows:

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Outputs):
        def __init__(self,
                    (...)
                    displace_mc = None,
                    swap_type = [None, None],

.. label:: end_MonteCarlo_class

and

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Outputs):
        def __init__(self,
            (...)
            self.displace_mc = displace_mc
            self.swap_type = swap_type

.. label:: end_MonteCarlo_class

Finally, the *monte_carlo_exchange()* method must be included in the run:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def run(self):
        (...)
            self.monte_carlo_move()
            self.monte_carlo_swap()

.. label:: end_MonteCarlo_class

Test the code
-------------

Let's test the Monte Carlo swap.

.. label:: start_test_9b_class

.. code-block:: python

    from MonteCarlo import MonteCarlo
    from pint import UnitRegistry
    ureg = UnitRegistry()
    import os

    # Define atom number of each group
    nmb_1 = 50
    nmb_2 = 50  # New group for testing swaps
    # Define LJ parameters (sigma)
    sig_1 = 3 * ureg.angstrom
    sig_2 = 4 * ureg.angstrom  # Different sigma for group 2
    # Define LJ parameters (epsilon)
    eps_1 = 0.1 * ureg.kcal / ureg.mol
    eps_2 = 0.15 * ureg.kcal / ureg.mol  # Different epsilon for group 2
    # Define atom mass
    mss_1 = 10 * ureg.gram / ureg.mol
    mss_2 = 12 * ureg.gram / ureg.mol  # Different mass for group 2
    # Define box size
    L = 20 * ureg.angstrom
    # Define a cut off
    rc = 2.5 * sig_1
    # Pick the desired temperature
    T = 300 * ureg.kelvin

    # Initialize the prepare object
    mc = MonteCarlo(
        ureg=ureg,
        maximum_steps=100,
        thermo_period=10,
        dumping_period=10,
        number_atoms=[nmb_1, nmb_2],  # Include two groups of atoms for swap
        epsilon=[eps_1, eps_2],  # kcal/mol
        sigma=[sig_1, sig_2],  # A
        atom_mass=[mss_1, mss_2],  # g/mol
        box_dimensions=[L, L, L],  # A
        cut_off=rc,
        thermo_outputs="Epot-press",
        desired_temperature=T,  # K
        neighbor=1,
        swap_type=[0, 1]  # Enable Monte Carlo swap between groups 1 and 2
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

    # Test the swap counters
    def test_swap_counters():
        assert mc.successful_swap + mc.failed_swap > 0, \
            "Test failed: No swaps were attempted"
        print("Swap test passed")

    # If the script is run directly, execute the tests
    if __name__ == "__main__":
        import pytest
        # Run pytest programmatically
        pytest.main(["-s", __file__])

.. label:: end_test_9b_class


