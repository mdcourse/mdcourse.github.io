.. _chapter8-label:

Monte Carlo insert
==================

Here, a *Monte Carlo* simulation is implemented in the Grand Canonical ensemble,
where the number of atoms in the system fluctuates. The principle of the
simulation resembles the Monte Carlo move from :ref:`chapter6-label`:

- 1) We start from a given initial configuration, and measure the potential
  energy, :math:`E_\text{pot}^\text{initial}`.
- 2) A random number is chosen, and depending on its value, either a new particle
  will try to be inserted, or an existing particle will try to be deleted.
- 3) The energy of the system after the insertion/deletion,
  :math:`E_\text{pot}^\text{trial}`, is measured.
- 4) We then decide to keep or reject the move by calculating
  the difference in energy between the trial and the initial configurations:
  :math:`\Delta E = E_\text{pot}^\text{trial} - E_\text{pot}^\text{initial}`.
- 5) Steps 1-4 are repeated a large number of times, generating a broad range of
     possible configurations.

Implementation
--------------

Let us add the following method called *monte_carlo_exchange* to the *MonteCarlo* class:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_exchange(self):
        if self.desired_mu is not None:
            # The first step is to make a copy of the previous state
            # Since atoms numbers are evolving, its also important to store the
            # neighbor, sigma, and epsilon lists
            self.Epot = self.compute_potential() # TOFIX: not necessary every time
            initial_positions = copy.deepcopy(self.atoms_positions)
            initial_number_atoms = copy.deepcopy(self.number_atoms)
            initial_neighbor_lists = copy.deepcopy(self.neighbor_lists)
            initial_sigma_lists = copy.deepcopy(self.sigma_ij_list)
            initial_epsilon_lists = copy.deepcopy(self.epsilon_ij_list)
            # Apply a 50-50 probability to insert or delete
            insert_or_delete = np.random.random()
            if np.random.random() < insert_or_delete:
                self.monte_carlo_insert()
            else:
                self.monte_carlo_delete()
            if np.random.random() < self.acceptation_probability: # accepted move
                # Update the success counters
                if np.random.random() < insert_or_delete:
                    self.successful_insert += 1
                else:
                    self.successful_delete += 1
            else:
                # Reject the new position, revert to inital position
                self.neighbor_lists = initial_neighbor_lists
                self.sigma_ij_list = initial_sigma_lists
                self.epsilon_ij_list = initial_epsilon_lists
                self.atoms_positions = initial_positions
                self.number_atoms = initial_number_atoms
                # Update the failed counters
                if np.random.random() < insert_or_delete:
                    self.failed_insert += 1
                else:
                    self.failed_delete += 1

.. label:: end_MonteCarlo_class

First, the potential energy is calculated, and the current state of the
simulation, such as positions and atom numbers, is saved. Then, an insertion
trial or a deletion trial is made, each with a probability of 0.5. The
*monte_carlo_insert* and *monte_carlo_delete* methods, which are implemented
below, both calculate the *acceptance_probability*. A random number is selected,
and if that number is larger than the acceptance probability, the system reverts
to its previous position. If the random number is lower than the acceptance
probability, nothing happens, meaning that the last trial is accepted.

Then, let us write the *monte_carlo_insert()* method:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_insert(self):
        self.number_atoms[self.inserted_type] += 1
        new_atom_position = np.random.random(3)*np.diff(self.box_boundaries).T \
            - np.diff(self.box_boundaries).T/2
        shift_id = 0 
        for N in self.number_atoms[:self.inserted_type]:
            shift_id += N
        self.atoms_positions = np.vstack([self.atoms_positions[:shift_id],
                                        new_atom_position,
                                        self.atoms_positions[shift_id:]])
        self.update_neighbor_lists()
        self.identify_atom_properties()
        self.update_cross_coefficients()
        trial_Epot = self.compute_potential()
        Lambda = self.calculate_Lambda(self.atom_mass[self.inserted_type])
        beta =  1/self.desired_temperature
        Nat = np.sum(self.number_atoms) # Number atoms, should it really be N? of N (type) ?
        Vol = np.prod(self.box_size[:3]) # box volume
        # dimension of 3 is enforced in the power of the Lambda
        self.acceptation_probability = np.min([1, Vol/(Lambda**3*Nat) \
            *np.exp(beta*(self.desired_mu-trial_Epot+self.Epot))])

.. label:: end_MonteCarlo_class

After trying to insert a new particle, neighbor lists and cross coefficients
must be re-evaluated. Then, the acceptance probability is calculated.

Let us add the very similar *monte_carlo_delete()* method:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_delete(self):
        # Pick one atom to delete randomly
        atom_id = np.random.randint(self.number_atoms[self.inserted_type])
        self.number_atoms[self.inserted_type] -= 1
        if self.number_atoms[self.inserted_type] > 0:
            shift_id = 0
            for N in self.number_atoms[:self.inserted_type]:
                shift_id += N
            self.atoms_positions = np.delete(self.atoms_positions, shift_id+atom_id, axis=0)
            self.update_neighbor_lists()
            self.identify_atom_properties()
            self.update_cross_coefficients()
            trial_Epot = self.compute_potential()
            Lambda = self.calculate_Lambda(self.atom_mass[self.inserted_type])
            beta =  1/self.desired_temperature
            Nat = np.sum(self.number_atoms) # Number atoms, should it really be N? of N (type) ?
            Vol = np.prod(self.box_size[:3]) # box volume
            # dimension of 3 is enforced in the power of the Lambda
            self.acceptation_probability = np.min([1, (Lambda**3 *(Nat-1)/Vol) \
                *np.exp(-beta*(self.desired_mu+trial_Epot-self.Epot))])
        else:
            print("Error: no more atoms to delete")

.. label:: end_MonteCarlo_class

Complete the *__init__* method as follows:

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Measurements):
        def __init__(self,
                    (...)
                    displace_mc = None,
                    desired_mu = None,
                    inserted_type = 0,

.. label:: end_MonteCarlo_class

and

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Measurements):
        def __init__(self,
            (...)
            self.displace_mc = displace_mc
            self.desired_mu = desired_mu
            self.inserted_type = inserted_type

.. label:: end_MonteCarlo_class

Let us also normalize the "desired_mu":

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Outputs):
        def __init__(self,
            (...)
            self.nondimensionalize_units(["desired_temperature", "displace_mc"])
            self.nondimensionalize_units(["desired_mu"])
            self.successful_move = 0
            self.failed_move = 0
            self.successful_insert = 0
            self.failed_insert = 0
            self.successful_delete = 0
            self.failed_delete = 0

.. label:: end_MonteCarlo_class

Finally, the *monte_carlo_exchange()* method must be included in the run:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def run(self):
        (...)
            self.monte_carlo_move()
            self.monte_carlo_exchange()
            self.wrap_in_box()
        (...)

.. label:: end_MonteCarlo_class

We need to calculate :math:`\Lambda`:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def calculate_Lambda(self, mass):
        """Estimate the de Broglie wavelength."""
        T = self.desired_temperature  # N
        return 1/np.sqrt(2*np.pi*mass*T)

.. label:: end_MonteCarlo_class

To output the density, let us add the following method to the *Measurements* class:

.. label:: start_Measurements_class

.. code-block:: python

    def calculate_density(self):
        """Calculate the mass density."""
        # TOFIX: not used yet
        volume = np.prod(self.box_size[:3])  # Unitless
        total_mass = np.sum(self.atoms_mass)  # Unitless
        return total_mass/volume  # Unitless

.. label:: end_Measurements_class

Test the code
-------------

One can use a similar test as previously, but with an imposed chemical
potential *desired_mu*:

.. label:: start_test_8a_class

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
    # choose the desired_mu
    desired_mu = -3*ureg.kcal/ureg.mol

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
        neighbor=1,
        desired_mu = desired_mu,
    )
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

.. label:: end_test_98a_class

The evolution of the potential energy as a function of the
number of steps is written in the *Outputs/Epot.dat*.
