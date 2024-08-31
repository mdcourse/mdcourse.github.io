.. _chapter8-label:

Monte Carlo insert
==================

Here, a *Monte Carlo* simulation is implemented in the Grand Canonical ensemble,
where the number of atoms in the system fluctuates. The principle of the
simulation resemble the Monte Carlo move from :ref:`chapter6-label`:

- 1) We start from a given intial configuration, and measure the potential
  energy, :math:`E_\text{pot}^\text{initial}`.
- 2) A random number is chosen, and depending on its value, either a new particle
  will tried to be inserted, or an existing particle will tried to be deleted.
- 3) The energy of the system after the insersion/deletion,
  :math:`E_\text{pot}^\text{trial}`, is measured.
- 4) We then decide to keep or reject the move by calculating
  the difference in energy between the trial and the initial configurations:
  :math:`\Delta E = E_\text{pot}^\text{trial} - E_\text{pot}^\text{initial}`.
  
  - If :math:`\Delta E < 0`, then the move is automatically accepted. 
  - If :math:`\Delta E > 0`, then the move is accepted with a probability given
    by the Boltzmann factor :math:`\exp{- \beta \Delta E}`, where
    :math:`\beta = 1 / k_\text{B} T` and :math:`T` is the imposed temperature.

- 5) Steps 1-4 are repeated a large number of times, generating a broad range of
     possible configurations.

Implementation
--------------

Let us add the following method called *monte_carlo_exchange* to the *MonteCarlo* class:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_exchange(self):
        # The first step is to make a copy of the previous state
        # Since atoms numbers are evolving, its also important to store the
        # neighbor, sigma, and epsilon lists
        self.Epot = self.compute_potential()
        initial_positions = copy.deepcopy(self.atoms_positions)
        initial_number_atoms = copy.deepcopy(self.number_atoms)
        initial_neighbor_lists = copy.deepcopy(self.neighbor_lists)
        initial_sigma_lists = copy.deepcopy(self.sigma_ij_list)
        initial_epsilon_lists = copy.deepcopy(self.epsilon_ij_list)
        # Apply a 50-50 probability to insert or delete
        if np.random.random() < 0.5:
            self.monte_carlo_insert()
        else:
            self.monte_carlo_delete()
        # If np.random.random() < self.acceptation_probability, do nothing
        if np.random.random() > self.acceptation_probability:
            # Reject the new position, revert to inital position
            self.neighbor_lists = initial_neighbor_lists
            self.sigma_ij_list = initial_sigma_lists
            self.epsilon_ij_list = initial_epsilon_lists
            self.atoms_positions = initial_positions
            self.number_atoms = initial_number_atoms

.. label:: end_MonteCarlo_class

First, the potential energy is calculated, and the current state of the
simulations, such as positions and atom numbers, are saved. Then, an insertion
try or a deletion trial is made, each with a probability of 0.5. The
*monte_carlo_insert* and *monte_carlo_delete*, that are implemented here below,
are both calculting the *acceptation_probability*. A random number is again selected,
and if that number is larger than the acceptation probability, then the system
revert to its previous position. If the random number is lower than the acceptation
probability, nothing appends, which means that the last trial is accepted.

Then, let us write the *monte_carlo_insert* method:

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
        self.total_number_atoms = np.sum(self.number_atoms)
        self.update_neighbor_lists()
        self.identify_atom_properties()
        self.update_cross_coefficients()
        trial_Epot = self.compute_potential()
        Lambda = self.calculate_Lambda(self.atom_mass[self.inserted_type])
        beta =  1/self.desired_temperature
        Nat = np.sum(self.number_atoms) # NUmber atoms, should it relly be N? of N (type) ?
        Vol = np.prod(self.box_size[:3]) # box volume
        # dimension of 3 is enforced in the power of the Lambda
        self.acceptation_probability = np.min([1, Vol/(Lambda**3*Nat) \
            *np.exp(beta*(self.desired_mu-trial_Epot+self.Epot))])

.. label:: end_MonteCarlo_class

as well as the *monte_carlo_delete* method:

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
            Nat = np.sum(self.number_atoms) # NUmber atoms, should it relly be N? of N (type) ?
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

    class MonteCarlo(Outputs):
        def __init__(self,
                    (...)
                    displace_mc = None,
                    desired_mu = None,
                    inserted_type = 0,
                    desired_temperature = 300,
                    (...)

.. label:: end_MonteCarlo_class

and

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Outputs):
        def __init__(self,
            (...)
            self.displace_mc = displace_mc
            self.desired_mu = desired_mu
            self.inserted_type = inserted_type
            self.desired_temperature = desired_temperature
            (...)

.. label:: end_MonteCarlo_class

Let us non-dimentionalize desired_mu by adding:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def nondimensionalize_units_3(self):
        (...)
            self.displace_mc /= self.reference_distance
        if self.desired_mu is not None:
            self.desired_mu /= self.reference_energy
        (...)

.. label:: end_MonteCarlo_class

Finally, the *monte_carlo_insert_delete()* method must be included in the run:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def run(self):
        (...)
            self.monte_carlo_move()
            self.monte_carlo_exchange()
            self.wrap_in_box()
        (...)

.. label:: end_MonteCarlo_class

We need to calculate Lambda:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def calculate_Lambda(self, mass):
        """Estimate the de Broglie wavelength."""
        # Is it normal that mass is unitless ???
        m = mass/cst.Avogadro*cst.milli  # kg
        kB = cst.Boltzmann  # J/K
        Na = cst.Avogadro
        kB *= Na/cst.calorie/cst.kilo #  kCal/mol/K
        T = self.desired_temperature  # N
        T_K = T*self.reference_energy/kB  # K
        Lambda = cst.h/np.sqrt(2*np.pi*kB*m*T_K)  # m
        Lambda /= cst.angstrom  # Angstrom
        return Lambda / self.reference_distance # dimensionless

.. label:: end_MonteCarlo_class

Test the code
-------------

One can use a similar test as previously. Let us use a displace distance of
0.5 Angstrom, and make 1000 steps.

.. label:: start_test_MonteCarlo_class

.. code-block:: python

    import os
    from MonteCarlo import MonteCarlo

    mc = MonteCarlo(maximum_steps=1000,
        dumping_period=100,
        thermo_period=100,
        desired_mu = -3,
        inserted_type = 0,
        number_atoms=[50],
        epsilon=[0.1], # kcal/mol
        sigma=[3], # A
        atom_mass=[1], # g/mol
        box_dimensions=[20, 20, 20], # A
        data_folder = "outputs/",
        )
    mc.run()

.. label:: end_test_MonteCarlo_class

The evolution of the potential energy as a function of the
number of steps are written in the *Outputs/Epot.dat* file
and can be plotted.
