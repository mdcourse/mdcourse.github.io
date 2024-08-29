.. _chapter6-label:

Monte Carlo move
================

.. figure:: ../projects/project1/avatar-dm.webp
    :alt: The fluid made of argon atoms and simulated using monte carlo and python.
    :height: 200
    :align: right
    :class: only-dark

.. figure:: ../projects/project1/avatar.webp
    :alt: The fluid made of argon atoms and simulated using monte carlo and python
    :height: 200
    :align: right
    :class: only-light

Here, a *Monte Carlo move* simulation is implemented. The principle of the
simulation is the following:

- 1) We start from a given intial configuration, and measure the potential
  energy, :math:`E_\text{pot}^\text{initial}`.
- 2) One of the particles is picked and moved in a random direction. This displacement
  is made over a distance lower than a certain parameter, :math:`d_\text{mc}`.
- 3) The energy of the system after the move, :math:`E_\text{pot}^\text{trial}`, is measured.
- 4) We then have to decide to keep or reject the move. This is done by calculating
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

Let us add a method named *monte_carlo_move* to the *MonteCarlo* class:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_move(self):
        """Monte Carlo move trial."""
        if self.displace_mc is not None: # only trigger if displace_mc was provided by the user
            try: # try using the last saved Epot, if it exists
                initial_Epot = self.Epot
            except: # If self.Epot does not exists yet, calculate it
                initial_Epot = self.compute_potential(output="potential")
            # Make a copy of the initial atoms positions
            initial_positions = copy.deepcopy(self.atoms_positions)
            # Pick an atom id randomly
            atom_id = np.random.randint(self.total_number_atoms)
            # Move the chosen atom in a random direction
            # The maximum displacement is set by self.displace_mc
            if self.dimensions == 3:
                move = (np.random.random(3)-0.5)*self.displace_mc 
            elif self.dimensions == 2: # the third value will be 0
                move = np.append((np.random.random(2) - 0.5) * self.displace_mc, 0)
            self.atoms_positions[atom_id] += move
            # Measure the optential energy of the new configuration
            trial_Epot = self.compute_potential(output="potential")
            # Evaluate whether the new configuration should be kept or not
            beta =  1/self.desired_temperature
            delta_E = trial_Epot-initial_Epot
            random_number = np.random.random() # random number between 0 and 1
            acceptation_probability = np.min([1, np.exp(-beta*delta_E)])
            if random_number <= acceptation_probability: # Accept new position
                self.Epot = trial_Epot
            else: # Reject new position
                self.Epot = initial_Epot
                self.atoms_positions = initial_positions # Revert to initial positions

.. label:: end_MonteCarlo_class

Parameters
----------

The *monte_carlo_move* method requires a few parameters to be selected by the
users, such as *displace_mc* (:math:`d_\text{mc}`), the maximum number of steps,
and the desired temperature (:math:`T`). Let us add these parameters to the
*__init__* method:

.. label:: start_MonteCarlo_class

.. code-block:: python

    class MonteCarlo(Measurements):
        def __init__(self,
                    maximum_steps,
                    cut_off = 9,
                    displace_mc = None,
                    neighbor = 1,
                    desired_temperature = 300,
                    thermo_outputs = "press",
                    data_folder = None,
                    *args,
                    **kwargs):
            self.maximum_steps = maximum_steps
            self.cut_off = cut_off
            self.displace_mc = displace_mc
            self.neighbor = neighbor
            self.desired_temperature = desired_temperature
            self.thermo_outputs = thermo_outputs
            self.data_folder = data_folder
            if self.data_folder is not None:
                if os.path.exists(self.data_folder) is False:
                    os.mkdir(self.data_folder)
            super().__init__(*args, **kwargs)
            self.nondimensionalize_units_3()

.. label:: end_MonteCarlo_class

Here, we anticipate that some of the parameters have to be nondimensionalized, which
is done with the *nondimensionalize_units_3* method that must also be added to
the *MonteCarlo* class:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def nondimensionalize_units_3(self):
        """Use LJ prefactors to convert units into non-dimensional."""
        self.cut_off = self.cut_off/self.reference_distance
        self.desired_temperature = self.desired_temperature \
            /self.reference_temperature
        if self.displace_mc is not None:
            self.displace_mc /= self.reference_distance

.. label:: end_MonteCarlo_class

Run method
----------

Finally, let us add a *run* method to the *MonteCarlo* class, that is used to
perform a loop over the desired number of steps *maximum_steps*:

.. label:: start_MonteCarlo_class

.. code-block:: python
        
    def run(self):
        """Perform the loop over time."""
        for self.step in range(0, self.maximum_steps+1):
            self.update_neighbor_lists()
            self.update_cross_coefficients()
            self.monte_carlo_move()
            self.wrap_in_box()

.. label:: end_MonteCarlo_class

At each step, the *monte_carlo_move* method is called. The previously defined
methods *update_neighbor_lists* and *wrap_in_box* are also called to ensure that
the neighbor lists are kept up to date despite the motion of the atoms, and that
the atoms remain inside the box, respectively.

Let us call *update_log_md_mc* from the run method of the MonteCarlo class.
Let us add a dump too:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def run(self):
        (...)
        for self.step in range(0, self.maximum_steps+1):
            (...)
            self.wrap_in_box()
            log_simulation_data(self)
            update_dump_file(self, "dump.mc.lammpstrj")

.. label:: end_MonteCarlo_class

To output the density, let us add the following method to the *Utilities* class:
# TOFIX: not used yet

.. label:: start_Utilities_class

.. code-block:: python

    def calculate_density(self):
        """Calculate the mass density."""
        volume = np.prod(self.box_size[:3])  # Unitless
        total_mass = np.sum(self.atoms_mass)  # Unitless
        return total_mass/volume  # Unitless

.. label:: end_Utilities_class

Test the code
-------------

One can use a similar test as previously. Let us use a displace distance of
0.5 Angstrom, and make 1000 steps.

.. label:: start_test_6a_class

.. code-block:: python

    import os
    from MonteCarlo import MonteCarlo

    mc = MonteCarlo(maximum_steps=1000,
        dumping_period=100,
        thermo_period=100,
        thermo_outputs = "Epot",
        displace_mc = 0.5,
        number_atoms=[50],
        epsilon=[0.1], # kcal/mol
        sigma=[3], # A
        atom_mass=[10], # g/mol
        box_dimensions=[20, 20, 20], # A
        data_folder="Outputs/",
        )
    mc.run()

.. label:: end_test_6a_class

The evolution of the potential energy as a function of the number of steps
are written in the *simulation.log* file. The data can be used to plot
the evolution of the system with time.
