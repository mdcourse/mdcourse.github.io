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
            # If needed, recalculate neighbor/coeff lists
            self.update_neighbor_lists()
            self.update_cross_coefficients()
            if hasattr(self, 'Epot') is False: # If self.Epot does not exists yet, calculate it
                self.Epot = self.compute_potential()
            initial_Epot = self.Epot
            # Make a copy of the initial atoms positions
            initial_positions = copy.deepcopy(self.atoms_positions)
            # Pick an atom id randomly
            atom_id = np.random.randint(np.sum(self.number_atoms))
            # Move the chosen atom in a random direction
            # The maximum displacement is set by self.displace_mc
            move = (np.random.random(3)-0.5)*self.displace_mc 
            self.atoms_positions[atom_id] += move
            # Measure the optential energy of the new configuration
            trial_Epot = self.compute_potential()
            # Evaluate whether the new configuration should be kept or not
            beta =  1/self.desired_temperature
            delta_E = trial_Epot-initial_Epot
            random_number = np.random.random() # random number between 0 and 1
            acceptation_probability = np.min([1, np.exp(-beta*delta_E)])
            if random_number <= acceptation_probability: # Accept new position
                self.Epot = trial_Epot
                self.successful_move += 1
            else: # Reject new position
                self.atoms_positions = initial_positions # Revert to initial positions
                self.failed_move += 1

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
                    desired_temperature,
                    displace_mc = None,
                    *args,
                    **kwargs):
            self.maximum_steps = maximum_steps
            self.displace_mc = displace_mc
            self.desired_temperature = desired_temperature
            super().__init__(*args, **kwargs)
            self.nondimensionalize_units(["desired_temperature", "displace_mc"])
            self.successful_move = 0
            self.failed_move = 0

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
            self.monte_carlo_move()
            self.wrap_in_box()

.. label:: end_MonteCarlo_class

At each step, the *monte_carlo_move* method is called. The previously defined
mthe *wrap_in_box* method is also called to ensure that
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
        maximum_steps=1000,
        thermo_period=100,
        dumping_period=100,
        number_atoms=[nmb_1],
        epsilon=[eps_1], # kcal/mol
        sigma=[sig_1], # A
        atom_mass=[mss_1], # g/mol
        box_dimensions=[L, L, L], # A
        cut_off=rc,
        thermo_outputs="Epot",
        desired_temperature=T, # K
        neighbor=20,
        displace_mc = displace_mc,
    )
    mc.run()

.. label:: end_test_6a_class

The evolution of the potential energy as a function of the number of steps
are written in the *simulation.log* file. The data can be used to plot
the evolution of the system with time.
