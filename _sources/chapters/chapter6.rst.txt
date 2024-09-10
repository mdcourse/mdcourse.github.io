.. _chapter6-label:

Monte Carlo dispace
===================

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

Here, a Monte Carlo simulation is implemented where the only allowed move
is a dispacement of the particles. The principle of such Monte Carlo
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
            # When needed, recalculate neighbor/coeff lists
            self.update_neighbor_lists()
            self.update_cross_coefficients()
            # If self.Epot does not exist yet, calculate it
            # It should only be necessary when step = 0
            if hasattr(self, 'Epot') is False:
                self.Epot = self.compute_potential()
            # Make a copy of the initial atom positions and initial energy
            initial_Epot = self.Epot
            initial_positions = copy.deepcopy(self.atoms_positions)
            # Pick an atom id randomly
            atom_id = np.random.randint(np.sum(self.number_atoms))
            # Move the chosen atom in a random direction
            # The maximum displacement is set by self.displace_mc
            move = (np.random.random(3)-0.5)*self.displace_mc 
            self.atoms_positions[atom_id] += move
            # Measure the potential energy of the new configuration
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

The counters *successful_move* and *failed_move* are incremented with each
successful and failed attempt, respectively.

Parameters
----------

The *monte_carlo_move* method requires a few parameters to be selected by
the users, such as *displace_mc* (:math:`d_\text{mc}`, in Ångströms), the
maximum number of steps, and the desired temperature (:math:`T`, in Kelvin).
Let us add these parameters to the *__init__()* method:

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

Run Method
----------

Finally, let us add a *run* method to the *MonteCarlo* class that performs
a loop over the desired number of steps, *maximum_steps*:

.. label:: start_MonteCarlo_class

.. code-block:: python
        
    def run(self):
        """Perform the loop over time."""
        for self.step in range(0, self.maximum_steps+1):
            self.monte_carlo_move()
            self.wrap_in_box()

.. label:: end_MonteCarlo_class

At each step, the *monte_carlo_move()* method is called. The previously
defined *wrap_in_box* method is also called to ensure that the atoms remain
inside the box. Additionally, let us call *log_simulation_data()* and
*update_dump_file()*:

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

Test the Code
-------------

Let us use a similar test as before. Set a displacement distance corresponding
to a quarter of sigma, and perform a very small number of steps:

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
        maximum_steps=100,
        thermo_period=10,
        dumping_period=10,
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

.. label:: end_test_6a_class

The evolution of the potential energy as a function of the number of steps
is recorded in the *simulation.log* file. The data in *simulation.log* can
be used to plot the evolution of the system over time.
