Monte Carlo displace
====================

Move
----

Let us add a new method callsed *monte_carlo_displacement* to the *MonteCarlo* class. 
The goal of this method is to implement a Monte Carlo move. The energy of the system is
first evaluated, an atom is randomly selected and moved. Then, the energy
of the system is re-measured, and the move is either accepted or rejected, based on the
Metropolis criteria.

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_displacement(self):
        if self.displace_mc is not None:
            initial_Epot = self.compute_potential(output="potential")
            initial_positions = copy.deepcopy(self.atoms_positions)
            atom_id = np.random.randint(self.total_number_atoms)
            # Test a new position
            self.atoms_positions[atom_id] += (np.random.random(3)-0.5)*self.displace_mc
            trial_Epot = self.compute_potential(output="potential")
            acceptation_probability = np.min([1, np.exp(-self.beta*(trial_Epot-initial_Epot))])
            if np.random.random() <= acceptation_probability: # Accept new position
                pass
            else: # Reject new position
                self.atoms_positions = initial_positions

.. label:: end_MonteCarlo_class

Parameters
----------

Let us non-dimentionalized the units, and improve the *__init__* method:

.. label:: start_MonteCarlo_class

.. code-block:: python

    from scipy import constants as cst
    import numpy as np
    import copy
    from Outputs import Outputs

    import warnings
    warnings.filterwarnings('ignore')

    class MonteCarlo(Outputs):
        def __init__(self,
                    maximum_steps,
                    cut_off = 9,
                    displace_mc = None,
                    neighbor = 1,
                    desired_temperature = 300,
                    *args,
                    **kwargs):
            self.maximum_steps = maximum_steps
            self.cut_off = cut_off
            self.displace_mc = displace_mc
            self.neighbor = neighbor
            self.desired_temperature = desired_temperature
            self.beta =  1/self.desired_temperature
            super().__init__(*args, **kwargs)
            self.nondimensionalize_units_3()

.. label:: end_MonteCarlo_class

where:

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

Finally, let us write a *run* method:

.. label:: start_MonteCarlo_class

.. code-block:: python
        
    def run(self):
        """Perform the loop over time."""
        for self.step in range(0, self.maximum_steps+1):
            self.update_neighbor_lists()
            self.monte_carlo_displacement()
            self.wrap_in_box()

.. label:: end_MonteCarlo_class

Output
------

To follow the simulation, let us create a new update function. Add the following
method to the Outputs class:

.. label:: start_Outputs_class

.. code-block:: python

    def update_log_md_mc(self, velocity=True):
        """Update the log file during MD or MC simulations"""
        if self.thermo_period is not None:
            if (self.step % self.thermo_period == 0) \
                    | (self.thermo_period == 0):
                # refresh values
                if velocity:
                    self.calculate_temperature()
                    self.calculate_kinetic_energy()
                    self.calculate_pressure()
                # convert the units
                if velocity:
                    temperature = self.temperature \
                        * self.reference_temperature  # K
                    pressure = self.pressure \
                        * self.reference_pressure  # Atm
                    Ekin = self.Ekin*self.reference_energy  # kcal/mol
                else:
                    temperature = self.desired_temperature \
                        * self.reference_temperature  # K
                    pressure = 0.0
                    Ekin = 0.0
                volume = np.prod(self.box_size[:3]) \
                    *self.reference_distance**3  # A3
                Epot = self.compute_potential(output="potential") \
                    * self.reference_energy  # kcal/mol
                density = self.calculate_density()  # Unitless
                density *= self.reference_mass \
                    /self.reference_distance**3  # g/mol/A3
                Na = cst.Avogadro
                density *= (cst.centi/cst.angstrom)**3 / Na  # g/cm3
                if self.step == 0:
                    characters = '{:<5} {:<5} {:<9} {:<9} {:<9} {:<13} {:<13} {:<13}'
                    print(characters.format(
                        '%s' % ("step"),
                        '%s' % ("N"),
                        '%s' % ("T (K)"),
                        '%s' % ("p (atm)"),
                        '%s' % ("V (A3)"),
                        '%s' % ("Ep (kcal/mol)"),
                        '%s' % ("Ek (kcal/mol)"),
                        '%s' % ("dens (g/cm3)"),
                        ))
                characters = '{:<5} {:<5} {:<9} {:<9} {:<9} {:<13} {:<13} {:<13}'
                print(characters.format(
                    '%s' % (self.step),
                    '%s' % (self.total_number_atoms),
                    '%s' % (f"{temperature:.3}"),
                    '%s' % (f"{pressure:.3}"),
                    '%s' % (f"{volume:.3}"),
                    '%s' % (f"{Epot:.3}"),
                    '%s' % (f"{Ekin:.3}"),
                    '%s' % (f"{density:.3}"),
                    ))
                for output_value, filename in zip([self.total_number_atoms,
                                                   Epot,
                                                   Ekin,
                                                   pressure,
                                                   temperature,
                                                   density,
                                                   volume],
                                                  ["atom_number.dat",
                                                   "Epot.dat",
                                                   "Ekin.dat",
                                                   "pressure.dat",
                                                   "temperature.dat",
                                                   "density.dat",
                                                   "volume.dat"]):
                    self.update_data_file(output_value, filename)

.. label:: end_Outputs_class

As well as

.. label:: start_Outputs_class

.. code-block:: python

    def update_data_file(self, output_value, filename):
        if self.step == 0:
            f = open(self.data_folder + filename, "w")
        else:
            f = open(self.data_folder + filename, "a")
        characters = "%d %.3f %s"
        v = [self.step, output_value]
        f.write(characters % (v[0], v[1], '\n'))
        f.close()

.. label:: end_Outputs_class

Let us call *update_log_md_mc* from the run method of the MonteCarlo class.
Let us add a dump too:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def run(self):
        (...)
        for self.step in range(0, self.maximum_steps+1):
            (...)
            self.wrap_in_box()
            self.update_log_md_mc(velocity=False)
            self.update_dump_file(filename="dump.mc.lammpstrj")

.. label:: end_MonteCarlo_class

To output the density, let us add the following method to the *Utilities* class:

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

.. label:: start_test_MonteCarlo_class

.. code-block:: python

    import os
    from MonteCarlo import MonteCarlo

    mc = MonteCarlo(maximum_steps=1000,
        dumping_period=100,
        thermo_period=100,
        displace_mc = 0.5,
        number_atoms=[50],
        epsilon=[0.1], # kcal/mol
        sigma=[3], # A
        atom_mass=[1], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    mc.run()

.. label:: end_test_MonteCarlo_class

The evolution of the potential energy as a function of the number of steps
are written in the *Outputs/Epot.dat* file and can be plotted. 