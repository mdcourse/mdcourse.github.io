Monte Carlo simulation
======================

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
                    neighbor = 10,
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
            # self.update_log_md_mc(velocity=False) --> we will need to add it too
            self.update_dump_file(filename="dump.mc.lammpstrj")

.. label:: end_MonteCarlo_class

Test the code
-------------

One can use the same test as previously, and ask the code to print information
every 10 steps in the dump files, as well as in the log:

.. label:: start_test_MonteCarlo_class

.. code-block:: python

    import os
    from MonteCarlo import MonteCarlo

    mc = MonteCarlo(maximum_steps=100,
        dumping_period=10,
        displace_mc = 0.5,
        number_atoms=[30],
        epsilon=[0.1], # kcal/mol
        sigma=[3], # A
        atom_mass=[1], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    mc.run()

.. label:: end_test_MonteCarlo_class
