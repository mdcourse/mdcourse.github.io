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
