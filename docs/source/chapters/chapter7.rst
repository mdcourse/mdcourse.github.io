Grand Canonical Monte Carlo
===========================

In the Grand Canonical ensemble, the number of atom in the
simulation box is not constant. 

Let us add the following method to the MonteCarlo class:

.. label:: start_MonteCarlo_class

.. code-block:: python

    def monte_carlo_insert_delete(self):
        if self.desired_mu is not None:
            initial_Epot = self.compute_potential(output="potential")
            initial_positions = copy.deepcopy(self.atoms_positions)
            initial_number_atoms = copy.deepcopy(self.number_atoms)
            initial_neighbor_lists = copy.deepcopy(self.neighbor_lists)
            volume = np.prod(np.diff(self.box_boundaries))
            if np.random.random() < 0.5:
                # Try adding an atom
                self.number_atoms[self.inserted_type] += 1
                atom_position = np.zeros((1, self.dimensions))
                for dim in np.arange(self.dimensions):
                    atom_position[:, dim] = np.random.random(1)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2
                shift_id = 0
                for N in self.number_atoms[:self.inserted_type]:
                    shift_id += N
                self.atoms_positions = np.vstack([self.atoms_positions[:shift_id], atom_position, self.atoms_positions[shift_id:]])
                self.update_neighbor_lists()
                self.calculate_cross_coefficients()
                trial_Epot = self.compute_potential(output="potential")
                Lambda = self.calculate_Lambda(self.atom_mass[self.inserted_type])
                acceptation_probability = np.min([1, volume/(Lambda**self.dimensions*(self.total_number_atoms))*np.exp(self.beta*(self.desired_mu-trial_Epot+initial_Epot))])
            else:
                # Pick one atom to delete randomly
                atom_id = np.random.randint(self.number_atoms[self.inserted_type])
                self.number_atoms[self.inserted_type] -= 1
                if self.number_atoms[self.inserted_type] > 0:
                    shift_id = 0
                    for N in self.number_atoms[:self.inserted_type]:
                        shift_id += N
                    self.atoms_positions = np.delete(self.atoms_positions, shift_id+atom_id, axis=0)
                    self.update_neighbor_lists()
                    self.calculate_cross_coefficients()
                    trial_Epot = self.compute_potential(output="potential")
                    Lambda = self.calculate_Lambda(self.atom_mass[self.inserted_type])
                    acceptation_probability = np.min([1, (Lambda**self.dimensions*(self.total_number_atoms-1)/volume)*np.exp(-self.beta*(self.desired_mu+trial_Epot-initial_Epot))])
                else:
                    acceptation_probability = 0
            if np.random.random() < acceptation_probability:
                # Accept the new position
                pass
            else:
                # Reject the new position
                self.neighbor_lists = initial_neighbor_lists
                self.atoms_positions = initial_positions
                self.number_atoms = initial_number_atoms
                self.calculate_cross_coefficients()

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

