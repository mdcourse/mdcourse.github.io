Measure the potential energy
============================

.. container:: justify

    The positions of the atoms are used to calculate
    the potential energy of the system.

Calculate the inter-atomic distances
------------------------------------

.. container:: justify

    Add the following method to the *Utilities* class.

.. code-block:: python

    def calculate_r(self, position_i, positions_j):
        rij = (np.remainder(position_i - positions_j
                            + self.box_size/2., self.box_size) - self.box_size/2.)
        return np.linalg.norm(rij, axis=1)

.. container:: justify

    The NumPy remainder function is used to calculate the shortest distance between
    *position_i* and *positions_j*, which account for the periodic boundary conditions.

Calculate the potential energy
-------------------------------

.. container:: justify

    For all the atoms, the inter-atomic distances are calculated using
    the previously defined function *calculate_r*, and 
    the potential energy is calculated by solving the Lennard-Jones potential
    for all pairs of atoms.

.. code-block:: python

    def calculate_potential_energy(self, atoms_positions):
        """Calculate potential energy from Lennard-Jones potential."""
        energy_potential = 0
        for position_i, sigma_i, epsilon_i in zip(atoms_positions,
                                                  self.atoms_sigma,
                                                  self.atoms_epsilon):
            r = self.calculate_r(position_i, atoms_positions)
            sigma_j = self.atoms_sigma
            epsilon_j = self.atoms_epsilon
            sigma_ij = np.array((sigma_i+sigma_j)/2)
            epsilon_ij = np.array((epsilon_i+epsilon_j)/2)
            energy_potential_i = np.sum(4*epsilon_ij[r>0]*(np.power(sigma_ij[r>0]/r[r>0], 12)-np.power(sigma_ij[r>0]/r[r>0], 6)))
            energy_potential += energy_potential_i
        return energy_potential/2
