Pressure measurement
====================

Extract the pressure
--------------------

Let us also measure the pressure. 

.. label:: start_Utilities_class

.. code-block:: python

    def calculate_pressure(self):
        """Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smith 2002)"""
        Ndof = self.dimensions*self.total_number_atoms-self.dimensions    
        volume = np.prod(self.box_size)
        try:
            self.calculate_temperature() # this is for later on, when velocities are computed
        except:
            self.temperature = self.desired_temperature # for MC, simply use the desired temperature
        p_ideal = (Ndof/self.dimensions)*self.temperature/volume
        p_non_ideal = 1/(volume*self.dimensions)*np.sum(self.compute_potential(output="force-matrix")*self.evaluate_rij_matrix())
        self.pressure = (p_ideal+p_non_ideal)

.. label:: end_Utilities_class

One needs to evaluate the vector between all the particles.

.. label:: start_Utilities_class

.. code-block:: python

    def evaluate_rij_matrix(self):
        """Matrix of vectors between all particles."""
        rij_matrix = np.zeros((self.total_number_atoms,self.total_number_atoms,3))
        for Ni in range(self.total_number_atoms-1):
            position_i = self.atoms_positions[Ni]
            positions_j = self.atoms_positions
            rij_xyz = (np.remainder(position_i - positions_j
                                    + self.box_size[:3]/2., self.box_size[:3])
                                    - self.box_size[:3]/2.)
            rij_matrix[Ni] = rij_xyz
        return rij_matrix

.. label:: end_Utilities_class

Test the code
-------------

One can use a similar test as previously. Let us use a displace distance of
0.5 Angstrom, and make 1000 steps. The pressure will be written in *pressure.dat*.

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