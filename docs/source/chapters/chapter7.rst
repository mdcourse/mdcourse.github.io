Pressure measurement
====================

In order to extract the equation of state in our simulation, we need to measure the
pressure of the system, :math:`p`. The pressure in a molecular simulation can be
calculated from the interactions between particles. The pressure can be measured as
the sum of the ideal contribution, :math:`p_\text{ideal} = N_\text{DOF} k_\text{B} T / V d`,
which comes from the ideal gas law, and a Virial term which accounts for the
pressure contribution from the forces between particles,
:math:`p_\text{non_ideal} = \left< \sum_i r_i \cdot F_i \right> / V d`. The final
expression reads:

.. math:: 

    p = \dfrac{1}{V d} \left[ N_\text{DOF} k_\text{B} T +  \left< \sum_i r_i \cdot F_i \right> \right]

:math:`N_\text{DOF}` is the number of degrees-of-freedom, which can be calculated
from the number of particles, :math:`N`, and the dimension of the system, :math:`d`, as 
:math:`N_\text{DOF} = d N - d` :cite:`frenkel2023understanding`.

The calculation of :math:`p_\text{ideal}` is straighforward. For Monte Carlo simulation,
as atoms do not have temperature, the *imposed* temperature will be used instead.
The calculation of :math:`p_\text{non_ideal}` requires the measurement of all the
force and distance between the atoms. The calculation of the forces was already
implemented in a previous chapter, but a new function that returns all the
vector direction between atoms pairs will have to be written here.

Implement the Virial equation
-----------------------------

Let us add the following method to the *Utilities* class.

.. label:: start_Utilities_class

.. code-block:: python

    def calculate_pressure(self):
        """Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smit,
        Understanding molecular simulation: from algorithms to applications, 2002)"""
        # Ideal contribution
        Ndof = self.dimensions*self.total_number_atoms-self.dimensions    
        volume = np.prod(self.box_size[:3])
        try:
            self.calculate_temperature() # this is for later on, when velocities are computed
            temperature = self.temperature
        except:
            temperature = self.desired_temperature # for MC, simply use the desired temperature
        p_ideal = Ndof*temperature/(volume*self.dimensions)
        # Non-ideal contribution
        distances_forces = np.sum(self.compute_potential(output="force-matrix")*self.evaluate_rij_matrix())
        p_nonideal = distances_forces/(volume*self.dimensions)
        # Final pressure
        self.pressure = p_ideal+p_nonideal

.. label:: end_Utilities_class

To evaluate all the vectors between all the particles, let us also add the
*evaluate_rij_matrix* method to the *Utilities* class:

.. label:: start_Utilities_class

.. code-block:: python

    def evaluate_rij_matrix(self):
        """Matrix of vectors between all particles."""
        box_size = self.box_size[:3]
        half_box_size = self.box_size[:3]/2.0
        rij_matrix = np.zeros((self.total_number_atoms,self.total_number_atoms,3))
        positions_j = self.atoms_positions
        for Ni in range(self.total_number_atoms-1):
            position_i = self.atoms_positions[Ni]
            rij_xyz = (np.remainder(position_i - positions_j + half_box_size, box_size) - half_box_size)
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