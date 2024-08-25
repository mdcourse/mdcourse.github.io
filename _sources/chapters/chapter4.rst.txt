Minimize the energy
===================

Now that the code for placing the atoms within the box has been written,
let us perform an energy minimization to ensure that there is no overlapping
between the atoms, before starting the actual simulation.

Prepare the minimization
------------------------

Let us start by importing NumPy and the copy libraries. Add the following
to the beginning of the *MinimizeEnergy.py* file:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    import numpy as np
    import copy

.. label:: end_MinimizeEnergy_class

Then, let us fill the *__init__()* method:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    class MinimizeEnergy(Outputs):
        def __init__(self,
                    maximum_steps,
                    cut_off=9,
                    neighbor=1,
                    displacement=0.01,
                    *args,
                    **kwargs):
            self.neighbor = neighbor
            self.cut_off = cut_off
            self.displacement = displacement
            self.maximum_steps = maximum_steps
            super().__init__(*args, **kwargs)

.. label:: end_MinimizeEnergy_class

An important parameter is *maximum_steps*, which sets the maximum number of
steps that the energy minimization will last. A *cut_off* with a default value
of 9 Ångströms is also defined. The *neighbor* parameter sets the period
between two recalculations of the neighbor lists, and the *displacement* with
a default value of 0.01 Ångström sets the initial value for atom displacement.

Nondimensionalize units
-----------------------

Two parameters from the *MinimizeEnergy* class must be nondimensionalized,
namely *cut_off* and *displacement*. Add the following method to the
*MinimizeEnergy* class:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def nondimensionalize_units_2(self):
        """Use LJ prefactors to convert units into non-dimensional."""
        self.cut_off = self.cut_off/self.reference_distance
        self.displacement = self.displacement/self.reference_distance

.. label:: end_MinimizeEnergy_class

Finally, let us call the *nondimensionalize_units_2()* method from the *__init__()*
method:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def __init__(self,
        (...)
        super().__init__(*args, **kwargs)
        self.nondimensionalize_units_2()

.. label:: end_MinimizeEnergy_class

Energy minimizer
----------------

Here, the steepest descent method is used. First, the initial energy of the
system and the maximum force are measured. Then, a set of new positions for the atoms
is tested. The energy of the new configuration is compared to the initial
energy, and the new position is either accepted or rejected based on energy criteria.  

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def run(self):
        for self.step in range(0, self.maximum_steps+1):
            # Measure the initial energy and max force
            self.update_neighbor_lists()
            init_Epot = self.compute_potential(output="potential")
            initial_positions = copy.deepcopy(self.atoms_positions)
            forces = self.compute_potential(output="force-vector")
            max_forces = np.max(np.abs(forces))
            # Test a new sets of positions
            self.atoms_positions = self.atoms_positions \
                + forces/max_forces*self.displacement
            trial_Epot = self.compute_potential(output="potential")
            # Keep the more favorable energy
            if trial_Epot < init_Epot:  # accept new position
                Epot = trial_Epot
                self.wrap_in_box()
                self.displacement *= 1.2
            else:  # reject new position
                Epot = init_Epot
                self.atoms_positions = initial_positions
                self.displacement *= 0.2

.. label:: end_MinimizeEnergy_class

The displacement, which has an initial value of 0.01, is adjusted through energy
minimization. When the trial is successful, its value is multiplied by 1.2. When
the trial is rejected, its value is multiplied by 0.2.

Build neighbor lists
--------------------

To save time, it is common in molecular simulations to detect which atoms are
neighbors. This way, only interactions between neighbors are recalculated.
This is the purpose of the *update_neighbor_lists()* method that must be
added to the *Utilities* class:

.. label:: start_Utilities_class

.. code-block:: python

    def update_neighbor_lists(self):
        if (self.step % self.neighbor == 0):

            matrix = distances.contact_matrix(self.atoms_positions,
                cutoff=self.cut_off, #+2,
                returntype="numpy",
                box=self.box_size)

            cpt = 0
            neighbor_lists = []
            for array in matrix[:-1]:
                list = np.where(array)[0].tolist()
                list = [ele for ele in list if ele > cpt]
                cpt += 1
                neighbor_lists.append(list)

            self.neighbor_lists = neighbor_lists

.. label:: end_Utilities_class

The *update_neighbor_lists()* method generates neighbor lists that are stored
as Python list named *neighbor_lists*. The following library
should also be imported within the *Utilities.py* file:

.. label:: start_Utilities_class

.. code-block:: python

    import numpy as np

    from MDAnalysis.analysis import distances

.. label:: end_Utilities_class

Compute_potential
-----------------

Computing the potential energy of the system is central to the energy minimizer,
as the value of the potential is used to decide if the trial is accepted or
rejected. Add the following method called *compute_potential()*  to the *Utilities*
class.

.. label:: start_Utilities_class

.. code-block:: python

    def compute_potential(self, output):
        if output == "force-vector":
            forces = np.zeros((self.total_number_atoms,3))
        elif output == "force-matrix":
            forces = np.zeros((self.total_number_atoms,self.total_number_atoms,3))
        energy_potential = 0
        for Ni in np.arange(self.total_number_atoms-1):
            # Read information about atom i
            position_i = self.atoms_positions[Ni]
            sigma_i = self.atoms_sigma[Ni]
            epsilon_i = self.atoms_epsilon[Ni]
            neighbor_of_i = self.neighbor_lists[Ni]
            # Read information about neighbors j
            positions_j = self.atoms_positions[neighbor_of_i]
            sigma_j = self.atoms_sigma[neighbor_of_i]
            epsilon_j = self.atoms_epsilon[neighbor_of_i]
            # Measure distances and other cross parameters
            rij_xyz = (np.remainder(position_i - positions_j
                                    + self.box_size[:3]/2., self.box_size[:3])
                                    - self.box_size[:3]/2.)
            rij = np.linalg.norm(rij_xyz, axis=1)
            sigma_ij = (sigma_i+sigma_j)/2
            epsilon_ij = (epsilon_i+epsilon_j)/2
            # Measure potential
            if output == "potential":
                energy_potential += np.sum(LJ_potential(epsilon_ij, sigma_ij, rij))
            else:
                derivative_potential = LJ_potential(epsilon_ij, sigma_ij, rij, derivative = True)
                if output == "force-vector":
                    forces[Ni] += np.sum((derivative_potential*rij_xyz.T/rij).T, axis=0)
                    forces[neighbor_of_i] -= (derivative_potential*rij_xyz.T/rij).T 
                elif output == "force-matrix":
                    forces[Ni][neighbor_of_i] += (derivative_potential*rij_xyz.T/rij).T
        if output=="potential":
            return energy_potential
        elif (output == "force-vector") | (output == "force-matrix"):
            return forces

.. label:: end_Utilities_class

Wrap in box
-----------

Every time atoms are being displaced, one has to ensure that they remain in
the box. This is done by the *wrap_in_box()* method that must be placed
within the *Utilities* class:

.. label:: start_Utilities_class

.. code-block:: python

    def wrap_in_box(self):
        for dim in np.arange(self.dimensions):
            out_ids = self.atoms_positions[:, dim] \
                > self.box_boundaries[dim][1]
            self.atoms_positions[:, dim][out_ids] \
                -= np.diff(self.box_boundaries[dim])[0]
            out_ids = self.atoms_positions[:, dim] \
                < self.box_boundaries[dim][0]
            self.atoms_positions[:, dim][out_ids] \
                += np.diff(self.box_boundaries[dim])[0]

.. label:: end_Utilities_class

Test the code
-------------

Let us test the *MinimizeEnergy* class to make sure that it does what
is expected, i.e. that it leads to a potential energy that is small, and
typically negative.

.. label:: start_test_MinimizeEnergy_class

.. code-block:: python

    from MinimizeEnergy import MinimizeEnergy

    min = MinimizeEnergy(maximum_steps=100,
        number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 1], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    min.run()

    assert min.compute_potential(output="potential") < 0

.. label:: end_test_MinimizeEnergy_class
