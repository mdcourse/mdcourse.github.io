Minimize the energy
===================

.. container:: justify

    Now that the code for placing the atoms within the box has been written,
    let us perform an energy minimization to ensure there is no overlapping
    between the atoms.

Prepare the minimization
------------------------

.. container:: justify

    Let us start by importing NumPy and the copy libraries. Add the following
    to the beginning of the *MinimizeEnergy.py* file:

.. code-block:: python

    import numpy as np
    import copy
    from Outputs import Outputs
    (...)

.. container:: justify

    Then, let us fill the *__init__()* method:

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

.. container:: justify

    An important parameter is *maximum_steps*, which sets the maximum number of
    steps that the energy minimization will last. A *cut_off* with a default value
    of 9 Ångströms is also defined. The *neighbor* parameter sets the period
    between two recalculations of the neighbor lists, and the *displacement* with
    a default value of 0.01 Ångström sets the initial value for atom displacement.

Nondimensionalize units
-----------------------

.. container:: justify

    Two parameters from the *MinimizeEnergy* class must be nondimensionalized,
    namely *cut_off* and *displacement*. Add the following method to the
    *MinimizeEnergy* class:

.. code-block:: python

        def nondimensionalize_units_2(self):
            """Use LJ prefactors to convert units into non-dimensional."""
            self.cut_off = self.cut_off/self.reference_distance
            self.displacement = self.displacement/self.reference_distance

.. container:: justify

    Let us call the *nondimensionalize_units_2()* method from the *__init__()*
    method:

.. code-block:: python

    (...)
    self.maximum_steps = maximum_steps
    super().__init__(*args, **kwargs)
    self.nondimensionalize_units_2()

Energy minimizer
----------------

.. container:: justify

    Here the steepest descent method is used. First, the initial energy of the
    system and maximum force are measured. Then, a set of new positions for the atoms
    is tested. The energy of the new configuration is compared to the initial
    energy, and the new position is either accepted or rejected based on
    energy criteria.  

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
            self.update_log_minimize(Epot, max_forces)
            self.update_dump_file(filename="dump.min.lammpstrj")

Build neighbor lists
--------------------

.. container:: justify

    To save time, it is common in molecular simulation to detect which atoms are
    neighbors. This way, only interactions between neighbors are recalculated.
    This is the purpose of the *update_neighbor_lists()* method that must be
    added to the *Utilities* class:

.. code-block:: python

    def update_neighbor_lists(self):
        """Update the neighbor lists."""
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

.. container:: justify

    At the start of the *Utilities.py*, the NumPy library must be imported as well
    as the *distances* function from MDAnalysis using:

.. code-block:: python

    import numpy as np
    from MDAnalysis.analysis import distances

    from Potentials import LJ_potential
    (...)

.. container:: justify

    The *update_neighbor_lists()* method generates neighbor lists that are stored
    as Python list named *neighbor_lists*.

Compute_potential
-----------------

.. container:: justify

    Computing the potential is central to the energy minimizer.
    Add the following method called *compute_potential()*  to the *Utilities*
    class. 

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

Final code
----------

.. label:: start_Utilities_class

.. code-block:: python

    import numpy as np

    from MDAnalysis.analysis import distances
    from Potentials import LJ_potential

    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)
            
        def update_neighbor_lists(self):
            """Update the neighbor lists."""
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

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    import numpy as np
    import copy
    from Outputs import Outputs


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
            self.nondimensionalize_units_2()

        def nondimensionalize_units_2(self):
            """Use LJ prefactors to convert units into non-dimensional."""
            self.cut_off = self.cut_off/self.reference_distance
            self.displacement = self.displacement/self.reference_distance

        def run(self):
            """Perform energy minimmization using the steepest descent method."""
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
                self.update_log_minimize(Epot, max_forces)
                self.update_dump_file(filename="dump.min.lammpstrj")

.. label:: end_MinimizeEnergy_class


Test the code
-------------

.. container:: justify

    Let us test the *MinimizeEnergy* class to make sure that it does what
    is expected.

.. label:: start_test_MinimizeEnergy_class

.. code-block:: python

    from MinimizeEnergy import MinimizeEnergy

    self = MinimizeEnergy(maximum_steps=100,
        number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 1], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    self.run()

.. label:: end_test_MinimizeEnergy_class