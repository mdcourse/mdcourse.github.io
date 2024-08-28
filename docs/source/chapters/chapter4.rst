.. _chapter4-label:

Minimize the energy
===================

Now that the code for placing the atoms within the box has been written,
let us proceed to write the code for performing energy minimization of the
system. This step helps ensure that there is no unreasonable overlapping
between the atoms.

The steepest descent method is used for energy minimization, following these steps:

- 1) Start with an initial configuration and measure the initial potential energy,
     :math:`E_\text{pot}^\text{initial}`.
- 2) Calculate the gradient of the potential energy with respect to atomic positions
     to determine the direction of the steepest ascent in energy. The magnitude
     of this gradient is used to compute the maximum force acting on the atoms.
     This maximum force indicates the direction in which the energy decreases most
     rapidly.
- 3) Move the atoms in the opposite direction of the maximum
     force to minimize the potential energy by a displacement step.
     The size of the step is adjusted iteratively based on the reduction in energy.
- 4) Compute the new potential energy after the displacement, :math:`E_\text{pot}^\text{trial}`.
- 5) Evaluate the change in energy: :math:`\Delta E = E_\text{pot}^\text{trial} - E_\text{pot}^\text{initial}`.
  
  - If :math:`\Delta E < 0`, the new configuration is accepted as it results in
    lower energy, and the step size is increased.
  - If :math:`\Delta E \geq 0`, the new configuration is rejected, and the step
    size is decreased.

The process is repeated until the maximum number of steps is reached.
The goal is to iteratively reduce the potential energy and reach a stable,
minimized energy state.

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

    class MinimizeEnergy(Measurements):
        def __init__(self,
                    maximum_steps,
                    cut_off=9,
                    neighbor=1,
                    displacement=0.01,
                    thermo_outputs="MaxF",
                    data_folder=None,
                    *args,
                    **kwargs):
            self.neighbor = neighbor
            self.cut_off = cut_off
            self.displacement = displacement
            self.maximum_steps = maximum_steps
            self.thermo_outputs = thermo_outputs
            self.data_folder = data_folder
            if self.data_folder is not None:
                if os.path.exists(self.data_folder) is False:
                    os.mkdir(self.data_folder)
            super().__init__(*args, **kwargs)

.. label:: end_MinimizeEnergy_class

An important parameter is *maximum_steps*, which sets the maximum number
of steps for the energy minimization process. A *cut_off* value with a
default of 9 Ångströms is also defined. The *neighbor* parameter determines
the interval between recalculations of the neighbor lists, and the *displacement*
parameter, with a default value of 0.01 Ångström, sets the initial atom
displacement value.

The *thermo_outputs* and *data_folder* parameters are used for printing data
to files. These two parameters will be useful in the next chapter, :ref:`chapter5-label`.

Nondimensionalize units
-----------------------

As was done previously, some parameters from the *MinimizeEnergy* class
must be non-dimensionalized: *cut_off* and *displacement*. Add the following
method to the *MinimizeEnergy* class:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def nondimensionalize_units_2(self):
        """Use LJ prefactors to convert units into non-dimensional."""
        self.cut_off = self.cut_off/self.reference_distance
        self.displacement = self.displacement/self.reference_distance

.. label:: end_MinimizeEnergy_class

Let us call the *nondimensionalize_units_2()* method from the *__init__()*
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

Let us implement the energy minimized described at the top of this page. Add the
following *run()* method to the *MinimizeEnergy* class:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def run(self):
        for self.step in range(0, self.maximum_steps+1): # *step* loops for 0 to *maximum_steps*+1
            # First, meevaluate the initial energy and max force
            self.update_neighbor_lists() # Rebuild neighbor list, if necessary
            self.update_cross_coefficients() # Recalculate the cross coefficients, if necessary
            try: # try using the last saved Epot and MaxF, if it exists
                init_Epot = self.Epot
                init_MaxF = self.MaxF
            except: # If Epot/MaxF do not exists yet, calculate them both
                init_Epot = self.compute_potential(output="potential")
                forces = self.compute_potential(output="force-vector")
                init_MaxF = np.max(np.abs(forces))
            # Save the current atom positions
            init_positions = copy.deepcopy(self.atoms_positions)
            # Move the atoms in the opposite direction of the maximum force
            self.atoms_positions = self.atoms_positions \
                + forces/init_MaxF*self.displacement
            # Recalculate the energy
            trial_Epot = self.compute_potential(output="potential")
            # Keep the more favorable energy
            if trial_Epot < init_Epot: # accept new position
                self.Epot = trial_Epot
                # calculate the new max force and save it
                forces = self.compute_potential(output="force-vector")
                self.MaxF = np.max(np.abs(forces))
                self.wrap_in_box()  # Wrap atoms in the box, if necessary
                self.displacement *= 1.2 # Multiply the displacement by a factor 1.2
            else: # reject new position
                self.Epot = init_Epot # Revert to old energy
                self.atoms_positions = init_positions # Revert to old positions
                self.displacement *= 0.2 # Multiply the displacement by a factor 0.2

.. label:: end_MinimizeEnergy_class

The displacement, which has an initial value of 0.01, is adjusted through energy
minimization. When the trial is successful, its value is multiplied by 1.2. When
the trial is rejected, its value is multiplied by 0.2.

Build neighbor lists
--------------------

In molecular simulations, it is common practice to identify neighboring atoms
to save computational time. By focusing only on interactions between
neighboring atoms, the simulation becomes more efficient. Add the following
*update_neighbor_lists()* method to the *Utilities*class:

.. label:: start_Utilities_class

.. code-block:: python

    def update_neighbor_lists(self):
        if (self.step % self.neighbor == 0):
            matrix = distances.contact_matrix(self.atoms_positions,
                cutoff=self.cut_off, #+2,
                returntype="numpy",
                box=self.box_size)
            neighbor_lists = []
            for cpt, array in enumerate(matrix[:-1]):
                list = np.where(array)[0].tolist()
                list = [ele for ele in list if ele > cpt]
                neighbor_lists.append(list)
            self.neighbor_lists = neighbor_lists

.. label:: end_Utilities_class

The *update_neighbor_lists()* method generates neighbor lists for each
atom, ensuring that only relevant interactions are considered in the
calculations. These lists will be recalculated at intervals specified by
the *neighbor* input parameter.

Update cross coefficients
-------------------------

At the same time as the neighbor lists are getting build up, let us also
pre-calculate the cross coefficients. This will make the force calculation
more practical (see below).

.. label:: start_Utilities_class

.. code-block:: python

    def update_cross_coefficients(self):
        if (self.step % self.neighbor == 0):
            # Precalculte LJ cross-coefficients
            sigma_ij_list = []
            epsilon_ij_list = []
            for Ni in np.arange(self.total_number_atoms-1): # tofix error for GCMC
                # Read information about atom i
                sigma_i = self.atoms_sigma[Ni]
                epsilon_i = self.atoms_epsilon[Ni]
                neighbor_of_i = self.neighbor_lists[Ni]
                # Read information about neighbors j
                sigma_j = self.atoms_sigma[neighbor_of_i]
                epsilon_j = self.atoms_epsilon[neighbor_of_i]
                # Calculare cross parameters
                sigma_ij_list.append((sigma_i+sigma_j)/2)
                epsilon_ij_list.append((epsilon_i+epsilon_j)/2)
            self.sigma_ij_list = sigma_ij_list
            self.epsilon_ij_list = epsilon_ij_list

.. label:: end_Utilities_class

Here, the values of the cross coefficients between atom of type 1 and 2,
:math:`\sigma_{12}` and :math:`\epsilon_{12}`, are assumed to follow the arithmetic mean:

.. math::

    \sigma_{12} = (\sigma_{11}+\sigma_{22})/2 \\
    \epsilon_{12} = (\epsilon_{11}+\epsilon_{22})/2

Finally, import the following library in the *Utilities.py* file:

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
        box_size = self.box_size[:3]
        half_box_size = self.box_size[:3]/2.0
        for Ni in np.arange(self.total_number_atoms-1):
            # Read information about atom i
            position_i = self.atoms_positions[Ni]
            neighbor_of_i = self.neighbor_lists[Ni]
            # Read information about neighbors j and cross coefficient
            positions_j = self.atoms_positions[neighbor_of_i]
            sigma_ij = self.sigma_ij_list[Ni]
            epsilon_ij = self.epsilon_ij_list[Ni]
            # Measure distances
            rij_xyz = (np.remainder(position_i - positions_j + half_box_size, box_size) - half_box_size)
            rij = np.linalg.norm(rij_xyz, axis=1)
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

Here, the method is a little bit complicated, because three types of outputs can
be requested by the user: *force-vector*, *force-matrix*, and *potential*. The last
one, *potential*, simply returns the value of the potential energy for the entire system.
If *force-vector* or *force-matrix* are selected instead, then the individual forces
between atoms are returned.

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

.. label:: start_test_4a_class

.. code-block:: python

    from MinimizeEnergy import MinimizeEnergy

    min = MinimizeEnergy(maximum_steps=100,
        number_atoms=[2, 3],
        epsilon=[0.2, 0.4], # kcal/mol
        sigma=[3, 4], # A
        atom_mass=[10, 20], # g/mol
        box_dimensions=[20, 20, 20], # A
        )
    min.run()

    Final_Epot = min.Epot
    Final_MaxF = min.MaxF
    assert Final_Epot < 0, f"Test failed: Final energy too large"
    assert Final_MaxF < 10, f"Test failed: Final max force too large"
    print("Test passed")

.. label:: end_test_4a_class

For such as low density in particle, we can reasonably expect the energy to be always
negative after 100 steps.