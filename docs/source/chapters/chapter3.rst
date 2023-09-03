Monte carlo simulations
=======================

The typical structure of a Monte Carlo simulation is the following:

.. code-block:: bw

    - Initialize the system
    - Loop for a desired number of steps
        * Measure the potential energy of the system
        * Make a random attempt (e.g. move a particle, insert a particle, rotate a molecule, ...)
        * Measure the potential energy of the system again
        * Decide wether to keep the attempted move based on Metropolis criteria

Measure the potential energy
----------------------------

One of the most important step of Monte Carlo simulation is to evaluate the potential energy
of the system. Potential energy needs to be measured before and after any random Monte Carlo
attempt.

Assuming that the only potential of interaction in the system is the Lennard-Jones potential, 
one needs to calculate, for every pair of atoms :math:`i` and :math:`j`, with :math:`i < j`,

.. math::

    V_\text{LJ} (r) = \sum_{i < j} 4 \epsilon \left[ \left(\dfrac{\sigma}{r_{ij}}\right)^{12} - \left(\dfrac{\sigma}{r_{ij}}\right)^{6} \right],

where :math:`r_{ij}` is the distance between atoms :math:`i` and :math:`j`. In non-dimensionalized units,
it reads

.. math::

    V_\text{LJ} (r) = \sum_{i < j} 4 \left[ r_{ij}^{-12} - r_{ij}^{-6} \right],

The condition :math:`i < j` ensures that interaction are counted only once, and that interaction of an 
atom with itself is excluded.

Let us create a Python function, within the *Utilities* class, named *calculate_potential_energy*:

.. code-block:: python

    def calculate_potential_energy(self, atoms_positions):
        energy_potential = 0
        for position_i in atoms_positions:
            r = self.calculate_r(position_i, atoms_positions)
            energy_potential += np.sum(4*(1/np.power(r[r>0], 12)-1/np.power(r[r>0], 6)))
        energy_potential /= 2
        return energy_potential

The function *calculate_potential_energy* takes the *atoms_positions*, which is an array
containing all x, y, z coordinates of the atoms, and return a scalar, the total potential
energy of the system. Let us create another function called *calculate_r*, still within
the *Utilities* class, for the calculation of the inter atom distances. 

.. code-block:: python

    def calculate_r(self, position_i, positions_j):
        """Calculate the shortest distance between position_i and positions_j."""
        rij2 = np.zeros(self.number_atoms)
        box_size = np.diff(self.box_boundaries).reshape(3)
        rij = (np.remainder(position_i - positions_j + box_size/2., box_size) - box_size/2.).T
        for dim in np.arange(self.dimensions):
            rij2 += np.power(rij[dim, :], 2)
        return np.sqrt(rij2)

In principle, inter atoms positions could simply be calculated as *position_i - positions_j*.
However, due to the periodic boundary condition, the simple direct evaluation *position_i - positions_j*
would not always return the shortest distance between the atom :math:`i` and all the others atoms :math:`j`,
which explain the more complex relation used. 

Monte Carlo step
----------------

.. code-block:: python

    def monte_carlo_displacement(self):
        beta =  1/self.desired_temperature
        Epot = self.calculate_potential_energy(self.atoms_positions)
        trial_atoms_positions = copy.deepcopy(self.atoms_positions)
        atom_id = np.random.randint(self.number_atoms)
        trial_atoms_positions[atom_id] += (np.random.random(3)-0.5)*self.displace_mc
        trial_Epot = self.calculate_potential_energy(trial_atoms_positions)
        acceptation_probability = np.min([1, np.exp(-beta*(trial_Epot-Epot))])
        if np.random.random() <= acceptation_probability:
            self.atoms_positions = trial_atoms_positions
            self.Epot = trial_Epot
        else:
            self.Epot = Epot 

Calculate Lambda
----------------

.. code-block:: python

    def calculate_Lambda(self, mass):
        """Estimate de Broglie wavelength in LJ units."""
        m_kg = mass/cst.Avogadro*cst.milli*self.reference_mass
        kB_kCal_mol_K = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo
        T_K = self.desired_temperature*self.reference_energy/kB_kCal_mol_K
        Lambda = cst.h/np.sqrt(2*np.pi*cst.Boltzmann*m_kg*T_K)/cst.angstrom
        return self.nondimensionalise_units(Lambda, "distance")

Main loop
----------

.. code-block:: python

    class MonteCarlo(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
            maximum_steps,
            cut_off = 10,
            displace_mc=None,
            *args,
            **kwargs,
            ):

            self.maximum_steps = maximum_steps
            self.cut_off = cut_off
            self.displace_mc = displace_mc
            super().__init__(*args, **kwargs)

            self.cut_off = self.nondimensionalise_units(self.cut_off, "distance")
            self.displace_mc = self.nondimensionalise_units(self.displace_mc, "distance")

        def run(self):
            """Perform the loop over time."""
            
            for self.step in range(0, self.maximum_steps+1):
                self.monte_carlo_displacement()
                self.wrap_in_box()
                self.update_log()
                self.update_dump(velocity=False)
            self.write_lammps_data(filename="final.data")

