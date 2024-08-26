.. _chapter7-label:

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

Let us test the outputed pressure. An interesting test is to contront the output from
our code with some data from the litterature. Let us used the same parameters
as in Ref. :cite:`woodMonteCarloEquation1957`, where Monte Carlo simulations
are used to simulate argon bulk phase. All we have to do is to apply our current
code using their parameters, i.e. :math:`\sigma = 3.405~\text{Å}`, :math:`\epsilon = 0.238~\text{kcal/mol}`,
and :math:`T = 55~^\circ\text{C}`. More details are given in the first illustration, :ref:`project1-label`.

On the side note, a relatively small cut-off as well as a small number of atoms were
chosen to make the calculation faster. 

.. label:: start_test_MonteCarloPressure_class

.. code-block:: python

    import numpy as np
    from MonteCarlo import MonteCarlo

    from scipy import constants as cst
    from pint import UnitRegistry
    ureg = UnitRegistry()

    # Constants
    kB = cst.Boltzmann*ureg.J/ureg.kelvin # boltzman constant
    Na = cst.Avogadro/ureg.mole # avogadro
    R = kB*Na # gas constant

    # Parameters taken from Wood1957
    tau = 2 # ratio between volume / reduced volume
    epsilon = (119.76*ureg.kelvin*kB*Na).to(ureg.kcal/ureg.mol) # kcal/mol
    r_star = 3.822*ureg.angstrom # angstrom
    sigma = r_star / 2**(1/6) # angstrom
    N_atom = 50 # no units
    m_argon = 39.948*ureg.gram/ureg.mol # g/mol
    T =  328.15 * ureg.degK # 328 K or 55°C
    volume_star = r_star**3 * Na * 2**(-0.5)
    cut_off = sigma*2.5 # angstrom
    displace_mc = sigma/5 # angstrom
    volume = N_atom*volume_star*tau/Na # angstrom**3
    box_size = volume**(1/3) # angstrom

    mc = MonteCarlo(maximum_steps=15000,
        dumping_period=1000,
        thermo_period=1000,
        neighbor=50,
        number_atoms=[N_atom],
        epsilon=[epsilon.magnitude],
        sigma=[sigma.magnitude],
        atom_mass=[m_argon.magnitude],
        box_dimensions=[box_size.magnitude, box_size.magnitude, box_size.magnitude],
        displace_mc = displace_mc.magnitude,
        desired_temperature = T.magnitude,
        cut_off = cut_off.magnitude,
        )
    mc.run()

    # Import the data and calculate p V / R T
    output = np.mean(np.loadtxt("Outputs/pressure.dat")[:,1][10:])
    pressure = (output*ureg.atm).to(ureg.pascal)
    volume = (volume_star * tau / Na).to(ureg.meter**3)
    pV_over_RT = np.round((pressure * volume / (R * T) * Na).magnitude,2)
    print("p v / R T =", pV_over_RT, " --- (The expected value from Wood1957 is 1.5)")

.. label:: end_test_MonteCarloPressure_class

Which should return a value for :math:`p V / R T` that is close to the expected value
of 1.5 by Wood and Parker for :math:`\tau = V/V^* = 2` (see Fig. 4 in Ref. :cite:`woodMonteCarloEquation1957`):

.. code-block:: bw

    (...)
    p v / R T = 1.56  --- (The expected value from Wood1957 is 1.5)

The exact value will varie from one simulation to the other due to noise.