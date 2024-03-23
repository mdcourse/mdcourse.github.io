LJ parameters and units
=======================

.. container:: justify

    The normalization of the basic parameters is defined here. 

Unit systems
------------

.. container:: justify

    Two separate unit systems are used. The first unit system is called
    *real*, and followed the convention from LAMMPS. All inputs will be 
    provided to the code in *real* units, for which:

    - mass in grams/mole
    - distance in Angstroms
    - time in femtoseconds
    - energy in kcal/mol
    - velocity in Angstroms/femtosecond
    - force in (kcal/mol)/Angstrom
    - temperature in Kelvin
    - pressure in atmospheres
    - density in g/cm^dim

.. container:: justify

    The *real* unit system is conventional in molecular simulations. However,
    it is not practical to perform the calculations with such complex units,
    as it would involve complex prefactors. Instead, the LJ (for Lennard-Jones)
    unit system will be used. With the LJ unit systems, all quantities are
    unitless. All masses, distances, and energies are specified as multiples 
    of :math:`m`, :math:`\sigma`, and :math:`\epsilon`, which are the mass and LJ
    parameters of the atoms. Other quantities are specified from these 3 quantities:

    - mass in :math:`m`
    - distance in :math:`\sigma`
    - time in :math:`\sqrt{m \sigma^2 / \epsilon}`
    - energy in :math:`\epsilon`
    - velocity in :math:`\sqrt{m \sigma^4 / \epsilon}`
    - force in :math:`\epsilon/\sigma`
    - temperature in :math:`\epsilon/k_\text{B}`
    - pressure in :math:`\epsilon/\sigma^3`
    - density in :math:`m/\sigma^3`

    where :math:`k_\text{B}` is the Boltzmann constant. 