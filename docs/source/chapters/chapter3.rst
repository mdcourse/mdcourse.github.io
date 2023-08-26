Monte carlo simulations
=======================

The typical structure of a Monte Carlo simulation is the following:

.. code-block:: bash

    - Initialize the system
    - Loop for a desired number of steps
        * Measure the potential energy of the system
        * Make a random attempt (e.g. move a particle, insert a particle, rotate a molecule, ...)
        * Measure the potential energy of the system again
        * Decide wether to keep the attempted move based on Metropolis criteria
