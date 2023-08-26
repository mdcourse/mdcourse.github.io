Molecular dynamics simulations
==============================

The typical structure of a molecular dynamics simulation is the following:

.. code-block:: bash

    - Initialize the system
    - Loop for a desired number of steps
        * Measure the force on each atoms
        * Integrate the Newton's equation of motion
        * Increment the time
