Molecular dynamics simulations
==============================

The typical structure of a molecular dynamics simulation is the following:

.. code-block:: bw

    - Initialize the system
    - Loop for a desired number of steps
        * Measure the force on each atoms
        * Integrate the Newton s equations of motion
        * Increment the time
