Code structure
==============

The code is divided into fives classes:

- 1. InitializeSimulation
- 2. Outputs
- 3. Utilities
- 4. MolecularDynamics
- 5. MonteCarlo

Classes 1, 2, and 3 are being inherited by classes 4 and 5. The purpose of
class 1 is to set up the system, place atoms in the box. Class 2 deals with outputs
and logging of information during the simulations. Class 3 contains some 
functionalities such as pressure or temperature measurements. Finally classes 4 
and 5 are the core Molecular Dynamics and Monte Carlo codes, respectively.

Create a blank python file, call it *ms_code.py*, and copy the following
classes into it:

.. code-block:: python

    class InitializeSimulation:
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs) 

            print("Initialize Simulation")


    class Outputs:
        def __init__(self,
                     *args,
                     **kwargs):
            super().__init__(*args, **kwargs)

            print("Outputs")


    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

            print("Utilities")


    class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)

            print("Start molecular dynamics simulation")


    class MonteCarlo(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs)

            print("Start Monte Carlo simulation")

The *args* and *kwargs* arguments ensure that arguments of classes
*InitializeSimulation*, *Outputs*, *Utilities* are inherited by
the classes *MolecularDynamics* and *MonteCarlo*.

Tests
-----

Although the code is currently mostly empty, one can test that classes
are being inherited as expected.

.. code-block:: python

    from ms_code import Utilities, MolecularDynamics, MonteCarlo

    x = Utilities()
    print()
    y = MolecularDynamics()
    print()
    z = MonteCarlo()