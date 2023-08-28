Structure of the code
=====================

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

The overall structure of the code is the following:

.. code-block:: python

    class InitializeSimulation:
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs) 


    class Outputs:
        def __init__(self,
                     *args,
                     **kwargs):
            super().__init__(*args, **kwargs)


    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)


    class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)


    class MonteCarlo(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs)
