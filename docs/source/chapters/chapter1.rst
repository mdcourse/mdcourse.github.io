Code structure
==============

Presentation
------------

.. container:: justify

    The Python code will be divided into five classes:

    - 1. InitializeSimulation
    - 2. Outputs
    - 3. Utilities
    - 4. MolecularDynamics
    - 5. MonteCarlo

    Here, classes 1, 2, and 3 are being inherited by classes 4 and 5. 
    The purpose of class 1 is to set up the system and place atoms in the
    box. Class 2 deals with outputs and logging of information during
    the simulations. Class 3 contains some functionalities such as pressure
    or temperature measurements. Finally classes 4 
    and 5 are the core Molecular Dynamics and Monte Carlo codes, respectively.

Start coding
-------------

.. container:: justify

    Inside a dedicated folder, create 5 blank Python scripts named:

    - InitializeSimulation.py
    - Outputs.py
    - Utilities.py
    - MolecularDynamics.py
    - MonteCarlo.py

.. container:: justify

    For each script, let us write the basis of each script.
    Let us start with InitializeSimulation.py: 

.. code-block:: python

    class InitializeSimulation:
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs) 

            print("Initialize Simulation")

.. container:: justify

    The Outputs.py and Utilities.py are similar:

.. code-block:: python

    class Outputs:
        def __init__(self,
                     *args,
                     **kwargs):
            super().__init__(*args, **kwargs)

            print("Outputs")

.. code-block:: python

    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

            print("Utilities")

.. container:: justify

    The MolecularDynamics.py and MonteCarlo.py are inheriting
    the 3 previously defined classes:

.. code-block:: python

    class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)

            print("Start molecular dynamics simulation")

.. code-block:: python

    class MonteCarlo(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                     *args,
                     **kwargs,
                     ):
            super().__init__(*args, **kwargs)

            print("Start Monte Carlo simulation")

.. container:: justify

    The *args* and *kwargs* arguments ensure that arguments of classes
    *InitializeSimulation*, *Outputs*, *Utilities* are inherited by
    the classes *MolecularDynamics* and *MonteCarlo*.

Tests
-----

.. container:: justify

    Although the code is currently mostly empty, one can test that classes
    are being inherited as expected.

.. code-block:: python

    from ms_code import Utilities, MolecularDynamics, MonteCarlo

    x = Utilities()
    print()
    y = MolecularDynamics()
    print()
    z = MonteCarlo()