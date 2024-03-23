Code structure
==============

.. container:: justify

    The main files and classes that constitute the basis of the code
    are defined.

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
                     **kwargs):
            super().__init__(*args, **kwargs) 

.. container:: justify

    The Utilities.py is similar:

.. code-block:: python

    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. container:: justify

    For the Outputs.py class, let us anticipate that the outputs
    from the code will be saved in a folder, which by default
    is named *results/*. If the folder does not exist, it will be
    created using *os.mkdir*:

.. code-block:: python

    import os

    class Outputs:
        def __init__(self,
                    data_folder = "results/",
                    *args,
                    **kwargs):
            self.data_folder = data_folder
            super().__init__(*args, **kwargs)

            if os.path.exists(self.data_folder) is False:
                os.mkdir(self.data_folder)

.. container:: justify

    The MolecularDynamics class is inheriting
    the 3 previously defined classes:

.. code-block:: python

    from InitializeSimulation import InitializeSimulation
    from Utilities import Utilities
    from Outputs import Outputs

    class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            pass

.. container:: justify

    The *run* method will be filled later. Let us do the same for the
    MonteCarlo class:

.. code-block:: python

    from InitializeSimulation import InitializeSimulation
    from Utilities import Utilities
    from Outputs import Outputs

    class MonteCarlo(InitializeSimulation, Utilities, Outputs):
        def __init__(self,
                     *args,
                     **kwargs):
            super().__init__(*args, **kwargs)

        def run(self):
            pass

.. container:: justify

    The *args* and *kwargs* arguments ensure that arguments of classes
    *InitializeSimulation*, *Outputs*, *Utilities* are inherited by
    the classes *MolecularDynamics* and *MonteCarlo*.

Test the code
-------------

.. container:: justify

    We can create a simple test to ensure that the classes
    are being inherited as expected. Within the same folder,
    create a new Jupyter notebook called *test.ipynb*, and copy
    the following lines into it:

.. code-block:: python

    from InitializeSimulation import InitializeSimulation
    from Utilities import Utilities
    from Outputs import Outputs
    from MolecularDynamics import MolecularDynamics
    from MonteCarlo import MonteCarlo

    md = MolecularDynamics(data_folder = "md-output/")
    md.run()
    mc = MolecularDynamics(data_folder = "mc-output/")
    mc.run()

.. container:: justify

    If everything is working just fine, two folders named *md-output/*
    and *mc-output/* must have been created.