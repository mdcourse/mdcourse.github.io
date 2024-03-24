Code structure
==============

.. container:: justify

    Here, the five files that constitute the basis of the code
    are created. These files will be progressively filled. 

Presentation
------------

.. container:: justify

    The Python code performing the molecular simulations will be
    divided into five classes:

    - *InitializeSimulation*
    - *Outputs*
    - *Utilities*
    - *MolecularDynamics*
    - *MonteCarlo*

    The three first classes named *InitializeSimulation*, *Outputs*, and *Utilities*
    are inherited by the two former classes *MolecularDynamics* and *MonteCarlo*. 

.. container:: justify

    The purpose of the class *InitializeSimulation* is to set up the
    system and place atoms in the box. Class *Outputs* deals with
    data files, dumps, and log files. Class *Utilities* contains
    some general-use functionalities. Finally classes *MolecularDynamics*
    and *MonteCarlo* are the core Molecular Dynamics and Monte Carlo
    codes, respectively.

Start coding
-------------

.. container:: justify

    Inside a dedicated folder, create 5 blank Python scripts named respectively:

    - *InitializeSimulation.py*
    - *Outputs.py*
    - *Utilities.py*
    - *MolecularDynamics.py*
    - *MonteCarlo.py*

.. container:: justify

    Let us start by filing some parts of the *InitializeSimulation*
    class within the *InitializeSimulation.py* file: 

.. code-block:: python

    class InitializeSimulation:
        def __init__(self,
                     *args,
                     **kwargs):
            super().__init__(*args, **kwargs) 

.. container:: justify

    Let us write something similar for the *Utilities* class 
    within the *Utilities.py* file:

.. code-block:: python

    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. container:: justify

    For the *Outputs* class, let us anticipate that the outputs
    from the code will be saved in a folder, which by default
    is named *results/*. If the folder does not exist, it will be
    created using *os.mkdir()*:

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

    Let us create the *MolecularDynamics* class which inherits
    the 3 previously defined classes.

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
    *MonteCarlo* class:

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

    from MonteCarlo import MonteCarlo
    from MolecularDynamics import MolecularDynamics

    md = MolecularDynamics(data_folder = "md-output/")
    md.run()
    mc = MonteCarlo(data_folder = "mc-output/")
    mc.run()

.. container:: justify

    If everything is working well two folders named *md-output/*
    and *mc-output/* must have been created, and no error message
    should appear.