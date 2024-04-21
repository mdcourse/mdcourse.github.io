Code structure
==============

.. container:: justify

    Here, the files that constitute the basis of the code
    are created. These files will be progressively filled. 

Presentation
------------

.. container:: justify

    The Python code performing the molecular simulations will be
    divided into seven classes:

    - *Prepare.py --* Methods preparing the non-dimensionalization of the units
    - *Utilities --* Methods of general purpose, inherited by all the other classes
    - *Outputs --* Methods of interest for printing information in log or data files, inherited by all the classes except *Utilities*
    - *InitializeSimulation --* Methods necessary to set up the system and prepare the simulation, inherited by all the classes except *Outputs* and *Utilities*
    - *MinimizeEnergy --* Methods for performing energy minimization, including 
    - *MonteCarlo --* Methods for performing Monte Carlo simulation in different ensembles (Grand canonical, canonical)
    - *MolecularDynamics --* Methods for performing molecular dynamics in different ensembles (NVE, NPT, NVT)

.. container:: justify

    Two additional files named *Potentials.py* and *Functions.py* will contain
    some additional files.

Start coding
-------------

.. container:: justify

    Inside a dedicated folder, create 9 blank Python scripts named respectively:

    - *Potentials.py*
    - *Prepare.py*
    - *Utilities.py*
    - *InitializeSimulation.py*
    - *Measurements.py*
    - *Outputs.py*
    - *MinimizeEnergy.py*
    - *MolecularDynamics.py*
    - *MonteCarlo.py*

Potential for inter-atomic interaction
--------------------------------------

.. container:: justify

    The first file name *Potentials.py* contains the only potential that is
    currently implemented: the Lennard-Jones potential (LJ).

.. label:: start_Potentials_class

.. code-block:: python

    def LJ_potential(epsilon, sigma, r, derivative = False):
        if derivative:
            return 48*epsilon*((sigma/r)**12-0.5*(sigma/r)**6)/r
        else:
            return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

.. label:: end_Potentials_class

.. container:: justify

    Depending on the value of *derivative*, which can be either *False* or *True*,
    the *LJ_potential()* function will return the LJ force

.. math::

    F_\text{LJ} = 48 \dfrac{\epsilon}{r} \left[ \left( \frac{\sigma}{r} \right)^{12}- \frac{1}{2} \left( \frac{\sigma}{r} \right)^6 \right],

.. container:: justify

    or the LJ potential

.. math::

    U_\text{LJ} = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{12}- \left( \frac{\sigma}{r} \right)^6 \right].

Create the classes
------------------

.. container:: justify

    Let us create all the classes and their inheritance. The classes will be
    filled progressively during the following chapters.

.. container:: justify

    The first class is the *Prepare* class which will serve the
    nondimensionalization of all the parameters. Within the *Prepare.py* file,
    copy the following lines:

.. label:: start_Prepare_class

.. code-block:: python

    class Prepare:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_Prepare_class

.. container:: justify

    The second class is named *Utilities*. Within the *Utilities.py* file,
    copy the following lines:

.. label:: start_Utilities_class

.. code-block:: python

    from Potentials import LJ_potential


    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_Utilities_class

.. container:: justify

    The *InitializeSimulation* class inherits the *Prepare* class. Within the
    *InitializeSimulation.py* file, copy the following lines:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    import numpy as np
    from Prepare import Prepare


    class InitializeSimulation(Prepare):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)

.. label:: end_InitializeSimulation_class

.. container:: justify

    The *Measurements* class inherits both *InitializeSimulation*  and
    *Utilities* classes. Within the *Measurements.py* file, copy the following lines:

.. label:: start_Measurements_class

.. code-block:: python

    from InitializeSimulation import InitializeSimulation
    from Utilities import Utilities


    class Measurements(InitializeSimulation, Utilities):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)
          
.. label:: end_Measurements_class

.. container:: justify

    The *Outputs* class inherits the *Measurements* class. Within the
    *Outputs.py* file, copy the following lines:

.. label:: start_Outputs_class

.. code-block:: python

    import os
    from Measurements import Measurements


    class Outputs(Measurements):
        def __init__(self,
                    data_folder="Outputs/",
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)
            self.data_folder = data_folder
            if os.path.exists(self.data_folder) is False:
                os.mkdir(self.data_folder)

.. label:: end_Outputs_class

.. container:: justify

    Here, we anticipate that the outputs
    from the code will be saved in a folder, which by default
    is named *results/*. If the folder does not exist, it will be
    created using *os.mkdir()* from the *os* module, which was previously
    imported.

.. container:: justify

    Finally, let us create the three remaining classes, named respectively *MinimizeEnergy*,
    *MonteCarlo*, and *MolecularDynamics*. Each class inherits
    the *Outputs* class. Within the *MinimizeEnergy.py* file, copy the
    following lines:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    from Outputs import Outputs


    class MinimizeEnergy(Outputs):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MinimizeEnergy_class

.. container:: justify

    Within the *MonteCarlo.py* file, copy the following lines:

.. label:: start_MonteCarlo_class

.. code-block:: python

    from Outputs import Outputs


    class MonteCarlo(Outputs):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MonteCarlo_class

.. container:: justify

    Finally, within the *MolecularDynamics.py* file, copy the following lines:

.. label:: start_MolecularDynamics_class

.. code-block:: python

    from Outputs import Outputs

    class MolecularDynamics(Outputs):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)

.. label:: end_MolecularDynamics_class

Test the code
-------------

.. container:: justify

    We can create a simple test to ensure that the classes
    are being inherited as expected. Within the same folder,
    create a new Jupyter notebook called *test.ipynb*, and copy
    the following lines into it:

.. label:: start_test_First_class

.. code-block:: python

    from MonteCarlo import MonteCarlo
    from MolecularDynamics import MolecularDynamics

    md = MolecularDynamics(data_folder = "md-output/")
    md.__init__()
    mc = MonteCarlo(data_folder = "mc-output/")
    mc.__init__()

.. label:: end_test_First_class

.. container:: justify

    If everything is working well two folders named *md-output/*
    and *mc-output/* must have been created, and no error message
    should appear.