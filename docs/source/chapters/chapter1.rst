Start coding
============

Presentation
------------

To perform the molecular simulations, three key steps will be followed here:

- First, the system will be initialized. This includes creating the simulation
  box, placing the atoms, and choosing the potential and parameters for their
  interactions.
- Second, an energy minimization of the system will be performed to place the
  atoms in reasonable relative positions and reduce the energy of the system.
- Third, the main algorithm will be executed: either molecular dynamics or
  Monte Carlo.

In addition to these three key steps, some additional tasks will be performed, 
such as units non-dimensionalization, data measurements, and the outputting of 
thermodynamic information during the simulation.

For better readability, the present code is split into separate files, with
each file containing either Python funtions or Python classes.

List of files
-------------

In total, 7 Python files containing classes will be written during this
course, as well as 2 files containing functions:


.. list-table::
   :widths: 40 60
   :header-rows: 1
   
   * - Fine name 
     - Content
   * - *Prepare.py* 
     - *Prepare* class
   * - *Utilities.py* 
     - *Utilities* class
   * - *InitializeSimulation.py*
     - *InitializeSimulation* class
   * - *Outputs.py*
     - *Outputs* class
   * - *MinimizeEnergy.py* 
     - *MinimizeEnergy* class
   * - *MolecularDynamics.py*
     - *MolecularDynamics* class
   * - *MonteCarlo.py*
     - *MonteCarlo* class
   * - *Measurements.py* 
     - Functions
   * - *Potentials.py* 
     - Functions

Each of these files will serve to perform specific tasks. Within these files,
the final Python code will be divided into seven classes:

- *Prepare --* Methods preparing the non-dimensionalization of the units
- *Utilities --* Methods of general purpose, inherited by all the other classes
- *Outputs --* Methods of interest for printing information in log or data
  files, inherited by all the classes except *Utilities*
- *InitializeSimulation --* Methods necessary to set up the system and prepare
  the simulation, inherited by all the classes except *Outputs* and *Utilities*
- *MinimizeEnergy --* Methods for performing energy minimization, including 
- *MonteCarlo --* Methods for performing Monte Carlo simulation in different
  ensembles (Grand canonical, canonical)
- *MolecularDynamics --* Methods for performing molecular dynamics in
  different ensembles (NVE, NPT, NVT)

Potential for inter-atomic interaction
--------------------------------------

In molecular simulations, potential functions are used to mimick the
interaction between atoms. Although there exists some more complicated
options, potentials are usually defined as functions of the
distance :math:`r` between two atoms. 

Within a dedicated folder, create the first file named *Potentials.py*. This
file will contain a functon named *LJ_potential* for the Lennard-Jones
potential (LJ). Copy the following into *Potentials.py*:

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

Let us create all the classes and their inheritance. The classes will be
filled progressively during the following chapters.

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

The *Outputs* class inherits the *Measurements* class. Within the
*Outputs.py* file, copy the following lines:

.. label:: start_Outputs_class

.. code-block:: python

    from scipy import constants as cst
    import numpy as np
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

Here, we anticipate that the outputs
from the code will be saved in a folder, which by default
is named *results/*. If the folder does not exist, it will be
created using *os.mkdir()* from the *os* module, which was previously
imported. NumPy and the *constants* module of SciPy were also imported.

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

We can create a simple test to ensure that the classes
are being inherited as expected. Within the same folder,
create a new Jupyter notebook called *test.ipynb*, and copy
the following lines into it:

.. label:: start_test_First_class

.. code-block:: python

    import os
    from MonteCarlo import MonteCarlo
    from MolecularDynamics import MolecularDynamics

    md = MolecularDynamics(data_folder = "md-output/")
    md.__init__()
    mc = MonteCarlo(data_folder = "mc-output/")
    mc.__init__()

    assert os.path.exists("mc-output"), """Error, missing mc-output/ folder"""
    assert os.path.exists("md-output"), """Error, missing md-output/ folder"""
    assert os.path.exists("Outputs"), """Error, missing Outputs/ folder"""

.. label:: end_test_First_class

If everything is working well, 3 folders named respectively *md-output/*, *mc-output/*,
and *Outputs/* must have been created. If not, running the test will generate an
error message.
