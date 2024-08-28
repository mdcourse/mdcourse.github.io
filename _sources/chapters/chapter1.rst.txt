Start coding
============

Presentation
------------

To perform the molecular simulations, three key steps will be followed:

1. **System Initialization:** This step involves creating the simulation box,
   placing the atoms, and selecting the parameters for their interactions.
  
2. **Energy Minimization:** An energy minimization of the system will be
   carried out to position the atoms in reasonable relative locations and
   reduce the system's energy.
  
3. **Execution of the Main Algorithm:** The main algorithm, either Monte Carlo
   or molecular dynamics, will then be executed.

In addition to these three key steps, some supplementary tasks will be
performed, such as non-dimensionalizing units, measuring data, and outputting
thermodynamic information during the simulation.

List of Files
-------------

For better readability, the code is split into separate files, with each file
containing either Python functions or classes:

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - File Name 
     - Content
   * - *Prepare.py* 
     - *Prepare* class: Methods for preparing the non-dimensionalization of the
       units
   * - *Utilities.py* 
     - *Utilities* class: General-purpose methods, inherited by all other classes
   * - *InitializeSimulation.py*
     - *InitializeSimulation* class: Methods necessary to set up the system and
       prepare the simulation, inherited by all the classes below
   * - *MinimizeEnergy.py* 
     - *MinimizeEnergy* class: Methods for performing energy minimization
   * - *MonteCarlo.py*
     - *MonteCarlo* class: Methods for performing Monte Carlo simulations in
       different ensembles (e.g., Grand Canonical, Canonical)
   * - *MolecularDynamics.py*
     - *MolecularDynamics* class: Methods for performing molecular dynamics in
       different ensembles (NVE, NPT, NVT)
   * - *measurements.py* 
     - Functions for performing specific measurements on the system
   * - *potentials.py* 
     - Functions for calculating the potentials and forces between atoms
   * - *tools.py*
     - Functions for outputting data into text files


Potential for inter-atomic interaction
--------------------------------------

In molecular simulations, potential functions are used to mimic the interaction
between atoms. Although more complicated options exist, potentials are usually
defined as functions of the distance :math:`r` between two atoms.

Within a dedicated folder, create the first file named *potentials.py*. This
file will contain a function named *LJ_potential* for the Lennard-Jones
potential (LJ). Copy the following into *potentials.py*:

.. label:: start_potentials_class

.. code-block:: python

    def LJ_potential(epsilon, sigma, r, derivative = False):
        if derivative:
            return 48*epsilon*((sigma/r)**12-0.5*(sigma/r)**6)/r
        else:
            return 4*epsilon*((sigma/r)**12-(sigma/r)**6)

.. label:: end_potentials_class

Depending on the value of the optional argument *derivative*, which can be
either *False* or *True*, the *LJ_potential* function will return the force:

.. math::

    F_\text{LJ} = 48 \dfrac{\epsilon}{r} \left[ \left( \frac{\sigma}{r} \right)^{12}- \frac{1}{2} \left( \frac{\sigma}{r} \right)^6 \right],

or the potential energy:

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

    from potentials import LJ_potential


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

Finally, let us create the three remaining classes, named respectively *MinimizeEnergy*,
*MonteCarlo*, and *MolecularDynamics*. Each class inherits
the *Measurements* class. Within the *MinimizeEnergy.py* file, copy the
following lines:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    from Measurements import Measurements
    import os


    class MinimizeEnergy(Measurements):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MinimizeEnergy_class

Within the *MonteCarlo.py* file, copy the following lines:

.. label:: start_MonteCarlo_class

.. code-block:: python

    from scipy import constants as cst
    import numpy as np
    import copy
    import os
    from Measurements import Measurements

    import warnings
    warnings.filterwarnings('ignore')


    class MonteCarlo(Measurements):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MonteCarlo_class

The *warnings* was placed to avoid the anoying message "*RuntimeWarning: overflow
encountered in exp*" that is sometimes triggered by the exponential of the
*acceptation_probability* (see :ref:`chapter6-label`).

Finally, within the *MolecularDynamics.py* file, copy the following lines:

.. label:: start_MolecularDynamics_class

.. code-block:: python

    import numpy as np
    from Measurements import Measurements


    class MolecularDynamics(Measurements):
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

    md = MolecularDynamics()
    md.__init__()
    mc = MonteCarlo()
    mc.__init__()

    #assert os.path.exists("mc-output"), """Error, missing mc-output/ folder"""
    #assert os.path.exists("md-output"), """Error, missing md-output/ folder"""

.. label:: end_test_First_class

No error should be returned.