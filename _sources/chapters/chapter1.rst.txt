.. _chapter1-label:

Start Coding
============

Let's start with the Python script. During this first chapter, some Python
files will be created and filled in a minimal fashion. At the end of this
chapter, a small test will be set up to ensure that the files were correctly
created.

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
   * - Prepare.py 
     - *Prepare* class: Methods for preparing the non-dimensionalization of the
       units
   * - Utilities.py
     - *Utilities* class: General-purpose methods, inherited by all other
       classes
   * - InitializeSimulation.py
     - *InitializeSimulation* class: Methods necessary to set up the system and
       prepare the simulation, inherited by all the classes below
   * - MinimizeEnergy.py
     - *MinimizeEnergy* class: Methods for performing energy minimization
   * - MonteCarlo.py
     - *MonteCarlo* class: Methods for performing Monte Carlo simulations in
       different ensembles (e.g., Grand Canonical, Canonical)
   * - MolecularDynamics.py
     - *MolecularDynamics* class: Methods for performing molecular dynamics in
       different ensembles (NVE, NPT, NVT)
   * - Measurements.py
     - *Measurements* class: Methods for for performing specific measurements on the system
   * - *potentials.py* 
     - Functions for calculating the potentials and forces between atoms
   * - logger.py
     - Functions for outputting data into text files
   * - dumper.py
     - Functions for outputting data into trajectory files for visualization
   * - reader.py
     - Functions for importing data from text files

Some of these files are created in this chapter; others will be created later
on. All of these files must be created within the same folder.

Potential for Inter-Atomic Interaction
--------------------------------------

In molecular simulations, potential functions are used to model the interaction
between atoms. Although more complex options exist, potentials are usually
defined as functions of the distance :math:`r` between two atoms.

Create a file named *potentials.py*. This file will contain a function called
*potentials*. For now, the only potential that can be returned by this function
is the Lennard-Jones (LJ) potential, but this may change in the future.

Copy the following lines into *potentials.py*:

.. label:: start_potentials_class

.. code-block:: python

    def potentials(epsilon, sigma, r, derivative=False):
        if derivative:
            return 48 * epsilon * ((sigma / r) ** 12 - 0.5 * (sigma / r) ** 6) / r
        else:
            return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

.. label:: end_potentials_class

Depending on the value of the optional argument *derivative*, which can be
either *False* or *True*, this function returns the derivative of the potential,
i.e., the force, :math:`F_\text{LJ} = - \mathrm{d} U_\text{LJ} / \mathrm{d} r`:

.. math::

    F_\text{LJ} = 48 \dfrac{\epsilon}{r} \left[ \left( \frac{\sigma}{r} \right)^{12}
    - \frac{1}{2} \left( \frac{\sigma}{r} \right)^6 \right], ~ \text{for} ~ r < r_\text{c},

or the potential energy:

.. math::

    U_\text{LJ} = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{12}
    - \left( \frac{\sigma}{r} \right)^6 \right], ~ \text{for} ~ r < r_\text{c}.

Here, :math:`\sigma` is the distance at which the potential :math:`U_\text{LJ}`
is zero, :math:`\epsilon` is the depth of the potential well, and
:math:`r_\text{c}` is a cutoff distance. For :math:`r > r_\text{c}`,
:math:`U_\text{LJ} = 0` and :math:`F_\text{LJ} = 0`.

Create the Classes
------------------

Let's create the files with some minimal details about the classes and their
inheritance. The classes will be developed progressively in the following
chapters.

The first class is the *Prepare* class, which will be used for the
nondimensionalization of the parameters. In the same folder as *potentials.py*,
create the *Prepare.py* file and copy the following lines into it:

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

    from potentials import potentials


    class Utilities:
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_Utilities_class

The line *from potentials import potentials* is used to import the
previously created *potentials* function.

Within the *InitializeSimulation.py* file, copy the following lines:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    import os
    import numpy as np
    from Prepare import Prepare
    from Utilities import Utilities


    class InitializeSimulation(Prepare, Utilities):
        def __init__(self,
                    *args,
                    **kwargs,
                    ):
            super().__init__(*args, **kwargs)

.. label:: end_InitializeSimulation_class

The *InitializeSimulation* class inherits from the previously created
*Prepare* and *Utilities* classes. Additionally, we anticipate that |NumPy|
will be required :cite:`harris2020array`. We also anticipate that the *os*
module, which provides a way to interact with the operating system, will
be required :cite:`Rossum2009Python3`.

.. |NumPy| raw:: html

   <a href="https://numpy.org/" target="_blank">NumPy</a>

Within the *Measurements.py* file, copy the following lines:

.. label:: start_Measurements_class

.. code-block:: python

    import numpy as np
    from InitializeSimulation import InitializeSimulation


    class Measurements(InitializeSimulation):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)
          
.. label:: end_Measurements_class

The *Measurements* class inherits from *InitializeSimulation* (and thus
also inherits from the *Prepare* and *Utilities* classes).

Finally, let us create the three remaining classes: *MinimizeEnergy*,
*MonteCarlo*, and *MolecularDynamics*. Each of these classes inherits
from the *Measurements* class (and thus also from the *Prepare*, *Utilities*,
and *InitializeSimulation* classes).

Within the *MinimizeEnergy.py* file, copy the following lines:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    from Measurements import Measurements
    import numpy as np
    import copy
    import os


    class MinimizeEnergy(Measurements):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MinimizeEnergy_class

The *copy* library, which provides functions to create shallow or deep copies of
objects, is imported, along with *NumPy* and *os*.

Within the *MonteCarlo.py* file, copy the following lines:

.. label:: start_MonteCarlo_class

.. code-block:: python

    import numpy as np
    import copy
    from Measurements import Measurements

    import warnings
    warnings.filterwarnings('ignore')


    class MonteCarlo(Measurements):
        def __init__(self,
                    *args,
                    **kwargs):
            super().__init__(*args, **kwargs)

.. label:: end_MonteCarlo_class

The *ignore warnings* commands are optional; they were added to avoid the
annoying message "*RuntimeWarning: overflow encountered in exp*" that is sometimes
triggered by the exponential function of *acceptation_probability* (see the
:ref:`chapter6-label` chapter).

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

We can create simple tests to ensure that the classes are inherited as
expected. Within the same folder, create a new Python file called *test_1a.py*,
and copy the following lines into it:

.. label:: start_test_1a_class

.. code-block:: python

    # Import the required modules
    from Utilities import Utilities
    from MonteCarlo import MonteCarlo

    # Make sure that MonteCarlo correctly inherits from Utilities
    def test_montecarlo_inherits_from_utilities():
        assert issubclass(MonteCarlo, Utilities), \
            "MonteCarlo should inherit from Utilities"
        print("MonteCarlo correctly inherits from Utilities")

    # Make sure that Utilities does not inherit from MonteCarlo
    def test_utilities_does_not_inherit_from_montecarlo():
        assert not issubclass(Utilities, MonteCarlo), \
            "Utilities should not inherit from MonteCarlo"
        print("Utilities does not inherit from MonteCarlo, as expected")

    # In the script is launched with Python, call Pytest
    if __name__ == "__main__":
        import pytest
        pytest.main(["-s", __file__])

.. label:: end_test_1a_class

When run with Python, this script should return the following messages without
any *AssertionError*:

.. code-block:: bw

    Utilities does not inherit from MonteCarlo, as expected
    MonteCarlo correctly inherits from Utilities

Alternatively, this test can also be launched using *Pytest* by typing in a terminal:

.. code-block:: bash

    pytest .

We can also test that calling the *__init__* method of the *MonteCarlo* class
does not return any error. In new Python file called *test_1b.py*, copy the
following lines:

.. label:: start_test_1b_class

.. code-block:: python

    # Import the MonteCarlo class
    from MonteCarlo import MonteCarlo

    # Define a function that try to call the *__init__()* method
    def test_init_method():
        try:
            MonteCarlo().__init__()  # Call the method
            print("Method call succeeded")
        except Exception as e:
            print(f"Method call raised an error: {e}")

    # In the script is launched with Python, call Pytest
    if __name__ == "__main__":
        import pytest
        pytest.main(["-s", __file__])

.. label:: end_test_1b_class

Running this second test with Python should return "Method call succeeded".