.. _chapter5-label:

Output the simulation state
===========================

Here, the code is modified to allow us to follow the evolution of the system
during the simulation. To do so, the *Output* class is modified.

Update the MinimizeEnergy class
-------------------------------

Let us start by calling two additional methods within the for loop of the
*MinimizeEnergy* class, within the *run()* method.

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    def run(self):
        (...)
        for self.step in range(0, self.maximum_steps+1):
            (...)
                self.displacement *= 0.2 # Multiply the displacement by a factor 0.2
            log_simulation_data(self)
            update_dump_file(self, "dump.min.lammpstrj")

.. label:: end_MinimizeEnergy_class

The two methods named *update_log_minimize()* and *update_dump_file()*, are used
to print the information in the terminal and in a LAMMPS-type data file, respectively.
These two methods will be written in the following.

Create logger
-------------

Let us create functions named *log_simulation_data* to a file named *logger.py*.
With the logger, some output are being printed in a file, as well as in the terminal.
The frequency of printing is set by *thermo_period*, see :ref:`chapter3-label`.
All quantities are re-dimensionalized before getting outputed.

.. label:: start_logger_class

.. code-block:: python

    import os
    import logging

    # Function to set up the logger
    def setup_logger(folder_name, overwrite=False):
        # Create a custom logger
        logger = logging.getLogger('simulation_logger')
        logger.setLevel(logging.INFO)
        logger.propagate = False  # Disable propagation to prevent double logging

        # Clear any existing handlers if this function is called again
        if logger.hasHandlers():
            logger.handlers.clear()

        # Create handlers for console and file
        console_handler = logging.StreamHandler()  # To log to the terminal
        log_file_path = os.path.join(folder_name, 'simulation.log')

        # Use 'w' mode to overwrite the log file if overwrite is True, else use 'a' mode to append
        file_mode = 'w' if overwrite else 'a'
        file_handler = logging.FileHandler(log_file_path, mode=file_mode)  # To log to a file

        # Create formatters and add them to the handlers
        formatter = logging.Formatter('%(message)s')
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)

        # Add handlers to the logger
        logger.addHandler(console_handler)
        logger.addHandler(file_handler)

        return logger

    def log_simulation_data(code):

        # Setup the logger with the folder name, overwriting the log if code.step is 0
        logger = setup_logger(code.data_folder, overwrite=(code.step == 0))

        # Logging the simulation data
        if code.thermo_period is not None:
            if code.step % code.thermo_period == 0:
                if code.step == 0:
                    Epot = code.compute_potential() \
                        * code.ref_energy  # kcal/mol
                else:
                    Epot = code.Epot * code.ref_energy  # kcal/mol
                if code.step == 0:
                    if code.thermo_outputs == "Epot":
                        logger.info(f"step Epot")
                    elif code.thermo_outputs == "Epot-MaxF":
                        logger.info(f"step Epot MaxF")
                    elif code.thermo_outputs == "Epot-press":
                        logger.info(f"step Epot press")
                if code.thermo_outputs == "Epot":
                    logger.info(f"{code.step} {Epot.magnitude:.2f}")
                elif code.thermo_outputs == "Epot-MaxF":
                    logger.info(f"{code.step} {Epot.magnitude:.2f} {code.MaxF:.2f}")
                elif code.thermo_outputs == "Epot-press":
                    code.calculate_pressure()
                    press = code.pressure * code.ref_pressure  # Atm
                    logger.info(f"{code.step} {Epot.magnitude:.2f} {press.magnitude:.2f}")

.. label:: end_logger_class

Create dumper
-------------

Let us create a function named *update_dump_file* to a file named 
*dumper.py*. The dumper will print a *.lammpstrj file*, which contains the box
dimensions and atom positions at every chosen frame (set by *dumping_period*,
see :ref:`chapter3-label`). All quantities are dimensionalized before getting outputed, and the file follows
a LAMMPS dump format, and can be read by molecular dynamics softwares like VMD.

.. label:: start_dumper_class

.. code-block:: python

    import numpy as np

    def update_dump_file(code, filename, velocity=False):
        if code.dumping_period is not None:
            if code.step % code.dumping_period == 0:
                # Convert units to the *real* dimensions
                box_boundaries = code.box_boundaries*code.ref_length # Angstrom
                atoms_positions = code.atoms_positions*code.ref_length # Angstrom
                atoms_types = code.atoms_type
                if velocity:
                    atoms_velocities = code.atoms_velocities \
                        * code.ref_length/code.ref_time # Angstrom/femtosecond
                # Start writting the file
                if code.step == 0: # Create new file
                    f = open(code.data_folder + filename, "w")
                else: # Append to excisting file
                    f = open(code.data_folder + filename, "a")
                f.write("ITEM: TIMESTEP\n")
                f.write(str(code.step) + "\n")
                f.write("ITEM: NUMBER OF ATOMS\n")
                f.write(str(np.sum(code.number_atoms)) + "\n")
                f.write("ITEM: BOX BOUNDS pp pp pp\n")
                for dim in np.arange(3):
                    f.write(str(box_boundaries[dim][0].magnitude) + " "
                            + str(box_boundaries[dim][1].magnitude) + "\n")
                cpt = 1
                if velocity:
                    f.write("ITEM: ATOMS id type x y z vx vy vz\n")
                    characters = "%d %d %.3f %.3f %.3f %.3f %.3f %.3f %s"
                    for type, xyz, vxyz in zip(atoms_types,
                                            atoms_positions.magnitude,
                                            atoms_velocities.magnitude):
                        v = [cpt, type, xyz[0], xyz[1], xyz[2],
                                vxyz[0], vxyz[1], vxyz[2]]
                        f.write(characters % (v[0], v[1], v[2], v[3], v[4],
                                            v[5], v[6], v[7], '\n'))
                        cpt += 1
                else:
                    f.write("ITEM: ATOMS id type x y z\n")
                    characters = "%d %d %.3f %.3f %.3f %s"
                    for type, xyz in zip(atoms_types,
                                        atoms_positions.magnitude):
                        v = [cpt, type, xyz[0], xyz[1], xyz[2]]
                        f.write(characters % (v[0], v[1], v[2],
                                            v[3], v[4], '\n'))
                        cpt += 1
                f.close()  

.. label:: end_dumper_class

Import the functions
--------------------

The Monte Carlo and the Minimize class must both import *update_dump_file*
and *log_simulation_data*. Add these lines at the top of the *MonteCarlo.py* file:

.. label:: start_MonteCarlo_class

.. code-block:: python

    from dumper import update_dump_file
    from logger import log_simulation_data

.. label:: end_MonteCarlo_class

Add the same lines at the top of the *MinimizeEnergy.py* file:

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    from dumper import update_dump_file
    from logger import log_simulation_data

.. label:: end_MinimizeEnergy_class

Finally, let us make sure that *thermo_period*, *dumping_period*, and *thermo_outputs*
parameters are passed the InitializeSimulation method:

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
            (...)
                neighbor=1, # Integer
                thermo_period = None,
                dumping_period = None,
                thermo_outputs = None,
                data_folder="Outputs/",

.. label:: end_InitializeSimulation_class

.. label:: start_InitializeSimulation_class

.. code-block:: python

    def __init__(self,
        (...)
        self.initial_positions = initial_positions
        self.thermo_period = thermo_period
        self.dumping_period = dumping_period
        self.thermo_outputs = thermo_outputs
        self.data_folder = data_folder
        if os.path.exists(self.data_folder) is False:
            os.mkdir(self.data_folder)

.. label:: end_InitializeSimulation_class

Test the code
-------------

One can use a test similar as the previous ones. Let us ask out code to print
information in the dump and the log files, and then let us assert the
files were indeed created without the *Outputs/* folder:

.. label:: start_test_5a_class

.. code-block:: python

    from MinimizeEnergy import MinimizeEnergy
    from pint import UnitRegistry
    ureg = UnitRegistry()
    import os

    # Define atom number of each group
    nmb_1, nmb_2= [2, 3]
    # Define LJ parameters (sigma)
    sig_1, sig_2 = [3, 4]*ureg.angstrom
    # Define LJ parameters (epsilon)
    eps_1, eps_2 = [0.2, 0.4]*ureg.kcal/ureg.mol
    # Define atom mass
    mss_1, mss_2 = [10, 20]*ureg.gram/ureg.mol
    # Define box size
    L = 20*ureg.angstrom
    # Define a cut off
    rc = 2.5*sig_1

    # Initialize the prepare object
    minimizer = MinimizeEnergy(
        ureg = ureg,
        maximum_steps=100,
        thermo_period=25,
        dumping_period=25,
        number_atoms=[nmb_1, nmb_2],
        epsilon=[eps_1, eps_2], # kcal/mol
        sigma=[sig_1, sig_2], # A
        atom_mass=[mss_1, mss_2], # g/mol
        box_dimensions=[L, L, L], # A
        cut_off=rc,
        data_folder="Outputs/",
        thermo_outputs="Epot-MaxF",
    )
    minimizer.run()

    # Test function using pytest
    def test_output_files():
        assert os.path.exists("Outputs/dump.min.lammpstrj"), "Test failed: dump file was not created"
        assert os.path.exists("Outputs/simulation.log"), "Test failed: log file was not created"
        print("Test passed")

    # If the script is run directly, execute the tests
    if __name__ == "__main__":
        import pytest
        # Run pytest programmatically
        pytest.main(["-s", __file__])

.. label:: end_test_5a_class

I addition to the files getting created, information must be printed in the terminal
during the simulation:

.. code-block:: bw

    step Epot MaxF
    0 -0.17 1.93
    25 -1.08 1.81
    50 -1.11 1.42
    75 -1.22 3.77
    100 -2.10 1.28

The data from the *simulation.log* can be used to generate plots using softwares
line XmGrace, GnuPlot, or Python/Pyplot. For the later, one can use a simple data
reader to import the data from *Outputs/simulation.log* into Python. Copy the
following lines in a file named *reader.py*:

.. label:: start_reader_class

.. code-block:: python

    import csv

    def import_data(file_path):
        """
        Imports a data file with a variable number of columns into a list
        of numerical arrays. The first line (header) is read as a string.

        Parameters:
        - file_path (str): Path to the data file.

        Returns:
        - header (str): The header line as a string.
        - data (list of lists): List where each sublist contains the numeric values of a row.
        """
        data = []
        header = ""
        with open(file_path, mode='r') as file:
            # Read the header as a string
            header = file.readline().strip()
            # Use csv.reader to process the remaining lines
            reader = csv.reader(file, delimiter=' ')
            for row in reader:
                # Filter out empty fields resulting from multiple spaces
                filtered_row = [float(value) for value in row if value]
                data.append(filtered_row)
        return header, data

.. label:: end_reader_class

The *import_data* function from *reader.py* can simply be used as follows:

.. label:: start_test_5b_class

    from reader import import_data

    file_path = "Outputs/simulation.log"
    header, data = import_data(file_path)

    print(header)
    for row in data:
        print(row)

.. label:: end_test_5b_class

Which must return:

.. code-block:: bw

    step Epot MaxF
    [0.0, 9.48, 1049.12]
    [25.0, -2.12, 1.22]
    [50.0, -2.19, 2.85]
    [75.0, -2.64, 0.99]
    [100.0, -2.64, 0.51]