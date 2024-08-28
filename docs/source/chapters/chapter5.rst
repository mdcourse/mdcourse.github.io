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
                self.displacement *= 0.2
            log_simulation_data(self)
            update_dump_file(self, "dump.min.lammpstrj")

.. label:: end_MinimizeEnergy_class

The two methods named *update_log_minimize()* and *update_dump_file()*, are used
to print the information in the terminal and in a LAMMPS-type data file, respectively.
These two methods will be written in the following.

Update the dump and log
-----------------------

Let us add the following functions named *update_dump_file* and *log_simulation_data*
to the *tools.py* file. Here, some variable are being printed in a file, such as box dimension and atom positions.
All quantities are dimensionalized before getting outputed, and the file follows
a LAMMPS dump format, and can be read by molecular dynamics softwares like VMD.

.. label:: start_tools_class

.. code-block:: python

    import numpy as np

    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s'
    )

    # Create a custom logger
    logger = logging.getLogger('simulation_logger')
    logger.setLevel(logging.INFO)
    # Disable propagation to prevent double logging
    logger.propagate = False

    console_handler = logging.StreamHandler()  # To log to the terminal
    file_handler = logging.FileHandler('simulation.log')  # To log to a file
    formatter = logging.Formatter('%(message)s')
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    def update_dump_file(code, filename, velocity=False):
        if code.dumping_period is not None:
            if code.step % code.dumping_period == 0:
                # Convert units to the *real* dimensions
                box_boundaries = code.box_boundaries\
                    * code.reference_distance # Angstrom
                atoms_positions = code.atoms_positions\
                    * code.reference_distance # Angstrom
                atoms_types = code.atoms_type
                if velocity:
                    atoms_velocities = code.atoms_velocities \
                        * code.reference_distance/code.reference_time # Angstrom/femtosecond
                # Start writting the file
                if code.step == 0: # Create new file
                    f = open(code.data_folder + filename, "w")
                else: # Append to excisting file
                    f = open(code.data_folder + filename, "a")
                f.write("ITEM: TIMESTEP\n")
                f.write(str(code.step) + "\n")
                f.write("ITEM: NUMBER OF ATOMS\n")
                f.write(str(code.total_number_atoms) + "\n")
                f.write("ITEM: BOX BOUNDS pp pp pp\n")
                for dim in np.arange(code.dimensions):
                    f.write(str(box_boundaries[dim][0]) + " "
                            + str(box_boundaries[dim][1]) + "\n")
                cpt = 1
                if velocity:
                    f.write("ITEM: ATOMS id type x y z vx vy vz\n")
                    characters = "%d %d %.3f %.3f %.3f %.3f %.3f %.3f %s"
                    for type, xyz, vxyz in zip(atoms_types,
                                            atoms_positions,
                                            atoms_velocities):
                        v = [cpt, type, xyz[0], xyz[1], xyz[2],
                                vxyz[0], vxyz[1], vxyz[2]]
                        f.write(characters % (v[0], v[1], v[2], v[3], v[4],
                                            v[5], v[6], v[7], '\n'))
                        cpt += 1
                else:
                    f.write("ITEM: ATOMS id type x y z\n")
                    characters = "%d %d %.3f %.3f %.3f %s"
                    for type, xyz in zip(atoms_types,
                                        atoms_positions):
                        v = [cpt, type, xyz[0], xyz[1], xyz[2]]
                        f.write(characters % (v[0], v[1], v[2],
                                            v[3], v[4], '\n'))
                        cpt += 1
                f.close()

    def log_simulation_data(code):
        if code.thermo_period is not None:
            if code.step % code.thermo_period == 0:
                try:
                    Epot = code.Epot * code.reference_energy  # kcal/mol
                except:
                    Epot = code.compute_potential(output="potential") \
                        * code.reference_energy  # kcal/mol
                if code.step == 0:
                    if code.thermo_outputs == "Epot":
                        logger.info(f"step Epot")
                    elif code.thermo_outputs == "Epot-MaxF":
                        logger.info(f"step Epot MaxF")
                    elif code.thermo_outputs == "Epot-press":
                        logger.info(f"step Epot press")
                if code.thermo_outputs == "Epot":
                    logger.info(f"{code.step} {Epot:.2f}")
                elif code.thermo_outputs == "Epot-MaxF":
                    logger.info(f"{code.step} {Epot:.2f} {code.MaxF:.2f}")
                elif code.thermo_outputs == "Epot-press":
                    code.calculate_pressure()
                    press = code.pressure \
                        * code.reference_pressure  # Atm
                    logger.info(f"{code.step} {Epot:.2f} {press:.2f}")    

.. label:: end_tools_class

Import the functions
--------------------

The Monte Carlo and the Minimize class must import *update_dump_file* and *log_simulation_data*.

.. label:: start_MonteCarlo_class

.. code-block:: python

    from tools import update_dump_file, log_simulation_data

.. label:: end_MonteCarlo_class

and

.. label:: start_MinimizeEnergy_class

.. code-block:: python

    from tools import update_dump_file, log_simulation_data

.. label:: end_MinimizeEnergy_class

Test the code
-------------

One can use the same test as previously, and ask the code to print information
every 10 steps in the dump files, as well as in the log:

.. label:: start_test_Outputs_class

.. code-block:: python

    import os
    from MinimizeEnergy import MinimizeEnergy

    min = MinimizeEnergy(maximum_steps=100,
        thermo_period=10,
        dumping_period=10,
        thermo_outputs = "Epot-MaxF",
        number_atoms=[2, 3],
        epsilon=[0.1, 1.0], # kcal/mol
        sigma=[3, 6], # A
        atom_mass=[1, 1], # g/mol
        box_dimensions=[20, 20, 20], # A
        data_folder="Outputs/",
        )
    min.run()

    assert os.path.exists("Outputs/dump.min.lammpstrj")

.. label:: end_test_Outputs_class

When running the simulation, information must be printed in the terminal:

.. code-block:: bw

    step Epot (kcal/mol) MaxF (kcal/A/mol)
    0, -1.40, 14.19
    10, -1.80, 17.98
    20, -2.37, 2.65
    30, -2.50, 4.00
    (...)

and a file named *dump.min.lammpstrj* must have appeared in the *Outputs/* folder.