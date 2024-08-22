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
            self.update_log_minimize(Epot, max_forces)
            self.update_dump_file(filename="dump.min.lammpstrj")

.. label:: end_MinimizeEnergy_class

The two methods named *update_log_minimize()* and *update_dump_file()*, are used
to print the information in the terminal and in a LAMMPS-type data file, respectively.
These two methods will be written in the following.

Improve the output class
------------------------

Modify the beginning of the *Outputs.py* file to import NumPy and the
constant module of SciPy:

.. code-block:: python
        
    from scipy import constants as cst
    import numpy as np
    import os
    from Measurements import Measurements


    class Outputs(Measurements):
        def __init__(self,
                    thermo_period=None,
                    dumping_period=None,
                    data_folder="Outputs/",
                    *args,
                    **kwargs):
            self.thermo_period = thermo_period
            self.dumping_period = dumping_period
            self.data_folder = data_folder
            super().__init__(*args, **kwargs)
            if os.path.exists(self.data_folder) is False:
                os.mkdir(self.data_folder)

Here, two additional variables were added: *thermo_period* which controls
the period at which information is printed during the run, and *dumping_period*
which controls the period at which atom positions are printed in the dump
file. 

Update the dump file
--------------------

Finally, add the following method named *update_dump_file()* to the
*Output* class. 

.. code-block:: python

    def update_dump_file(self, filename, velocity=False):
        if self.dumping_period is not None:
            box_boundaries = self.box_boundaries\
                * self.reference_distance
            atoms_positions = self.atoms_positions\
                * self.reference_distance
            atoms_types = self.atoms_type
            if velocity:
                atoms_velocities = self.atoms_velocities \
                    * self.reference_distance/self.reference_time
            if self.step % self.dumping_period == 0:
                if self.step == 0:
                    f = open(self.data_folder + filename, "w")
                else:
                    f = open(self.data_folder + filename, "a")
                f.write("ITEM: TIMESTEP\n")
                f.write(str(self.step) + "\n")
                f.write("ITEM: NUMBER OF ATOMS\n")
                f.write(str(self.total_number_atoms) + "\n")
                f.write("ITEM: BOX BOUNDS pp pp pp\n")
                for dim in np.arange(self.dimensions):
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

Update the log file
--------------------

Finally, add the following method to the *Output* class. 

.. code-block:: python

    def update_log_minimize(self, Epot, maxForce):
        if (self.thermo_period is not None):
            if ((self.step % self.thermo_period == 0)
                    | (self.thermo_period == 0)):
                epot_kcalmol = Epot * self.reference_energy
                max_force_kcalmolA = maxForce \
                    * self.reference_energy / self.reference_distance
                if self.step == 0:
                    characters = "%s %s %s"
                    print(characters % ("step",
                                        "epot",
                                        "maxF"))
                characters = "%d %.3f %.3f"
                print(characters % (self.step,
                                    epot_kcalmol,
                                    max_force_kcalmolA))
