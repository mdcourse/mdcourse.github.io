from scipy import constants as cst
import numpy as np
import logging
import copy

# to remove the runtime warning
import warnings
warnings.filterwarnings('ignore')

class InitializeSimulation:
    def __init__(self,
                 number_atoms,
                 Lx,
                 dimensions=3,
                 Ly = None,
                 Lz = None,
                 epsilon=0.1,
                 sigma=1,
                 atom_mass=1,
                 seed=None,
                 desired_temperature=300,
                 *args,
                 **kwargs,
                 ):
        self.number_atoms = number_atoms
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = dimensions
        self.epsilon = epsilon
        self.sigma = sigma
        self.atom_mass = atom_mass
        self.seed = seed
        self.desired_temperature = desired_temperature

        if self.seed is not None:
            np.random.seed(self.seed)

        self.reference_distance = self.sigma
        self.reference_energy = self.epsilon
        self.reference_mass = self.atom_mass
        self.reference_time = np.sqrt((self.reference_mass/cst.kilo/cst.Avogadro)*(self.reference_distance*cst.angstrom)**2/(self.reference_energy*cst.calorie*cst.kilo/cst.Avogadro))/cst.femto

        self.Lx = self.nondimensionalise_units(self.Lx, "distance")
        self.Ly = self.nondimensionalise_units(self.Ly, "distance")
        self.Lz = self.nondimensionalise_units(self.Lz, "distance")
        self.epsilon = self.nondimensionalise_units(self.epsilon, "energy")
        self.sigma = self.nondimensionalise_units(self.sigma, "distance")
        self.atom_mass = self.nondimensionalise_units(self.atom_mass, "mass")
        self.desired_temperature = self.nondimensionalise_units(self.desired_temperature, "temperature")

        self.initialize_box()
        self.populate_box()
        self.give_velocity()
        self.write_lammps_data(filename="initial.data")

    def nondimensionalise_units(self, variable, type):
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo # kCal/mol/K
        if variable is not None:
            if type == "distance":
                variable /= self.reference_distance
            elif type == "energy":
                variable /= self.reference_energy
            elif type == "mass":
                variable /= self.reference_mass
            elif type == "temperature":
                variable /= self.reference_energy/kB
            elif type == "time":
                variable /= self.reference_time
            else:
                print("Unknown variable type", type)
        return variable

    def initialize_box(self):
        """Define box boundaries based on Lx, Ly, and Lz.

        If Ly or Lz or both are None, then Lx is used instead"""
        box_boundaries = np.zeros((self.dimensions, 2))
        for dim, L in zip(range(self.dimensions), [self.Lx, self.Ly, self.Lz]):
            if L is not None:
                box_boundaries[dim] = -L/2, L/2
            else:
                box_boundaries[dim] = -self.Lx/2, self.Lx/2
        box_boundaries = self.nondimensionalise_units(box_boundaries, "distance")
        self.box_boundaries = box_boundaries

    def populate_box(self):
        """Place atoms at random positions within the box."""
        atoms_positions = np.zeros((self.number_atoms, self.dimensions))
        for dim in np.arange(self.dimensions):
            atoms_positions[:, dim] = np.random.random(self.number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2    
        self.atoms_positions = atoms_positions

    def give_velocity(self):
        """Give velocity to atoms so that the initial temperature is the desired one."""
        atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
        for dim in np.arange(self.dimensions):  
            atoms_velocities[:, dim] = np.random.normal(size=self.number_atoms)
        atoms_velocities *= np.sqrt(self.desired_temperature/self.atom_mass/self.dimensions)
        self.atoms_velocities = atoms_velocities

    def wrap_in_box(self):
        for dim in np.arange(self.dimensions):
            out_ids = self.atoms_positions[:, dim] > self.box_boundaries[dim][1]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] -= np.diff(self.box_boundaries[dim])[0]
            out_ids = self.atoms_positions[:, dim] < self.box_boundaries[dim][0]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] += np.diff(self.box_boundaries[dim])[0]

    def write_lammps_data(self, filename="lammps.data"):
        """Write a LAMMPS data file containing atoms positions and velocities"""
        f = open(filename, "w")
        f.write('# LAMMPS data file \n\n')
        f.write(str(self.number_atoms)+' atoms\n')
        f.write('1 atom types\n')
        f.write('\n')
        for LminLmax, dim in zip(self.box_boundaries*self.reference_distance, ["x", "y", "z"]):
            f.write(str(LminLmax[0])+' '+str(LminLmax[1])+' '+dim+'lo ' + dim +  'hi\n')
        f.write('\n')
        f.write('Atoms\n')
        f.write('\n')
        cpt = 1
        for xyz in self.atoms_positions*self.reference_distance:
            f.write(str(cpt)+ ' 1 ' + str(xyz[0]) + ' ' + str(xyz[1]) + ' ' + str(xyz[2]) +'\n')
            cpt += 1
        f.write('\n')
        f.write('Velocities\n')
        f.write('\n')
        cpt = 1
        for vxyz in self.atoms_velocities*self.reference_distance/self.reference_time:
            f.write(str(cpt) + ' ' + str(vxyz[0]) + ' ' + str(vxyz[1]) + ' ' + str(vxyz[2]) +'\n')
            cpt += 1
        f.close()

class Outputs:
    def __init__(self,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)

    def log(self):
        # Evaluate temperature (in K)
        self.calculate_temperature()
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo # kCal/mol/K
        temperature = np.round(self.temperature*self.reference_energy/kB,2)
        # Evaluate Pressure (in atm)
        self.calculate_pressure()
        pressure = np.round(self.pressure*self.reference_energy/self.reference_distance**3*cst.calorie*cst.kilo/cst.Avogadro/cst.angstrom**3/cst.atm, 3)
        if self.step ==1:
            print("step temp press")
        print(self.step, temperature, pressure)

class MolecularDynamics(InitializeSimulation, Outputs):
    def __init__(self,
                maximum_steps,
                tau_temp = None,
                tau_press = None,
                cut_off = 10,
                time_step=1,
                *args,
                **kwargs,
                ):
        self.maximum_steps = maximum_steps
        self.tau_temp = tau_temp  
        self.tau_press = tau_press
        self.cut_off = cut_off
        self.time_step = time_step
        super().__init__(*args, **kwargs)

        self.cut_off = self.nondimensionalise_units(self.cut_off, "distance")
        self.time_step = self.nondimensionalise_units(self.time_step, "time")

    def run(self):
        """Perform the loop over time."""
        for self.step in range(1, self.maximum_steps+1):
            self.integrate_equation_of_motion()
            self.wrap_in_box()
            if self.tau_temp is not None:
                self.apply_berendsen_thermostat()
            if self.tau_press is not None:
                self.apply_berendsen_barostat()
            self.log()
        self.write_lammps_data(filename="final.data")

    def evaluate_LJ_force(self):
        """Evaluate force based on LJ potential derivative."""
        forces = np.zeros((self.number_atoms,3))
        for Ni in range(self.number_atoms-1):
            position_i = self.atoms_positions[Ni]
            for Nj in np.arange(Ni+1,self.number_atoms):
                position_j = self.atoms_positions[Nj]
                box_size = np.diff(self.box_boundaries).reshape(3)
                rij_xyz = (np.remainder(position_i - position_j + box_size/2., box_size) - box_size/2.).T
                rij = np.sqrt(np.sum(rij_xyz**2))
                if rij < self.cut_off:
                    dU_dr = 48/rij*(1/rij**12-0.5/rij**6)
                    forces[Ni] += dU_dr*rij_xyz/rij
                    forces[Nj] -= dU_dr*rij_xyz/rij
        return forces

    def integrate_equation_of_motion(self):
        """Integrate equation of motion using half-step velocity"""
        if self.step == 1:
            self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
        atoms_velocity_Dt2 = self.atoms_velocities + self.atoms_accelerations*self.time_step/2
        self.atoms_positions = self.atoms_positions + atoms_velocity_Dt2*self.time_step
        self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
        self.atoms_velocities = atoms_velocity_Dt2 + self.atoms_accelerations*self.time_step/2
        if self.tau_temp is not None:
            self.apply_berendsen_thermostat()
        if self.tau_press is not None:
            self.apply_berendsen_barostat()

    def apply_berendsen_thermostat(self):
        """Rescale velocities based on Berendsten thermostat."""
        self.calculate_temperature()
        scale = np.sqrt(1+self.time_step*((self.desired_temperature/self.temperature)-1)/self.tau_temp)
        self.atoms_velocities *= scale

    def apply_berendsen_barostat(self):
        """Rescale box size based on Berendsten barostat."""
        self.calculate_pressure()
        scale = np.sqrt(1+self.time_step*((self.pressure/self.desired_pressure)-1)/self.tau_press)
        self.box_boundaries *= scale
        self.atoms_positions *= scale

    def calculate_temperature(self):
        """Follow the expression given in the LAMMPS documentation"""
        self.calculate_kinetic_energy()
        Ndof = self.dimensions*self.number_atoms-self.dimensions
        self.temperature = 2*self.Ekin/Ndof

    def calculate_kinetic_energy(self):
        self.Ekin = np.sum(self.atom_mass*np.sum(self.atoms_velocities**2, axis=1)/2)

    def calculate_pressure(self):
        """Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smith 2002)"""
        Ndof = self.dimensions*self.number_atoms-self.dimensions    
        volume = np.prod(np.diff(self.box_boundaries))
        p_ideal = (Ndof/self.dimensions)*self.temperature/volume
        p_non_ideal = 1/(volume*self.dimensions)*np.sum(self.atoms_positions*self.evaluate_LJ_force())
        self.pressure = (p_ideal+p_non_ideal)



















class MolecularSimulation:
    def __init__(self,
                 number_atoms,
                 Lx,
                 maximum_steps,
                 dimensions=3,
                 Ly = None,
                 Lz = None,
                 desired_temperature=300,
                 desired_pressure=None,
                 time_step=0.1,
                 epsilon=0.1,
                 sigma=1,
                 atom_mass = 1,
                 displace=None,
                 monte_carlo=False,
                 molecular_dynamics=False,
                 seed=None,
                 thermo=10,
                 dump=10,
                 mu = None,
                 tau_temp = 100,
                 tau_press = 1000,
                 LJ_cut_off = 10,
                 debug=False,
                 ):
        
        self.number_atoms = number_atoms
        self.maximum_steps = maximum_steps
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = dimensions
        self.time_step = time_step
        self.desired_temperature = desired_temperature
        self.desired_pressure= desired_pressure
        self.temperature = self.desired_temperature
        self.epsilon = epsilon # kcal/mol
        self.sigma = sigma # Angstrom
        self.atom_mass = atom_mass
        self.displace = displace
        self.monte_carlo = monte_carlo
        self.molecular_dynamics = molecular_dynamics
        self.thermo = thermo
        self.dump = dump
        self.mu = mu
        self.tau_temp = tau_temp
        self.tau_press = tau_press
        self.LJ_cut_off = LJ_cut_off
        self.debug = debug
        if seed is None:
            self.seed = np.random.randint(10000)
        else:
            self.seed = seed
        np.random.seed(self.seed)
        self.step = 0

        self.atoms_positions = np.zeros((self.number_atoms, self.dimensions))
        self.atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
        self.atoms_accelerations = np.zeros((self.number_atoms, self.dimensions))

        if (self.monte_carlo is None) & (self.mu is None) & (self.displace is None):
            print("Warning: wrong Monte Carlo parameters. Pick either mu or displace")

    def non_dimensionalize(self):
        """Non-dimensionalize all input value"""

        self.d0 = self.sigma # todo : deal with different sigma
        self.m0 = self.atom_mass # todo : deal with different mass
        self.e0 = self.epsilon # todo : deal with different epsilon
        self.m0_kg = self.m0/cst.kilo/cst.Avogadro
        d0_m = self.d0*cst.angstrom
        e0_J = self.e0*cst.calorie*cst.kilo/cst.Avogadro
        t0_s = np.sqrt(self.m0_kg*d0_m**2/e0_J)
        self.t0 = t0_s/cst.femto
        self.kB_kCal_mol_K = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo

        self.Lx /= self.d0
        if self.Ly is not None:
            self.Ly /= self.d0
        if self.Lz is not None:
            self.Lz /= self.d0
        self.time_step /= self.t0
        self.desired_temperature /= self.e0/self.kB_kCal_mol_K
        self.epsilon /= self.e0
        self.sigma /= self.d0
        self.atom_mass /= self.m0
        if self.displace_mc is not None:
            self.displace_mc /= self.d0
        if self.mu is not None:
            self.mu /= self.e0
        self.tau_temp /= self.t0
        self.tau_press /= self.t0
        self.LJ_cut_off /= self.d0

        if self.desired_pressure is not None:
            self.desired_pressure *= cst.atm*cst.angstrom**3*cst.Avogadro/cst.calorie/cst.kilo/self.e0*self.d0**3
    
        self.beta =  1/self.desired_temperature

    def run(self):
        self.non_dimensionalize()
        self.initialize_box()
        self.initialize_positions()
        self.calculate_pressure()
        self.initialize_log()
        self.update_outputs()
        for step in range(1, self.maximum_steps):
            self.step = step
            if self.monte_carlo_move:
                self.monte_carlo_displacement()
            if self.monte_carlo_insert:
                self.monte_carlo_insert_delete()
            if self.molecular_dynamics:
                self.molecular_dynamics_displacement()
            self.wrap_in_box()
            self.update_outputs()
        self.close_log()

    def update_outputs(self):
        self.update_log()   
        self.update_dump()
        self.update_data_files()

    def update_dump(self, filename="dump.lammpstrj"):
        if self.step % self.dump == 0:
            if self.step==0:
                f = open(filename, "w")
            else:
                f = open(filename, "a")
            f.write("ITEM: TIMESTEP\n")
            f.write(str(self.step) + "\n")
            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(str(self.number_atoms) + "\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            for dim in np.arange(self.dimensions):
                f.write(str(self.box_boundaries[dim][0]*self.d0) + " " + str(self.box_boundaries[dim][1]*self.d0) + "\n")
            f.write("ITEM: ATOMS id type x y z vx vy vz\n")
            cpt = 1
            atoms_positions = copy.deepcopy(self.atoms_positions)
            if self.molecular_dynamics:
                atoms_velocities = copy.deepcopy(self.atoms_velocities)
            else:
                atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
            for xyz, vxyz in zip(atoms_positions, atoms_velocities):
                #f.write(str(cpt+1) + " 1 " +str(xyz[0])+" "+str(xyz[1])+" "+str(xyz[2])+"\n") 
                f.write(str(cpt) + " " + str(1) + " " +str(xyz[0]*self.d0)+" "+str(xyz[1]*self.d0)+" "+str(xyz[2]*self.d0) + " " +str(vxyz[0]*self.d0/self.t0)+" "+str(vxyz[1]*self.d0/self.t0)+" "+str(vxyz[2]*self.d0/self.t0)+"\n") 
                cpt += 1
            f.close()

    def update_data_files(self):
        if self.step % self.dump == 0:
            if self.step==0:
                epot = open("Epot.dat", "w")
                ekin = open("Ekin.dat", "w")
                etot = open("Etot.dat", "w")
                press = open("pressure.dat", "w")
                temperature = open("temperature.dat", "w")
                density = open("density.dat", "w")
                volume = open("volume.dat", "w")
            else:
                epot = open("Epot.dat", "a")
                ekin = open("Ekin.dat", "a")
                etot = open("Etot.dat", "a")
                press = open("pressure.dat", "a")
                temperature = open("temperature.dat", "a")
                density = open("density.dat", "a")
                volume = open("volume.dat", "a")
            epot.write(str(self.step) + " " + str(self.Epot*self.e0) + "\n")
            ekin.write(str(self.step) + " " + str(self.Ekin*self.e0) + "\n")
            etot.write(str(self.step) + " " + str(self.Ekin*self.e0+self.Epot*self.e0) + "\n")
            self.calculate_pressure()
            pressure = self.pressure*self.e0/self.d0**3 # kcal/mol/A3
            press_atm = pressure*cst.calorie*cst.kilo/cst.Avogadro/cst.angstrom**3/cst.atm # atm
            press.write(str(self.step) + " " + str(press_atm) + "\n")
            T_K = self.temperature*(self.e0/self.kB_kCal_mol_K)
            temperature.write(str(self.step) + " " + str(T_K) + "\n")
            density.write(str(self.step) + " " + str(self.number_atoms/self.volume) + "\n")
            volume.write(str(self.step) + " " + str(self.volume) + "\n")
            epot.close()
            ekin.close()
            etot.close()
            press.close()
            density.close()
            volume.close()

    def initialize_log(self):
        
        # Create and configure logger 
        logging.basicConfig(filename="log-file.txt", 
                            format='%(message)s', 
                            filemode='w') 

        # Create a logger object 
        self.logger=logging.getLogger() 
        self.logger.setLevel(logging.DEBUG)
        self.logger.info("Initialization of the system") 
        self.logger.info("----------------------------\n")
        self.logger.info("initial number of atoms = " + str(self.number_atoms)) 
        if self.molecular_dynamics:
            self.logger.info("initial velocity of CoM = " + str(np.round(self.velocity_com,3).tolist()) + " (NO UNIT SO FAR)") 
            self.logger.info("initial kinetic energy = " + str(np.round(self.Ekin,3)) + " (NO UNIT SO FAR)") 
        self.logger.info("initial dimensions = " + str(self.dimensions)) 
        self.logger.info("\n")
        self.logger.info("Box boundaries")
        self.logger.info("xlow = " + str(self.box_boundaries[0][0]) + " Å, xhigh = " + str(self.box_boundaries[0][1]) + " Å") 
        if self.dimensions>1:
            self.logger.info("ylow = " + str(self.box_boundaries[1][0]) + " Å, yhigh = " + str(self.box_boundaries[1][1]) + " Å") 
        if self.dimensions>2:
            self.logger.info("zlow = " + str(self.box_boundaries[2][0]) + " Å, zhigh = " + str(self.box_boundaries[2][1]) + " Å") 
        self.logger.info("\n") 
        self.logger.info("Interaction parameters") 
        self.logger.info("epsilon = " + str(self.epsilon) + " Kcal/mol") 
        self.logger.info("sigma = " + str(self.sigma) + " Å")
        self.logger.info("atom mass = " + str(self.atom_mass) + " Å\n")
        self.logger.info("Monte Carlo parameters") 
        self.logger.info("displace_mc = " + str(self.displace_mc) + " Å\n")

    def update_log(self):
        if self.step == 0:
            self.logger.info("---------------------------------")
            self.logger.info("Starting the molecular simulation")
            self.logger.info("---------------------------------\n")
            self.logger.info("step n-atoms Epot Ekin press T")
        if self.step % self.thermo == 0:
            self.Epot = self.calculate_potential_energy(self.atoms_positions)
            Epot = np.round(self.Epot,2)*self.e0
            if self.molecular_dynamics:
                Ekin = np.round(self.Ekin,2)*self.e0
            else:
                Ekin = 0
            press = np.round(self.pressure,2)*self.e0/self.d0**3 # kcal/mol/A3
            press_MPa = press*cst.calorie*cst.kilo/cst.Avogadro/cst.angstrom**3/cst.mega # MPa
            T_K = self.temperature*(self.e0/self.kB_kCal_mol_K)
            self.logger.info(str(self.step) + " " + str(self.number_atoms) + " " + str(Epot) + " " + str(Ekin)+ " " + str(press_MPa) + " " + str(T_K))

    def close_log(self):
        handlers = self.logger.handlers[:]
        for handler in handlers:
            handler.close()
            self.logger.removeHandler(handler)

    def initialize_box(self):
        """Define box boundaries based on Lx, Ly, and Lz.

        If Ly or Lz are None, then Lx is used instead"""
        box_boundaries = np.zeros((self.dimensions, 2))
        box_size = np.zeros(self.dimensions)
        box_boundaries[0] = -self.Lx/2, self.Lx/2
        box_size[0] = np.diff(box_boundaries[0])
        if self.dimensions > 1:
            if self.Ly is None:
                box_boundaries[1] = -self.Lx/2, self.Lx/2
            else:
                box_boundaries[1] = -self.Ly/2, self.Ly/2
            box_size[1] = np.diff(box_boundaries[1])
        if self.dimensions > 2:
            if self.Lz is None:
                box_boundaries[2] = -self.Lx/2, self.Lx/2
            else:
                box_boundaries[2] = -self.Lz/2, self.Lz/2
            box_size[2] = np.diff(box_boundaries[2])
        self.box_boundaries = box_boundaries
        self.box_size = box_size
        self.volume = np.prod(box_size)

    def initialize_positions(self):
        """Randomly pick positions and velocities."""

        atoms_positions = np.zeros((self.number_atoms, self.dimensions))
        for dim in np.arange(self.dimensions):
            atoms_positions[:, dim] = np.random.random(self.number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2
        
        if self.debug:
            atoms_positions[0] = np.array([0, 0, 1])
            atoms_positions[1] = np.array([0, 0, 5])
        
        self.atoms_positions = atoms_positions

        if self.molecular_dynamics:
            atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
            for dim in np.arange(self.dimensions):  
                atoms_velocities[:, dim] = np.random.normal(size=self.number_atoms)

            if self.debug:
                atoms_velocities[0] = np.array([0, 0, 0.5])
                atoms_velocities[1] = np.array([0, 0, -0.5])
            self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
            self.atoms_velocities = atoms_velocities * np.sqrt(self.desired_temperature/self.atom_mass/self.dimensions) # doto verif
        
        #self.previous_positions = self.atoms_positions - self.atoms_velocities*self.time_step
        self.calculate_velocity_com()
        self.calculate_kinetic_energy()
        self.calculate_temperature()
    
        # to rescale temperature (attention units)
        # self.atoms_velocities = self.atoms_velocities*(self.desired_temperature/self.temperature)**0.5
        # self.calculate_temperature()

    def calculate_temperature(self):
        # from LAMMPS documentation
        self.calculate_kinetic_energy()
        Ndof = self.dimensions*self.number_atoms-self.dimensions
        self.temperature = 2*self.Ekin/Ndof

    def calculate_velocity_com(self):
        self.velocity_com = np.sum(self.atoms_velocities, axis=0)/self.number_atoms
    
    def calculate_kinetic_energy(self):
        self.Ekin = np.sum(self.atom_mass*np.sum(self.atoms_velocities**2, axis=1)/2)

    def calculate_r(self, position_i, positions_j, number_atoms = None):
        """Calculate the shortest distance between position_i and atoms_positions."""
        if number_atoms is None:
            rij2 = np.zeros(self.number_atoms)
        else:
            rij2 = np.zeros(number_atoms)
        rij = (np.remainder(position_i - positions_j + self.box_size/2., self.box_size) - self.box_size/2.).T
        for dim in np.arange(self.dimensions):
            rij2 += np.power(rij[dim, :], 2)
        return np.sqrt(rij2)

    def calculate_potential_energy(self, atoms_positions, number_atoms = None):
        """Calculate potential energy assuming Lennard-Jones potential interaction."""
        energy_potential = 0
        for position_i in atoms_positions:
            r = self.calculate_r(position_i, atoms_positions, number_atoms)
            energy_potential_i = np.sum(4*(1/np.power(r[r>0], 12)-1/np.power(r[r>0], 6)))
            energy_potential += energy_potential_i
        energy_potential /= 2 # Avoid counting potential energy twice
        return energy_potential
    
    def monte_carlo_displacement(self):
        Epot = self.calculate_potential_energy(self.atoms_positions)
        trial_atoms_positions = copy.deepcopy(self.atoms_positions)
        atom_id = np.random.randint(self.number_atoms)
        trial_atoms_positions[atom_id] += (np.random.random(3)-0.5)*self.displace_mc
        trial_Epot = self.calculate_potential_energy(trial_atoms_positions)
        acceptation_probability = np.min([1, np.exp(-self.beta*(trial_Epot-Epot))])
        if np.random.random() <= acceptation_probability:
            self.atoms_positions = trial_atoms_positions
            self.Epot = trial_Epot
        else:
            self.Epot = Epot 

    def molecular_dynamics_displacement(self):
        atoms_velocity_Dt2 = self.atoms_velocities + self.atoms_accelerations*self.time_step/2
        self.atoms_positions = self.atoms_positions + atoms_velocity_Dt2*self.time_step
        self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
        self.atoms_velocities = atoms_velocity_Dt2 + self.atoms_accelerations*self.time_step/2
        self.apply_berendsen_thermostat()
        if self.desired_pressure is not None:
            self.apply_berendsen_barostat()

    def monte_carlo_insert_delete(self):
        Epot = self.calculate_potential_energy(self.atoms_positions)
        trial_atoms_positions = copy.deepcopy(self.atoms_positions)
        if np.random.random() < 0.5:
            number_atoms = self.number_atoms + 1
            atom_position = np.zeros((1, self.dimensions))
            for dim in np.arange(self.dimensions):
                atom_position[:, dim] = np.random.random(1)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2
            trial_atoms_positions = np.vstack([trial_atoms_positions, atom_position])
            trial_Epot = self.calculate_potential_energy(trial_atoms_positions, number_atoms = number_atoms)
            Lambda = self.calculate_Lambda(self.atom_mass)
            acceptation_probability = np.min([1, self.volume/(Lambda**self.dimensions*(self.number_atoms + 1))*np.exp(self.beta*(self.mu-trial_Epot+Epot))])
        else:
            number_atoms = self.number_atoms - 1
            if number_atoms > 0:
                atom_id = np.random.randint(self.number_atoms)
                trial_atoms_positions = np.delete(trial_atoms_positions, atom_id, axis=0)
                trial_Epot = self.calculate_potential_energy(trial_atoms_positions, number_atoms = number_atoms)
                Lambda = self.calculate_Lambda(self.atom_mass)
                acceptation_probability = np.min([1, (Lambda**self.dimensions*(self.number_atoms)/self.volume)*np.exp(-self.beta*(self.mu+trial_Epot-Epot))])
            else:
                acceptation_probability = 0
        if np.random.random() < acceptation_probability:
            self.atoms_positions = trial_atoms_positions
            self.Epot = trial_Epot
            self.number_atoms = number_atoms
        else:
            self.Epot = Epot

    def calculate_Lambda(self, mass):
        m_kg = mass/cst.Avogadro*cst.milli*self.m0 # kg
        T_K = self.desired_temperature*self.e0/self.kB_kCal_mol_K
        Lambda = cst.h/np.sqrt(2*np.pi*cst.Boltzmann*m_kg*T_K)/cst.angstrom # de Broglie wavelenght, Angstrom
        return Lambda/self.d0
    
    def wrap_in_box(self):
        for dim in np.arange(self.dimensions):
            out_ids = self.atoms_positions[:, dim] > self.box_boundaries[dim][1]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] -= self.box_size[dim]
            out_ids = self.atoms_positions[:, dim] < self.box_boundaries[dim][0]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] += self.box_size[dim]

    def calculate_pressure(self):
        "Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smith 2002)"
        Ndof = self.dimensions*self.number_atoms-self.dimensions
        N = Ndof / self.dimensions            
        p_ideal = N*self.temperature/self.volume
        p_non_ideal = 1/(self.volume*self.dimensions)*np.sum(self.atoms_positions*self.evaluate_LJ_force())
        self.pressure = (p_ideal+p_non_ideal)

    def evaluate_LJ_force(self):
        """Evaluate force based on LJ potential derivative."""
        forces = np.zeros((self.number_atoms,3))
        for Ni in range(self.number_atoms-1):
            position_i = self.atoms_positions[Ni]
            for Nj in np.arange(Ni+1,self.number_atoms):
                position_j = self.atoms_positions[Nj]
                rij_xyz = (np.remainder(position_i - position_j + self.box_size/2., self.box_size) - self.box_size/2.).T
                rij = np.sqrt(np.sum(rij_xyz**2))
                if rij < self.LJ_cut_off:
                    dU_dr = 48/rij*(1/rij**12-0.5/rij**6)
                    forces[Ni] += dU_dr*rij_xyz/rij
                    forces[Nj] -= dU_dr*rij_xyz/rij
        return forces

    def apply_berendsen_thermostat(self):
        """Rescale velocities based on Berendsten thermostat"""
        self.calculate_temperature()
        scale = np.sqrt(1+self.time_step*((self.desired_temperature/self.temperature)-1)/self.tau_temp)
        self.atoms_velocities *= scale
        self.calculate_temperature()

    def apply_berendsen_barostat(self):
        """Rescale box size based on Berendsten barostat"""
        self.calculate_pressure()
        scale = np.sqrt(1+self.time_step*((self.pressure/self.desired_pressure)-1)/self.tau_press)
        self.volume *= scale
        self.box_boundaries *= scale
        self.box_size *= scale
        self.atoms_positions *= scale
        self.calculate_pressure()


    #def _molecular_dynamics_displacement(self):
    #    print("here")
    #    force = self.evaluate_LJ_force()
    #    atoms_positions = 2*self.atoms_positions-self.previous_positions + self.time_step**2*force/self.atom_mass # todo check periodic boundary condition issues        
    #    self.previous_positions = self.atoms_positions
    #    self.atoms_positions = atoms_positions
    #    self.atoms_velocities = (self.atoms_positions - self.previous_positions)/self.time_step
    #    self.apply_berendsen_thermostat()    

    #def __molecular_dynamics_displacement(self):
    #    print("here")
    #    self.atoms_positions = self.atoms_positions + self.atoms_velocities*self.time_step + self.atoms_accelerations*self.time_step**2/2  
    #    atoms_accelerations_Dt = self.evaluate_LJ_force()/self.atom_mass
    #    self.atoms_velocities = self.atoms_velocities + atoms_accelerations_Dt*self.time_step/2 + self.atoms_accelerations*self.time_step/2   
    #    self.atoms_accelerations = atoms_accelerations_Dt
    #    #self.apply_berendsen_thermostat()   