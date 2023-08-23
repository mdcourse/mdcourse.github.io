from scipy import constants as cst
import numpy as np
import logging
import copy
import sys
import matplotlib.pyplot as plt

# to remove the runtime warning
import warnings
warnings.filterwarnings('ignore')

class MolecularSimulation:
    def __init__(self,
                 number_atoms,
                 Lx,
                 maximum_steps,
                 dimensions=3,
                 Ly = None,
                 Lz = None,
                 desired_temperature=300,
                 time_step=0.1,
                 epsilon=0.1,
                 sigma=3,
                 atom_mass = 1,
                 displace_mc=0.5,
                 monte_carlo_move=False,
                 monte_carlo_insert=False,
                 molecular_dynamics=False,
                 seed=None,
                 thermo=10,
                 dump=10,
                 mu = -0.2,
                 ):
        
        self.number_atoms = number_atoms
        self.maximum_steps = maximum_steps
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dimensions = dimensions
        self.time_step = time_step
        self.desired_temperature = desired_temperature
        self.epsilon = epsilon # kcal/mol
        self.sigma = sigma # Angstrom
        self.atom_mass = atom_mass
        self.displace_mc = displace_mc
        self.monte_carlo_move = monte_carlo_move
        self.monte_carlo_insert = monte_carlo_insert
        self.molecular_dynamics = molecular_dynamics
        self.thermo = thermo
        self.dump = dump
        self.mu = mu
        if seed is None:
            self.seed = np.random.randint(10000)
        else:
            self.seed = seed
        np.random.seed(self.seed)

        self.beta = 1/(cst.Boltzmann*self.desired_temperature/cst.calorie/cst.kilo*cst.Avogadro)  # mol/kCal

    def run(self):
        self.initialize_box()
        self.initialize_positions()
        self.print_log()
        for step in range(self.maximum_steps):
            self.step = step
            if self.monte_carlo_move:
                self.monte_carlo_displacement()
            if self.monte_carlo_insert:
                self.monte_carlo_insert_delete()
            if self.molecular_dynamics:
                self.molecular_dynamics_displacement()
            self.wrap_in_box()
            self.update_log()
            self.update_data_files()
            self.update_dump()
        self.close_log()

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
                f.write(str(self.box_boundaries[dim][0]) + " " + str(self.box_boundaries[dim][1]) + "\n")
            f.write("ITEM: ATOMS id type x y z\n")
            for cpt, xyz in enumerate(self.atoms_positions):
                #f.write(str(cpt+1) + " 1 " +str(xyz[0])+" "+str(xyz[1])+" "+str(xyz[2])+"\n") 
                f.write(str(cpt+1) + " " + str(cpt+1) + " " +str(xyz[0])+" "+str(xyz[1])+" "+str(xyz[2])+"\n") 
            f.close()

    def update_data_files(self):
        if self.step % self.dump == 0:
            if self.step==0:
                epot = open("Epot.dat", "w")
                ekin = open("Ekin.dat", "w")
                etot = open("Etot.dat", "w")
            else:
                epot = open("Epot.dat", "a")
                ekin = open("Ekin.dat", "a")
                etot = open("Etot.dat", "a")
            epot.write(str(self.step) + " " + str(self.Epot) + "\n")
            ekin.write(str(self.step) + " " + str(self.Ekin) + "\n")
            etot.write(str(self.step) + " " + str(self.Ekin+self.Epot) + "\n")
            epot.close()
            ekin.close()
            etot.close()

    def print_log(self):
        
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
            self.logger.info("step n-atoms Epot Ekin press")
        if self.step % self.thermo == 0:
            self.calculate_pressure()
            if (self.monte_carlo_move is False) & (self.monte_carlo_insert is False):
                self.Epot = self.calculate_potential_energy(self.atoms_positions)
            if self.molecular_dynamics:
                self.calculate_kinetic_energy()
            else:
                self.Ekin = 0
            self.logger.info(str(self.step) + " " + str(self.number_atoms) + " " + str(np.round(self.Epot,2)) + " " + str(np.round(self.Ekin,2))+ " " + str(np.round(self.pressure,2)))

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
        
        #atoms_positions[0] = np.array([0, 0, 1])
        #atoms_positions[1] = np.array([0, 0, 5])
        
        self.atoms_positions = atoms_positions
        if self.molecular_dynamics:
            atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
            for dim in np.arange(self.dimensions):  
                atoms_velocities[:, dim] = np.random.random(self.number_atoms)-0.5 # todo - must be rescaled

            print(np.max(atoms_velocities))

            #atoms_velocities[0] = np.array([0, 0, 0.02])
            #atoms_velocities[1] = np.array([0, 0, -0.02])

            self.atoms_velocities = atoms_velocities
            self.previous_positions = self.atoms_positions - self.atoms_velocities*self.time_step
            self.calculate_velocity_com()
            self.calculate_kinetic_energy()
            self.calculate_temperature()
            
            # to rescale temperature (attention units)
            # self.atoms_velocities = self.atoms_velocities*(self.desired_temperature/self.temperature)**0.5
            # self.calculate_temperature()

    def calculate_temperature(self):
        """Equation 4.2.2 in Frenkel-Smith 2002"""
        self.temperature = np.sum(self.atom_mass*np.linalg.norm(self.atoms_velocities, axis=1)**2)/self.number_atoms/cst.Boltzmann

    def calculate_velocity_com(self):
        self.velocity_com = np.sum(self.atoms_velocities, axis=0)/self.number_atoms
    
    def calculate_kinetic_energy(self):
        kinetic_energy = 0
        for dim in np.arange(self.dimensions): 
            kinetic_energy += self.atom_mass*np.sum(self.atoms_velocities[:, dim]**2)
        # convert g/mol * A2/fs2 into kcal/mol
        kinetic_energy *= cst.angstrom**2/cst.femto**2/cst.kilo/cst.calorie/cst.kilo
        self.Ekin = kinetic_energy/self.number_atoms

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
        """Calculate potential energy assuming Lennard-Jones potential interaction.
        
        Output units are kcal/mol."""
        energy_potential = 0
        for position_i in atoms_positions:
            r = self.calculate_r(position_i, atoms_positions, number_atoms)
            energy_potential_i = np.sum(4*self.epsilon*(np.power(self.sigma/r[r>0], 12)-np.power(self.sigma/r[r>0], 6)))
            energy_potential += energy_potential_i
        # Avoid counting potential energy twice
        energy_potential /= 2
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

        force_kcalmolA = self.evaluate_LJ_force()
        force_SI = force_kcalmolA*cst.calorie*cst.kilo/cst.angstrom/cst.Avogadro
        mass_SI = self.atom_mass/cst.kilo/cst.Avogadro # todo : allow different masses
        acc_SI = force_SI/mass_SI
        acc_Afs2 = acc_SI/cst.angstrom*cst.femto**2
        atoms_positions = 2*self.atoms_positions-self.previous_positions + self.time_step**2*acc_Afs2 # todo check periodic boundary condition issues        
        self.previous_positions = self.atoms_positions
        self.atoms_positions = atoms_positions
        self.atoms_velocities = (self.atoms_positions - self.previous_positions)/2/self.time_step
            
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
        m_kg = mass/cst.Avogadro*cst.milli # kg
        Lambda = cst.h/np.sqrt(2*np.pi*cst.Boltzmann*m_kg*self.desired_temperature)/cst.angstrom # de Broglie wavelenght, Angstrom
        return Lambda
    
    def wrap_in_box(self):
        for dim in np.arange(self.dimensions):
            out_ids = self.atoms_positions[:, dim] > self.box_boundaries[dim][1]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] -= self.box_size[dim]
                if self.molecular_dynamics:
                    self.previous_positions[:, dim][out_ids] -= self.box_size[dim]
            out_ids = self.atoms_positions[:, dim] < self.box_boundaries[dim][0]
            if np.sum(out_ids) > 0:
                self.atoms_positions[:, dim][out_ids] += self.box_size[dim]
                if self.molecular_dynamics:
                    self.previous_positions[:, dim][out_ids] += self.box_size[dim]

    def calculate_pressure(self):
        "Evaluate p based on the Virial equation (Eq. 4.4.2 in Frenkel-Smith 2002)"
        V_m = self.volume*cst.angstrom**3
        p_ideal = self.number_atoms*cst.Boltzmann*self.desired_temperature/V_m
        p_non_ideal = 1/(V_m*self.dimensions)*np.sum(self.atoms_positions*cst.angstrom*self.evaluate_LJ_force()/cst.calorie/cst.kilo/cst.Avogadro/cst.angstrom)
        self.pressure = (p_ideal+p_non_ideal)/cst.mega

    def evaluate_LJ_force(self):
        """Evaluate force based on LJ potential derivative.
        
        Output has unit of kcal/mol/Angstrom"""
        forces = np.zeros((self.number_atoms,3))
        for Ni in range(self.number_atoms-1):
            position_i = self.atoms_positions[Ni]
            for Nj in np.arange(Ni+1,self.number_atoms):
                position_j = self.atoms_positions[Nj]
                rij_xyz = (np.remainder(position_i - position_j + self.box_size/2., self.box_size) - self.box_size/2.).T
                rij = np.sqrt(np.sum(rij_xyz**2))
                dU_dr = 48*self.epsilon/rij*((self.sigma/rij)**12-0.5*(self.sigma/rij)**6)
                forces[Ni] += dU_dr*rij_xyz/rij
                forces[Nj] -= dU_dr*rij_xyz/rij
        return forces
