import numpy as np
import logging
import copy
import sys

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
                 time_step=1,
                 epsilon=0.1,
                 sigma=3,
                 atom_mass = 1,
                 displace_mc=0.5,
                 monte_carlo=False,
                 molecular_dynamics=False,
                 seed=None,
                 thermo=10,
                 dump=10,
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
        self.monte_carlo = monte_carlo
        self.molecular_dynamics = molecular_dynamics
        self.thermo = thermo
        self.dump = dump
        if seed is None:
            self.seed = np.random.randint(10000)
        else:
            self.seed = seed
        np.random.seed(self.seed)

        self.kB = 1.38e-23 # J/K
        self.Na = 6.022e23 # mol-1
        self.J_to_Kcal = 0.0002389
        self.beta = 1/(self.kB*self.desired_temperature*self.J_to_Kcal*self.Na)  # mol/kCal

    def run(self):

        self.initialize_box()
        self.initialize_positions()
        self.print_log()
        for step in range(self.maximum_steps):
            self.step = step
            if self.monte_carlo:
                self.monte_carlo_displacement()
            self.update_log()
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
                f.write(str(cpt+1) + " 1 " +str(xyz[0])+" "+str(xyz[1])+" "+str(xyz[2])+"\n") 
            f.close()

    def print_log(self):
        
        # Create and configure logger 
        logging.basicConfig(filename="log.txt", 
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
            self.logger.info("step Epot")
        if self.step % self.thermo == 0:
            self.logger.info(str(self.step) + " " + str(np.round(self.Epot,2)))

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

    def initialize_positions(self):
        """Randomly pick positions and velocities."""

        atoms_positions = np.zeros((self.number_atoms, self.dimensions))
        for dim in np.arange(self.dimensions):
            atoms_positions[:, dim] = np.random.random(self.number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2
        self.atoms_positions = atoms_positions
        if self.molecular_dynamics:
            atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
            for dim in np.arange(self.dimensions):  
                atoms_velocities[:, dim] = np.random.random(self.number_atoms)-0.5 # todo - must be rescalled
            self.atoms_velocities = atoms_velocities
            self.calculate_velocity_com()
            self.calculate_kinetic_energy()
            self.calculate_temperature()
            
            # to rescale temperature (attention units)
            # self.atoms_velocities = self.atoms_velocities*(self.desired_temperature/self.temperature)**0.5
            # self.calculate_temperature()

    def calculate_temperature(self):
        """Equation 4.2.2 in Frenkel-Smith 2002"""
        self.temperature = np.sum(self.atom_mass*np.linalg.norm(self.atoms_velocities, axis=1)**2)/self.number_atoms/self.kB

    def calculate_velocity_com(self):
        self.velocity_com = np.sum(self.atoms_velocities, axis=0)/self.number_atoms
    
    def calculate_kinetic_energy(self):
        kinetic_energy = 0
        for dim in np.arange(self.dimensions): 
            kinetic_energy += np.sum(self.atoms_velocities[:, dim]**2)
        self.Ekin = kinetic_energy/self.number_atoms

    def calculate_r(self, position_i, positions_j):
        """Calculate the shortest distance between position_i and atoms_positions."""
        rij2 = np.zeros(self.number_atoms)
        rij = (np.remainder(position_i - positions_j + self.box_size/2., self.box_size) - self.box_size/2.).T
        for dim in np.arange(self.dimensions):
            rij2 += np.power(rij[dim, :], 2)
        return np.sqrt(rij2)

    def calculate_potential_energy(self, atoms_positions):
        """Calculate potential energy assuming Lennard-Jones potential interaction.
        
        Output units are kcal/mol."""
        energy_potential = 0
        for position_i in atoms_positions:
            r = self.calculate_r(position_i, atoms_positions)
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
            