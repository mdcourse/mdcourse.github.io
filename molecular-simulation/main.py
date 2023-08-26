from scipy import constants as cst
from decimal import Decimal
import numpy as np
import logging
import copy

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
                 desired_pressure=1,
                 provided_positions=None,
                 provided_velocities=None,
                 *args,
                 **kwargs,
                 ):
        super().__init__(*args, **kwargs) 
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
        self.desired_pressure = desired_pressure
        self.provided_positions = provided_positions
        self.provided_velocities = provided_velocities

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
        self.desired_pressure = self.nondimensionalise_units(self.desired_pressure, "pressure")

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
            elif type == "pressure":
                variable *= cst.atm*cst.angstrom**3*cst.Avogadro/cst.calorie/cst.kilo/self.reference_energy*self.reference_distance**3
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
        if self.provided_positions is not None:
            atoms_positions = self.provided_positions/self.reference_distance
        else:
            for dim in np.arange(self.dimensions):
                atoms_positions[:, dim] = np.random.random(self.number_atoms)*np.diff(self.box_boundaries[dim]) - np.diff(self.box_boundaries[dim])/2    
        self.atoms_positions = atoms_positions

    def give_velocity(self):
        """Give velocity to atoms so that the initial temperature is the desired one."""
        atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
        if self.provided_velocities is not None:
            atoms_velocities = self.provided_velocities/self.reference_distance*self.reference_time
        else:
            for dim in np.arange(self.dimensions):  
                atoms_velocities[:, dim] = np.random.normal(size=self.number_atoms)
        self.atoms_velocities = atoms_velocities
        self.calculate_temperature()
        scale = np.sqrt(1+((self.desired_temperature/self.temperature)-1))
        self.atoms_velocities *= scale

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
                 thermo = None,
                 dump = None,
                 *args,
                 **kwargs):
        self.thermo = thermo
        self.dump = dump
        super().__init__(*args, **kwargs)

    def evaluate_temperature(self):
        """Measure temperature and convert in Kelvin."""
        # Evaluate temperature (in K)
        self.calculate_temperature()
        kB = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo # kCal/mol/K
        return np.round(self.temperature*self.reference_energy/kB,2)

    def evaluate_pressure(self):
        """Measure pressure and convert in atmosphere."""
        self.calculate_pressure()
        return np.round(self.pressure*self.reference_energy/self.reference_distance**3*cst.calorie*cst.kilo/cst.Avogadro/cst.angstrom**3/cst.atm, 2)

    def evaluate_volume(self):
        """Measure volume and convert in Angstrom 3"""
        return np.round(np.prod(np.diff(self.box_boundaries))*self.reference_distance**3)

    def evaluate_potential_energy(self):
        """Measure energy and convert in kcal/mol"""
        Epot = self.calculate_potential_energy(self.atoms_positions)
        return Epot*self.reference_energy
    
    def evaluate_kinetic_energy(self):
        """Measure energy and convert in kcal/mol"""
        self.calculate_kinetic_energy()
        return self.Ekin*self.reference_energy
    
    def evaluate_density(self):
        return self.number_atoms/np.round(np.prod(np.diff(self.box_boundaries))*self.reference_distance**3)

    def update_log(self):
        if self.thermo is not None:
            if (self.step % self.thermo == 0) | (self.step == 0):
                temperature = self.evaluate_temperature()
                pressure = self.evaluate_pressure()
                volume = self.evaluate_volume()
                epot = self.evaluate_potential_energy()
                ekin = self.evaluate_kinetic_energy()
                density = self.evaluate_density()
                if self.step == 0:
                    print("step N temp epot ekin press vol")
                print(self.step,
                      self.number_atoms,
                      temperature,
                      '%.2E' % Decimal(epot),
                      '%.2E' % Decimal(ekin),
                      '%.2E' % Decimal(pressure),
                      '%.2E' % Decimal(volume),
                      )
                self.write_data_file(temperature, "temperature.dat")
                self.write_data_file(pressure, "pressure.dat")
                self.write_data_file(volume, "volume.dat")
                self.write_data_file(epot, "Epot.dat")
                self.write_data_file(ekin, "Ekin.dat")
                self.write_data_file(density, "density.dat")

    def write_data_file(self, output_value, filename):
        if self.step == 0:
            myfile = open(filename, "w")
        else:
            myfile = open(filename, "a")
        myfile.write(str(self.step) + " " + str(output_value) + "\n")
        myfile.close()

    def update_dump(self, filename="dump.lammpstrj", velocity=True):
        if self.dump is not None:
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
                    f.write(str(self.box_boundaries[dim][0]*self.reference_distance) + " " + str(self.box_boundaries[dim][1]*self.reference_distance) + "\n")
                f.write("ITEM: ATOMS id type x y z vx vy vz\n")
                cpt = 1
                atoms_positions = copy.deepcopy(self.atoms_positions)
                if velocity:
                    atoms_velocities = copy.deepcopy(self.atoms_velocities)
                else:
                    atoms_velocities = np.zeros((self.number_atoms, self.dimensions))
                for xyz, vxyz in zip(atoms_positions, atoms_velocities):
                    f.write(str(cpt) + " " + str(1) + " " +str(xyz[0]*self.reference_distance)+" "+str(xyz[1]*self.reference_distance)+" "+str(xyz[2]*self.reference_distance) + " " +str(vxyz[0]*self.reference_distance/self.reference_time)+" "+str(vxyz[1]*self.reference_distance/self.reference_time)+" "+str(vxyz[2]*self.reference_distance/self.reference_time)+"\n") 
                    cpt += 1
                f.close()


class Utilities:
    def __init__(self,
                 *args,
                 **kwargs):
        super().__init__(*args, **kwargs)

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
        self.calculate_temperature()
        p_ideal = (Ndof/self.dimensions)*self.temperature/volume
        p_non_ideal = 1/(volume*self.dimensions)*np.sum(self.evaluate_LJ_matrix()*self.evaluate_rij_matrix())
        self.pressure = (p_ideal+p_non_ideal)

    def evaluate_rij_matrix(self):
        """Evaluate vector rij between particles."""
        rij = np.zeros((self.number_atoms,self.number_atoms,3))
        for Ni in range(self.number_atoms-1):
            position_i = self.atoms_positions[Ni]
            for Nj in np.arange(Ni+1,self.number_atoms):
                position_j = self.atoms_positions[Nj]
                box_size = np.diff(self.box_boundaries).reshape(3)
                rij_xyz = (np.remainder(position_i - position_j + box_size/2., box_size) - box_size/2.).T
                r = np.sqrt(np.sum(rij_xyz**2))
                if r < self.cut_off:
                    rij[Ni][Nj] = rij_xyz
        return rij

    def calculate_r(self, position_i, positions_j, number_atoms = None):
        """Calculate the shortest distance between position_i and atoms_positions."""
        if number_atoms is None:
            rij2 = np.zeros(self.number_atoms)
        else:
            rij2 = np.zeros(number_atoms)
        box_size = np.diff(self.box_boundaries).reshape(3)
        rij = (np.remainder(position_i - positions_j + box_size/2., box_size) - box_size/2.).T
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
    
    def evaluate_LJ_matrix(self):
        """Evaluate force based on LJ potential derivative."""
        forces = np.zeros((self.number_atoms,self.number_atoms,3))
        for Ni in range(self.number_atoms-1):
            position_i = self.atoms_positions[Ni]
            for Nj in np.arange(Ni+1,self.number_atoms):
                position_j = self.atoms_positions[Nj]
                box_size = np.diff(self.box_boundaries).reshape(3)
                rij_xyz = (np.remainder(position_i - position_j + box_size/2., box_size) - box_size/2.).T
                rij = np.sqrt(np.sum(rij_xyz**2))
                if rij < self.cut_off:
                    dU_dr = 48/rij*(1/rij**12-0.5/rij**6)
                    forces[Ni][Nj] += dU_dr*rij_xyz/rij
        return forces

class MolecularDynamics(InitializeSimulation, Utilities, Outputs):
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
        self.tau_temp = self.nondimensionalise_units(self.tau_temp, "time")
        self.tau_press = self.nondimensionalise_units(self.tau_press, "time")

    def run(self):
        """Perform the loop over time."""
        for self.step in range(0, self.maximum_steps+1):
            self.integrate_equation_of_motion()
            self.wrap_in_box()
            if self.tau_temp is not None:
                self.apply_berendsen_thermostat()
            if self.tau_press is not None:
                self.apply_berendsen_barostat()
            self.update_log()
            self.update_dump()
        self.write_lammps_data(filename="final.data")

    def integrate_equation_of_motion(self):
        """Integrate equation of motion using half-step velocity"""
        if self.step == 0:
            self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
        atoms_velocity_Dt2 = self.atoms_velocities + self.atoms_accelerations*self.time_step/2
        self.atoms_positions = self.atoms_positions + atoms_velocity_Dt2*self.time_step
        self.atoms_accelerations = self.evaluate_LJ_force()/self.atom_mass
        self.atoms_velocities = atoms_velocity_Dt2 + self.atoms_accelerations*self.time_step/2

    def apply_berendsen_thermostat(self):
        """Rescale velocities based on Berendsen thermostat."""
        self.calculate_temperature()
        scale = np.sqrt(1+self.time_step*((self.desired_temperature/self.temperature)-1)/self.tau_temp)
        self.atoms_velocities *= scale

    def apply_berendsen_barostat(self):
        """Rescale box size based on Berendsten barostat."""
        self.calculate_pressure()
        scale = np.sqrt(1+self.time_step*((self.pressure/self.desired_pressure)-1)/self.tau_press)
        self.box_boundaries *= scale
        self.atoms_positions *= scale


class MonteCarlo(InitializeSimulation, Utilities, Outputs):
    def __init__(self,
        maximum_steps,
        cut_off = 10,
        displace_mc=None,
        mu = None,
        *args,
        **kwargs,
        ):

        self.maximum_steps = maximum_steps
        self.cut_off = cut_off
        self.displace_mc = displace_mc
        self.mu = mu
        super().__init__(*args, **kwargs)

        self.cut_off = self.nondimensionalise_units(self.cut_off, "distance")
        self.displace_mc = self.nondimensionalise_units(self.displace_mc, "distance")
        self.mu = self.nondimensionalise_units(self.mu, "energy")

    def run(self):
        """Perform the loop over time."""
        
        for self.step in range(0, self.maximum_steps+1):
            if self.displace_mc is not None:
                self.monte_carlo_displacement()
            if self.mu is not None:
                self.monte_carlo_insert_delete()
            self.wrap_in_box()
            self.update_log()
            self.update_dump(velocity=False)
        self.write_lammps_data(filename="final.data")

    def monte_carlo_displacement(self):
        beta =  1/self.desired_temperature
        Epot = self.calculate_potential_energy(self.atoms_positions)
        trial_atoms_positions = copy.deepcopy(self.atoms_positions)
        atom_id = np.random.randint(self.number_atoms)
        trial_atoms_positions[atom_id] += (np.random.random(3)-0.5)*self.displace_mc
        trial_Epot = self.calculate_potential_energy(trial_atoms_positions)
        acceptation_probability = np.min([1, np.exp(-beta*(trial_Epot-Epot))])
        if np.random.random() <= acceptation_probability:
            self.atoms_positions = trial_atoms_positions
            self.Epot = trial_Epot
        else:
            self.Epot = Epot 

    def calculate_Lambda(self, mass):
        """Estimate de Broglie wavelength in LJ units."""
        m_kg = mass/cst.Avogadro*cst.milli*self.reference_mass
        kB_kCal_mol_K = cst.Boltzmann*cst.Avogadro/cst.calorie/cst.kilo
        T_K = self.desired_temperature*self.reference_energy/kB_kCal_mol_K
        Lambda = cst.h/np.sqrt(2*np.pi*cst.Boltzmann*m_kg*T_K)/cst.angstrom
        return self.nondimensionalise_units(Lambda, "distance")

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
            volume = np.prod(np.diff(self.box_boundaries))
            beta = 1/self.desired_temperature
            acceptation_probability = np.min([1, volume/(Lambda**self.dimensions*(self.number_atoms + 1))*np.exp(beta*(self.mu-trial_Epot+Epot))])
        else:
            number_atoms = self.number_atoms - 1
            if number_atoms > 0:
                atom_id = np.random.randint(self.number_atoms)
                trial_atoms_positions = np.delete(trial_atoms_positions, atom_id, axis=0)
                trial_Epot = self.calculate_potential_energy(trial_atoms_positions, number_atoms = number_atoms)
                Lambda = self.calculate_Lambda(self.atom_mass)
                volume = np.prod(np.diff(self.box_boundaries))
                beta = 1/self.desired_temperature
                acceptation_probability = np.min([1, (Lambda**self.dimensions*(self.number_atoms)/volume)*np.exp(-beta*(self.mu+trial_Epot-Epot))])
            else:
                acceptation_probability = 0
        if np.random.random() < acceptation_probability:
            self.atoms_positions = trial_atoms_positions
            self.Epot = trial_Epot
            self.number_atoms = number_atoms
        else:
            self.Epot = Epot
