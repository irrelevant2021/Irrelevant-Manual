from sys import stdout, exit, stderr
import os, math, fnmatch
import sys

import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
import pytraj as pt

#!pip install --upgrade MDAnalysis 2>&1 1>/dev/null
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms

Google_Drive_Path = '/media/zhennan/78C2BC59C2BC1CF6/openmm-cookbook-main/notebooks/tutorials/makeitrain/8WD4-20250611-pdb_ligand' #@param {type:"string"}
workDir = Google_Drive_Path




#@title ### **Parameters for MD Equilibration protocol:**
import os
# remove whitespaces
Jobname = 'prot_lig_equil' #@param {type:"string"}

Ligand_Force_field = "GAFF2" #@param ["GAFF2", "OpenFF 2.0.0 (Sage)"]

if Ligand_Force_field == "OpenFF 2.0.0 (Sage)":
  top = os.path.join(workDir, "SYS_openff.prmtop")
  crd = os.path.join(workDir, "SYS_openff.inpcrd")
  pdb = os.path.join(workDir, "SYS.pdb")
else:
  top = os.path.join(workDir, "SYS_gaff2.prmtop")
  crd = os.path.join(workDir, "SYS_gaff2.crd")
  pdb = os.path.join(workDir, "SYS.pdb")


Minimization_steps = "5000" #@param ["1000", "5000", "10000", "20000", "50000", "100000"]

#@markdown Simulation time (in nanoseconds) and integration time (in femtoseconds): 
Time = "1" #@param {type:"string"}
stride_time_eq = Time
Integration_timestep = "2" #@param ["0.5", "1", "2", "3", "4"]
dt_eq = Integration_timestep

#@markdown Temperature (in Kelvin) and Pressure (in bar)
Temperature = 298 #@param {type:"string"}
temperature_eq = Temperature
Pressure = 1 #@param {type:"string"}
pressure_eq = Pressure

#@markdown Position restraints force constant (in kJ/mol): 
Force_constant = 700 #@param {type:"slider", min:0, max:2000, step:100}

#@markdown Frequency to write the trajectory file (in picoseconds): 

Write_the_trajectory = "100" #@param ["10", "100", "200", "500", "1000"]
write_the_trajectory_eq = Write_the_trajectory
#@markdown Frequency to write the log file (in picoseconds): 

Write_the_log = "100" #@param ["10", "100", "200", "500", "1000"]
write_the_log_eq = Write_the_log


#@markdown ---




#@title **Runs an Equilibration MD simulation (NPT ensemble)**
#@markdown Now, let's equilibrate our system!

###########################################
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
import pytraj as pt

from sys import stdout, exit, stderr
import os, math, fnmatch

#############################################
# Defining MD simulation parameters

jobname = os.path.join(workDir, Jobname)
coordinatefile = crd
pdbfile = pdb
topologyfile = top

time_ps = float(Time)*1000
simulation_time = float(time_ps)*picosecond		# in ps
dt = int(dt_eq)*femtosecond					
temperature = float(temperature_eq)*kelvin
savcrd_freq = int(write_the_trajectory_eq)*picosecond
print_freq  = int(write_the_log_eq)*picosecond

pressure	= float(pressure_eq)*bar

restraint_fc = int(Force_constant) # kJ/mol

nsteps  = int(simulation_time.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nprint  = int(print_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))
nsavcrd = int(savcrd_freq.value_in_unit(picosecond)/dt.value_in_unit(picosecond))

#############################################
# Defining functions to use below:
def backup_old_log(pattern, string):
	result = []
	for root, dirs, files in os.walk("./"):
		for name in files:
			if fnmatch.fnmatch(name, pattern):

				try:
					number = int(name[-2])
					avail = isinstance(number, int)
					#print(name,avail)
					if avail == True:
						result.append(number)
				except:
					pass

	if len(result) > 0:
		maxnumber = max(result)
	else:
		maxnumber = 0

	backup_file = "\#" + string + "." + str(maxnumber + 1) + "#"
	os.system("mv " + string + " " + backup_file)
	return backup_file

def restraints(system, crd, fc, restraint_array):

	boxlx = system.getDefaultPeriodicBoxVectors()[0][0].value_in_unit(nanometers)
	boxly = system.getDefaultPeriodicBoxVectors()[1][1].value_in_unit(nanometers)
	boxlz = system.getDefaultPeriodicBoxVectors()[2][2].value_in_unit(nanometers)

	if fc > 0:
		# positional restraints for all heavy-atoms
		posresPROT = CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2;')
		posresPROT.addPerParticleParameter('k')
		posresPROT.addPerParticleParameter('x0')
		posresPROT.addPerParticleParameter('y0')
		posresPROT.addPerParticleParameter('z0')
  
		for atom1 in restraint_array:
			atom1 = int(atom1)
               
			xpos  = crd.positions[atom1].value_in_unit(nanometers)[0]
			ypos  = crd.positions[atom1].value_in_unit(nanometers)[1]
			zpos  = crd.positions[atom1].value_in_unit(nanometers)[2]

			posresPROT.addParticle(atom1, [fc, xpos, ypos, zpos])
    
		system.addForce(posresPROT)

	return system
##############################################

#############################################
print("\n> Simulation details:\n")
print("\tJob name = " + jobname)
print("\tCoordinate file = " + str(coordinatefile))
print("\tPDB file = " + str(pdbfile))
print("\tTopology file = " + str(topologyfile))

print("\n\tSimulation_time = " + str(simulation_time))
print("\tIntegration timestep = " + str(dt))
print("\tTotal number of steps = " +  str(nsteps))

print("\n\tSave coordinates each " + str(savcrd_freq))
print("\tPrint in log file each " + str(print_freq))

print("\n\tTemperature = " + str(temperature))
print("\tPressure = " + str(pressure))
#############################################

print("\n> Setting the system:\n")

if Ligand_Force_field == "OpenFF 2.0.0 (Sage)":
  print("\t- Reading topology and structure file...")
  prmtop = pmd.load_file(topologyfile)
  inpcrd = AmberInpcrdFile(coordinatefile)

  print("\t- Creating system and setting parameters...")
  nonbondedMethod = PME
  nonbondedCutoff = 1.0*nanometers
  ewaldErrorTolerance = 0.0005
  constraints = HBonds
  rigidWater = True
  constraintTolerance = 0.000001
  friction = 1.0
  system = complex_structure.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
                                          constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)
else:
  print("\t- Reading topology and structure file...")
  prmtop = AmberPrmtopFile(topologyfile)
  inpcrd = AmberInpcrdFile(coordinatefile)

  print("\t- Creating system and setting parameters...")
  nonbondedMethod = PME
  nonbondedCutoff = 1.0*nanometers
  ewaldErrorTolerance = 0.0005
  constraints = HBonds
  rigidWater = True
  constraintTolerance = 0.000001
  friction = 1.0
  system = prmtop.createSystem(nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
                                          constraints=constraints, rigidWater=rigidWater, ewaldErrorTolerance=ewaldErrorTolerance)


print("\t- Applying restraints. Force Constant = " + str(Force_constant) + "kJ/mol")
pt_system = pt.iterload(coordinatefile, topologyfile)
pt_topology = pt_system.top
#restraint_array = pt.select_atoms('!(:H*) & !(:WAT) & !(:Na+) & !(:Cl-) & !(:Mg+) & !(:K+)', pt_topology)

#system = restraints(system, inpcrd, restraint_fc, restraint_array)

print("\t- Setting barostat...")
system.addForce(MonteCarloBarostat(pressure, temperature))

print("\t- Setting integrator...")
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)




#直接用pytraj选择需要加热的原子，AMBER的拓扑文件是分为两个的，prmtop(拓扑)+inpcrd(坐标)
import copy
rest_atoms = pt.select_atoms(':LIG', pt_topology)
len(rest_atoms)




# Create REST system
rest_system = openmm.System()

# Create dict of vanilla system forces (for easy retrieval of force objects later)
system_forces = {type(force).__name__ : force for force in system.getForces()}

# Add particles
for particle_idx in range(system.getNumParticles()):
    particle_mass = system.getParticleMass(particle_idx)
    rest_system.addParticle(particle_mass)

# Copy barostat
if "MonteCarloBarostat" in system_forces:
    barostat = copy.deepcopy(system_forces["MonteCarloBarostat"])
    rest_system.addForce(barostat)

# Copy box vectors
box_vectors = system.getDefaultPeriodicBoxVectors()
rest_system.setDefaultPeriodicBoxVectors(*box_vectors)

# Copy constraints
for constraint_idx in range(system.getNumConstraints()):
    atom1, atom2, length = system.getConstraintParameters(constraint_idx)
    rest_system.addConstraint(atom1, atom2, length)


# Define the custom expression
bond_expression = "rest_scale * (K / 2) * (r - length)^2;"
bond_expression += "rest_scale = is_rest * lambda_rest_bonds * lambda_rest_bonds " \
                   "+ is_inter * lambda_rest_bonds " \
                   "+ is_nonrest;"

# Create custom force
rest_bond_force = openmm.CustomBondForce(bond_expression)
rest_system.addForce(rest_bond_force)

# Add global parameters
rest_bond_force.addGlobalParameter("lambda_rest_bonds", 1.0)

# Add per-bond parameters for rest scaling
rest_bond_force.addPerBondParameter("is_rest")
rest_bond_force.addPerBondParameter("is_inter")
rest_bond_force.addPerBondParameter("is_nonrest")

# Add per-bond parameters for defining bond energy
rest_bond_force.addPerBondParameter('length')  # equilibrium bond length
rest_bond_force.addPerBondParameter('K')  # force constant


def get_rest_identifier(atoms, rest_atoms):
    """
    For a given atom or set of atoms, get the rest_id which is a list of binary ints that defines which
    (mutually exclusive) set the atom(s) belong to.
    
    If there is a single atom, the sets are: is_rest, is_nonrest
    If there is a set of atoms, the sets are: is_rest, is_inter, is_nonrest
    
    Example: if there is a single atom that is in the nonrest set, the rest_id is [0, 1]
    
    Arguments
    ---------
    atoms : set or int
        a set of hybrid atom indices or single atom
    rest_atoms : set or list
        a list (or list-like) of atoms whose interactions will be scaled by REST
    Returns
    -------
    rest_id : list
        list of binaries indicating which set the atom(s) belong to
    """

    if isinstance(atoms, int):
        rest_id = [0, 1] # Set the default rest_id to non-REST
        if atoms in rest_atoms:
            rest_id = [1, 0]
        return rest_id

    elif isinstance(atoms, set):
        rest_id = [0, 0, 1] # Set the default rest_id to non-REST
        if atoms.intersection(rest_atoms) != set(): # At least one of the atoms is REST
            if atoms.issubset(rest_atoms): # All atoms are REST
                rest_id = [1, 0, 0]
            else: # At least one (but not all) of the atoms is are REST
                rest_id = [0, 1, 0]
        return rest_id

    else:
        raise Exception(f"atoms is of type {type(atoms)}, but only `int` and `set` are allowable")


# Get vanilla system bond force
bond_force = system_forces['HarmonicBondForce']

# Set periodicity
if bond_force.usesPeriodicBoundaryConditions():
    rest_bond_force.setUsesPeriodicBoundaryConditions(True)

# Add bonds to rest_system
for term_idx in range(bond_force.getNumBonds()):
    # Get the bond parameters and rest id
    p1, p2, r0, k = bond_force.getBondParameters(term_idx)
    idx_set = set([p1, p2])
    rest_id = get_rest_identifier(idx_set, rest_atoms)

    # Add the bond
    bond_term = (p1, p2, rest_id + [r0, k])
    rest_bond_force.addBond(*bond_term)



# Define the custom expression
angle_expression = "rest_scale * (K / 2) * (theta - theta0)^2;"
angle_expression += "rest_scale = is_rest * lambda_rest_angles * lambda_rest_angles " \
                    "+ is_inter * lambda_rest_angles " \
                    "+ is_nonrest;"

# Create custom force
rest_angle_force = openmm.CustomAngleForce(angle_expression)
rest_system.addForce(rest_angle_force)

# Add global parameters
rest_angle_force.addGlobalParameter("lambda_rest_angles", 1.0)

# Add per-angle parameters for rest scaling
rest_angle_force.addPerAngleParameter("is_rest")
rest_angle_force.addPerAngleParameter("is_inter")
rest_angle_force.addPerAngleParameter("is_nonrest")

# Add per-angle parameters for defining angle energy
rest_angle_force.addPerAngleParameter('theta0')  # equilibrium angle 
rest_angle_force.addPerAngleParameter('K')  # force constant

# Get vanilla system angle force
angle_force = system_forces['HarmonicAngleForce']

# Set periodicity
if angle_force.usesPeriodicBoundaryConditions():
    rest_angle_force.setUsesPeriodicBoundaryConditions(True)

# Add angles to rest_system
for term_idx in range(angle_force.getNumAngles()):
    # Get the angle parameters and rest id
    p1, p2, p3, theta0, k = angle_force.getAngleParameters(term_idx)
    idx_set = set([p1, p2, p3])
    rest_id = get_rest_identifier(idx_set, rest_atoms)

    # Add the angle
    angle_term = (p1, p2, p3, rest_id + [theta0, k])
    rest_angle_force.addAngle(*angle_term)
# Define the custom expression
torsion_expression = "rest_scale * U;"
torsion_expression += "rest_scale = is_rest * lambda_rest_torsions * lambda_rest_torsions " \
                      "+ is_inter * lambda_rest_torsions " \
                      "+ is_nonrest;"
torsion_expression += "U = (K * (1 + cos(periodicity * theta - phase)));"

# Create custom force
rest_torsion_force = openmm.CustomTorsionForce(torsion_expression)
rest_system.addForce(rest_torsion_force)

# Add global parameters
rest_torsion_force.addGlobalParameter("lambda_rest_torsions", 1.0)

# Add per-torsion parameters for rest scaling
rest_torsion_force.addPerTorsionParameter("is_rest")
rest_torsion_force.addPerTorsionParameter("is_inter")
rest_torsion_force.addPerTorsionParameter("is_nonrest")

# Add per-torsion parameters for defining torsion energy
rest_torsion_force.addPerTorsionParameter('periodicity')
rest_torsion_force.addPerTorsionParameter('phase') # phase offset
rest_torsion_force.addPerTorsionParameter('K') # force constant

# Get vanilla system torsion force
torsion_force = system_forces['PeriodicTorsionForce']

# Set periodicity
if torsion_force.usesPeriodicBoundaryConditions():
    rest_torsion_force.setUsesPeriodicBoundaryConditions(True)

# Add torsions to rest_system
for torsion_idx in range(torsion_force.getNumTorsions()):
    # Get the torsion parameters and rest id
    p1, p2, p3, p4, periodicity, phase, K = torsion_force.getTorsionParameters(torsion_idx)
    idx_set = set([p1, p2, p3, p4])
    rest_id = get_rest_identifier(idx_set, rest_atoms)

    # Add torsion
    torsion_term = (p1, p2, p3, p4, rest_id + [periodicity, phase, K])
    rest_torsion_force.addTorsion(*torsion_term)


# Create nonbonded force
rest_nonbonded_force = openmm.NonbondedForce()
rest_system.addForce(rest_nonbonded_force)

# Get vanilla system nonbonded force
nonbonded_force = system_forces['NonbondedForce']

# Set the nonbonded method and related parameters
nonbonded_method = nonbonded_force.getNonbondedMethod()
rest_nonbonded_force.setNonbondedMethod(nonbonded_method)
if nonbonded_method != openmm.NonbondedForce.NoCutoff:
    epsilon_solvent = nonbonded_force.getReactionFieldDielectric()
    cutoff = nonbonded_force.getCutoffDistance()
    rest_nonbonded_force.setReactionFieldDielectric(epsilon_solvent)
    rest_nonbonded_force.setCutoffDistance(cutoff)
if nonbonded_method in [openmm.NonbondedForce.PME, openmm.NonbondedForce.Ewald]:
    [alpha_ewald, nx, ny, nz] = nonbonded_force.getPMEParameters()
    delta = nonbonded_force.getEwaldErrorTolerance()
    rest_nonbonded_force.setPMEParameters(alpha_ewald, nx, ny, nz)
    rest_nonbonded_force.setEwaldErrorTolerance(delta)

# Copy switching function from vanilla system
switch_bool = nonbonded_force.getUseSwitchingFunction()
rest_nonbonded_force.setUseSwitchingFunction(switch_bool)
if switch_bool:
    switching_distance = nonbonded_force.getSwitchingDistance()
    rest_nonbonded_force.setSwitchingDistance(switching_distance)

# Copy dispersion correction
dispersion_bool = nonbonded_force.getUseDispersionCorrection()
rest_nonbonded_force.setUseDispersionCorrection(dispersion_bool)

# Add global parameters
rest_nonbonded_force.addGlobalParameter('lambda_rest_electrostatics', 0.)
rest_nonbonded_force.addGlobalParameter('lambda_rest_sterics', 0.)

# Add nonbondeds to rest_system
for particle_idx in range(nonbonded_force.getNumParticles()):
    # Get the nonbonded parameters and rest id
    q, sigma, epsilon = nonbonded_force.getParticleParameters(particle_idx)
    rest_id = get_rest_identifier(particle_idx, rest_atoms)
    
    # Add particles and offsets
    if rest_id == [0, 1]: # nonrest
        rest_nonbonded_force.addParticle(q, sigma, epsilon)
    
    else: # rest
        rest_nonbonded_force.addParticle(q, sigma, epsilon)
        rest_nonbonded_force.addParticleParameterOffset('lambda_rest_electrostatics', particle_idx, q, 0.0*sigma, epsilon*0.0)
        rest_nonbonded_force.addParticleParameterOffset('lambda_rest_sterics', particle_idx, q*0.0, 0.0*sigma, epsilon)

# Handle exceptions
for exception_idx in range(nonbonded_force.getNumExceptions()):
    # Get exception parameters and rest id
    p1, p2, chargeProd, sigma, epsilon = nonbonded_force.getExceptionParameters(exception_idx)
    idx_set = set([p1, p2])
    rest_id = get_rest_identifier(idx_set, rest_atoms)
    
    # Add exceptions and offsets
    exc_idx = rest_nonbonded_force.addException(p1, p2, chargeProd, sigma, epsilon)
    if rest_id == [0, 0, 1]: # nonrest
        pass
    
    elif rest_id == [1, 0, 0]: # rest
        rest_nonbonded_force.addExceptionParameterOffset('lambda_rest_sterics', exc_idx, chargeProd, 0.0*sigma, epsilon)
        
    elif rest_id == [0, 1, 0]: # inter
        rest_nonbonded_force.addExceptionParameterOffset('lambda_rest_electrostatics', exc_idx, chargeProd, 0.0*sigma, epsilon)
















import math
import logging
import numpy as np
from openmmtools.constants import kB
from openmmtools import cache, mcmc, multistate
from openmmtools.multistate import ReplicaExchangeSampler
from openmmtools.states import GlobalParameterState, SamplerState, ThermodynamicState, CompoundThermodynamicState



class RESTState(GlobalParameterState):
    lambda_rest_bonds = GlobalParameterState.GlobalParameter('lambda_rest_bonds', standard_value=1.0)
    lambda_rest_angles = GlobalParameterState.GlobalParameter('lambda_rest_angles', standard_value=1.0)
    lambda_rest_torsions = GlobalParameterState.GlobalParameter('lambda_rest_torsions', standard_value=1.0)
    lambda_rest_electrostatics = GlobalParameterState.GlobalParameter('lambda_rest_electrostatics', standard_value=0.0)
    lambda_rest_sterics = GlobalParameterState.GlobalParameter('lambda_rest_sterics', standard_value=0.0)

    def set_rest_parameters(self, beta_m, beta_0):
        """Set all defined lambda parameters to the given value.

        The undefined parameters (i.e. those being set to None) remain undefined.

        Parameters
        ----------
        new_value : float
            The new value for all defined parameters.
        """
        lambda_functions = {'lambda_rest_bonds': lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_angles' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_torsions' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0),
                 'lambda_rest_electrostatics' : lambda beta_m, beta_0 : np.sqrt(beta_m / beta_0) - 1,
                 'lambda_rest_sterics' : lambda beta_m, beta_0 : beta_m / beta_0 - 1
                 }

        for parameter_name in self._parameters:
            if self._parameters[parameter_name] is not None:
                new_value = lambda_functions[parameter_name](beta_m, beta_0)
                setattr(self, parameter_name, new_value)

# Set temperatures for each thermodynamic state
n_replicas = 12  # Number of temperature replicas
T_min = 300 * kelvin  # Minimum temperature (i.e., temperature of desired distribution)
T_max = 600 * kelvin  # Maximum temperature
temperatures = [T_min + (T_max - T_min) * (math.exp(float(i) / float(n_replicas-1)) - 1.0) / (math.e - 1.0)
                for i in range(n_replicas)]

# Create reference thermodynamic state
rest_state = RESTState.from_system(rest_system)
thermostate = ThermodynamicState(rest_system, temperature=T_min)
compound_thermodynamic_state = CompoundThermodynamicState(thermostate, composable_states=[rest_state])

sampler_state =  SamplerState(inpcrd.positions, box_vectors=rest_system.getDefaultPeriodicBoxVectors())

#根据ReplicaExchangeSampler的构造函数，你需要通过MCMC moves的context_cache参数来传递平台属
#cache.global_context_cache.empty() 
#platform = Platform.getPlatformByName('CUDA')
#platform_properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
#cache.global_context_cache.set_platform(platform, platform_properties)

beta_0 = 1/(kB*T_min)
thermodynamic_state_list = []
sampler_state_list = []
for temperature in temperatures:
    # Create a thermodynamic state with REST interactions scaled to the given temperature
    beta_m = 1/(kB*temperature)
    compound_thermodynamic_state_copy = copy.deepcopy(compound_thermodynamic_state)
    compound_thermodynamic_state_copy.set_rest_parameters(beta_m, beta_0)
    thermodynamic_state_list.append(compound_thermodynamic_state_copy)

    # Generate a sampler_state with minimized positions
    context, integrator = cache.global_context_cache.get_context(compound_thermodynamic_state_copy)
    sampler_state.apply_to_context(context, ignore_velocities=True)
    openmm.LocalEnergyMinimizer.minimize(context)
    sampler_state.update_from_context(context)
    sampler_state_list.append(copy.deepcopy(sampler_state))

# Set up sampler
_logger = logging.getLogger()
_logger.setLevel(logging.DEBUG)
move = mcmc.LangevinDynamicsMove(timestep=2*femtoseconds, n_steps=500) # each move is 1 ps
simulation = ReplicaExchangeSampler(mcmc_moves=move, number_of_iterations=500, online_analysis_interval=10)

# Run repex
reporter_file = "test.nc"
reporter = multistate.MultiStateReporter(reporter_file, checkpoint_interval=10)
simulation.create(thermodynamic_states=thermodynamic_state_list,
                  sampler_states=sampler_state_list,
                  storage=reporter)
simulation.run()


import numpy as np 

reporter = multistate.MultiStateReporter("test.nc", open_mode="r") 

replica_states = reporter.read_replica_thermodynamic_states()
iterations = np.where(replica_states == 0)[0]
replicas = np.where(replica_states == 0)[1]

positions_300K = [] 
#available_checkpoints = [0, 100, 200]

for iter_idx, replica_idx in zip(iterations, replicas):  
    if iter_idx % 10 == 0: #"10"取自online_analysis_interval=10以及checkpoint_interval=10
#        print("num_exchange", iter_idx)
#        print("num_replica", replica_idx)
        positions_300K.append(reporter.read_sampler_states(iter_idx) [replica_idx].positions)

import mdtraj as md  
  
# 将positions转换为MDTraj轨迹  
positions_array = np.array([pos.value_in_unit(nanometers) for pos in positions_300K])  
traj = md.Trajectory(positions_array, topology=prmtop.topology)  
traj.save_dcd('rest_300K_trajectory.dcd')