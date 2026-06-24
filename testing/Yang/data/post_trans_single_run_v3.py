#!/usr/bin/env python3
try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
from sys import stdout, exit, stderr
import os, time, traceback
import parmed as pmd
import numpy as np
import mdtraj as mdt

usage = '''python post_trans_single_run_v2.py <psf file> <ncrst file> <prm file> <temperature> 
                                          <# CPUs> <outname> <random seed> <simulation step> 
                                          <2nd structure> <Q_threshold> <native cor> <gpu> 
                                          <restraint radius>
'''

###### convert time seconds to hours ######
def convert_time(seconds):
    return seconds/3600
###### END convert time seconds to hours ######



############## MAIN #################

psffile = sys.argv[1]
ncrstfile = sys.argv[2]
prmfile = sys.argv[3]
temp = float(sys.argv[4])
ppn = sys.argv[5]
outname = sys.argv[6]
rand = int(sys.argv[7])
sim_step = int(sys.argv[8])
secondary_structure_def = sys.argv[9]
Q_threshold = float(sys.argv[10])
native_cor = sys.argv[11]
    
cpfile = outname+'.ncrst'

timestep = 0.015*picoseconds
fbsolu = 0.05/picosecond
temp = temp*kelvin
nsteps_save = 5000



### contact map and distance map for start structure ###
native_cor = CharmmCrdFile(native_cor)
native_cor = native_cor.positions.value_in_unit(angstrom)
native_contact_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
native_distance_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
sec_strc_def = []


psf = CharmmPsfFile(psffile)
psf_pmd = pmd.load_file(psffile)
rst = pmd.load_file(ncrstfile)
forcefield = ForceField(prmfile)
top = psf.topology
# re-name residues that are changed by openmm
for resid, res in enumerate(top.residues()):
    if res.name != psf_pmd.residues[resid].name:
        res.name = psf_pmd.residues[resid].name
templete_map = {}
for chain in top.chains():
    for res in chain.residues():
        templete_map[res] = res.name
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=2.0*nanometer, constraints=AllBonds, 
        removeCMMotion=False, ignoreExternalBonds=True, 
        residueTemplates=templete_map)
for force in system.getForces():
    if force.getName() == 'CustomNonbondedForce':
        custom_nb_force = force
        break
custom_nb_force.setUseSwitchingFunction(True)
custom_nb_force.setSwitchingDistance(1.8*nanometer)


integrator = LangevinIntegrator(temp, fbsolu, timestep)
integrator.setConstraintTolerance(0.00001)
integrator.setRandomNumberSeed(rand)

# prepare simulation
dev_index = 0
properties = {'CudaPrecision': 'mixed', "DeviceIndex": "%d" % dev_index}
platform = Platform.getPlatformByName('CUDA')


simulation = Simulation(top, system, integrator, platform, properties)
simulation.context.setPositions(rst.coordinates[0]*angstrom)
simulation.context.setVelocitiesToTemperature(temp)


# append reporters
simulation.reporters = []
simulation.reporters.append(DCDReporter(outname+'.dcd', nsteps_save, append=False))
# simulation.reporters.append(pmd.openmm.reporters.RestartReporter(cpfile, nsteps_save, netcdf=True))
simulation.reporters.append(
    StateDataReporter(f"{outname}.log", nsteps_save, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                        totalEnergy=True, temperature=True, speed=True, separator='\t'))


# run production simulation

simulation.step(sim_step)
