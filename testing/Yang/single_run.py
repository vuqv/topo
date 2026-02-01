#!/usr/bin/env python3
import getopt
import time
from distutils.util import strtobool

import numpy as np
import parmed as pmd
from openmm import *
from openmm.app import *
from openmm.unit import *

usage = '\nUsage: python single_run.py -f control_file\n'


###### convert time seconds to hours ######
def convert_time(seconds):
    return seconds / 3600


###### END convert time seconds to hours ######


###### calculate native contact fraction ######
def calc_Q(current_cor):
    global native_contact_map, native_distance_map, native_contact_num, sec_strc_def, sdist
    current_contact_num = 0
    for i in range(len(current_cor) - 4):
        tag_i = 0
        for rs in sec_strc_def:
            if rs[0] - 1 <= i <= rs[1] - 1:
                tag_i = 1
                break
        if tag_i == 0:
            continue
        for j in range(i + 4, len(current_cor)):
            tag_j = 0
            for rs in sec_strc_def:
                if rs[0] - 1 <= j <= rs[1] - 1:
                    tag_j = 1
                    break
            if tag_j == 0:
                continue
            if native_contact_map[i][j] == 1:
                dist = pow(pow(current_cor[i][0] - current_cor[j][0], 2) + pow(current_cor[i][1] - current_cor[j][1], 2)
                           + pow(current_cor[i][2] - current_cor[j][2], 2), 0.5)
                if dist <= sdist * native_distance_map[i][j]:
                    current_contact_num += 1
    return current_contact_num / native_contact_num


###### END calculate native contact fraction ######


###### Q mod filter ######
def calc_Q_mod(Q_ts):
    edges = np.arange(0, 1.02, 0.02)
    (N, be) = np.histogram(Q_ts, bins=edges)
    idx = np.argwhere(N == np.max(N))[0]
    Q_mod = (edges[idx] + edges[idx + 1]) / 2
    return Q_mod


###### END Q mod filter ######


############## MAIN #################
ctrlfile = ''
if len(sys.argv) == 1:
    print(usage)
    sys.exit()

try:
    opts, args = getopt.getopt(sys.argv[1:], "hf:", ["ctrlfile="])
except getopt.GetoptError:
    print(usage)
    sys.exit()
for opt, arg in opts:
    if opt == '-h':
        print(usage)
        sys.exit()
    elif opt in ("-f", "--ctrlfile"):
        ctrlfile = arg
"""
Reading control file to load setup
"""
if not os.path.exists(ctrlfile):
    print('Error: cannot find control file ' + ctrlfile + '.')
    sys.exit()

file_object = open(ctrlfile, 'r')
try:
    for line in file_object:
        line = line.strip()
        if not line:
            # This is a blank line
            continue
        if line.startswith('#'):
            # This is a comment line
            continue
        if line.startswith('psffile'):
            words = line.split('=')
            psffile = words[1].strip()
            continue
        if line.startswith('corfile'):
            words = line.split('=')
            corfile = words[1].strip()
            continue
        if line.startswith('prmfile'):
            words = line.split('=')
            prmfile = words[1].strip()
            continue
        if line.startswith('secondary_structure_def'):
            words = line.split('=')
            secondary_structure_def = words[1].strip()
            continue
        if line.startswith('temp_heating'):
            words = line.split('=')
            temp_heating = float(words[1].strip()) * kelvin
            continue
        if line.startswith('temp_prod'):
            words = line.split('=')
            temp_prod = float(words[1].strip()) * kelvin
            continue
        if line.startswith('ppn'):
            words = line.split('=')
            ppn = str(words[1].strip())
            continue
        if line.startswith('outname'):
            words = line.split('=')
            outname = words[1].strip()
            continue
        if line.startswith('mdsteps'):
            words = line.split('=')
            mdsteps = int(words[1].strip())
            continue
        if line.startswith('heating_steps'):
            words = line.split('=')
            heating_steps = int(words[1].strip())
            continue
        if line.startswith('Q_threshold'):
            words = line.split('=')
            Q_threshold = float(words[1].strip())
            continue
        if line.startswith('native_cor'):
            words = line.split('=')
            native_cor = words[1].strip()
            continue
        if line.startswith('nstxout'):
            words = line.split('=')
            nstxout = int(words[1].strip())
            continue
        if line.startswith('nstlog'):
            words = line.split('=')
            nstlog = int(words[1].strip())
            continue
        if line.startswith('nsteps_calQ'):
            words = line.split('=')
            nsteps_calQ = int(words[1].strip())
            continue
        if line.startswith('restart'):
            words = line.split('=')
            restart = bool(strtobool(words[1].strip()))
            continue
        if line.startswith('use_gpu'):
            words = line.split('=')
            # use_gpu = bool(words[1].strip())
            use_gpu = bool(strtobool(words[1].strip()))
            continue
finally:
    file_object.close()

# checkpoint file
cpfile = outname + '.chk'

timestep = 0.015 * picoseconds
fbsolu = 0.05 / picosecond
half_window = 100

dist_cutoff = 8  # distance cutoff for finding native contact
sdist = 1.2  # multiple factor of native distance to determine native contact in trajectory
fold_nframe = 100  # number of frames to determine folding status

Q_list = []

### contact map and distance map for start structure ###
native_cor = CharmmCrdFile(native_cor)
native_cor = native_cor.positions.value_in_unit(angstrom)
native_contact_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
native_distance_map = [[0 for j in range(len(native_cor))] for i in range(len(native_cor))]
sec_strc_def = []
sec_def_object = open(secondary_structure_def, 'r')
for line in sec_def_object.readlines():
    line = line.strip()
    if line != '':
        words = line.split()
        sec_strc_def.append([int(words[1]), int(words[2])])
for i in range(len(native_cor) - 4):
    for j in range(i + 4, len(native_cor)):
        dist = pow(pow(native_cor[i][0] - native_cor[j][0], 2) + pow(native_cor[i][1] - native_cor[j][1], 2)
                   + pow(native_cor[i][2] - native_cor[j][2], 2), 0.5)
        if dist <= dist_cutoff:
            native_contact_map[i][j] = 1
            native_distance_map[i][j] = dist
            native_contact_map[j][i] = 1
            native_distance_map[j][i] = dist

sec_def_object.close()
native_contact_num = 0
for i in range(len(native_cor) - 4):
    tag_i = 0
    for rs in sec_strc_def:
        if rs[0] - 1 <= i <= rs[1] - 1:
            tag_i = 1
            break
    if tag_i == 0:
        continue
    for j in range(i + 4, len(native_cor)):
        tag_j = 0
        for rs in sec_strc_def:
            if rs[0] - 1 <= j <= rs[1] - 1:
                tag_j = 1
                break
        if tag_j == 0:
            continue
        if native_contact_map[i][j] == 1:
            native_contact_num += 1
### END contact map and distance map for start structure ###

psf = CharmmPsfFile(psffile)
psf_pmd = pmd.load_file(psffile)
cor = CharmmCrdFile(corfile)
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
system = forcefield.createSystem(top, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=2.0 * nanometer,
                                 constraints=AllBonds, removeCMMotion=False,
                                 ignoreExternalBonds=True, residueTemplates=templete_map)
# custom_nb_force = system.getForce(4)
for force in system.getForces():
    if force.getName() == 'CustomNonbondedForce':
        # custom_nb_force = force
        force.setUseSwitchingFunction(True)
        force.setSwitchingDistance(1.8 * nanometer)
# Done for system initialization

# prepare simulation
if use_gpu:
    print("Running simulation on CUDA device")
    dev_index = 0
    properties = {'CudaPrecision': 'mixed', "DeviceIndex": "%d" % dev_index}
    platform = Platform.getPlatformByName('CUDA')
else:
    print("Running simulation on CPU")
    properties = {'Threads': ppn}
    platform = Platform.getPlatformByName('CPU')

if restart:

    integrator = LangevinIntegrator(temp_prod, fbsolu, timestep)
    integrator.setConstraintTolerance(0.00001)
    simulation = Simulation(top, system, integrator, platform, properties)
    # load checkpoint
    simulation.loadCheckpoint(cpfile)
    timeCount = simulation.context.getState().getTime()
    stepCount = simulation.context.getState().getStepCount()
    print(f"Restart from checkpoint, Time = {timeCount}, Step= {stepCount}")
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(cpfile, nstxout))
    simulation.reporters.append(DCDReporter(outname + '_prod.dcd', nstxout, append=True))
    simulation.reporters.append(
        StateDataReporter(f'{outname}.log', nstlog, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                          totalEnergy=True, temperature=True, speed=True, separator='\t',
                          append=True))

    out_mode = 'a'

else:
    # if not restart= run from beginning => need to run heating phase
    print(f"Heating the system to {temp_heating} K for {heating_steps} steps")
    integrator = LangevinIntegrator(temp_heating, fbsolu, timestep)
    integrator.setConstraintTolerance(0.00001)
    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(cor.positions)
    simulation.context.setVelocitiesToTemperature(temp_heating)
    simulation.reporters = []
    simulation.reporters.append(DCDReporter(outname + '_heating.dcd', nstxout, append=False))
    simulation.reporters.append(
        StateDataReporter(f'{outname}_heating.log', nstlog, step=True, time=True, potentialEnergy=True,
                          kineticEnergy=True,
                          totalEnergy=True, temperature=True, speed=True, separator='\t'))
    simulation.step(heating_steps)
    last_heating_positions = simulation.context.getState(getPositions=True).getPositions()
    # save checkpoint at the end of heating phase
    # simulation.saveCheckpoint(outname + '_heating.chk')
    # Done for heating phase, move to production phase
    print(f"Production run for folding at {temp_prod} K for {mdsteps} steps")
    integrator = LangevinIntegrator(temp_prod, fbsolu, timestep)
    integrator.setConstraintTolerance(0.00001)
    simulation = Simulation(top, system, integrator, platform, properties)
    simulation.context.setPositions(last_heating_positions)
    simulation.context.setVelocitiesToTemperature(temp_prod)
    # simulation.loadCheckpoint(outname + '_heating.chk')  # load state from checkpoint in heating phase
    simulation.reporters = []
    simulation.reporters.append(CheckpointReporter(cpfile, nstxout))
    simulation.reporters.append(DCDReporter(outname + '_prod.dcd', nstxout, append=False))
    simulation.reporters.append(
        StateDataReporter(f'{outname}.log', nstlog, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
                          totalEnergy=True, temperature=True, speed=True, separator='\t'))
    out_mode = 'w'

# run production simulation
start_time = time.time()

# production run
print(
    "###########################################\n\nProduction Phase... \n###########################################")
with open('fQ.dat', out_mode) as f:
    if not restart:
        f.write('#    Time (ns)             step      fQ         Q_mod Folding\n')

    while simulation.context.getState().getStepCount() < mdsteps:
        simulation.step(nsteps_calQ)
        current_cor = simulation.context.getState(getPositions=True).getPositions().value_in_unit(angstrom)
        Q = calc_Q(current_cor)
        if len(Q_list) < 2 * half_window + 1:
            Q_list.append(Q)
            Q_mod_value = -1.0
            folding = False
        else:
            Q_list.pop(0)
            Q_list.append(Q)
            Q_mod = calc_Q_mod(np.array(Q_list))
            # Q_mod = '%f' % Q_mod
            Q_mod_value = Q_mod.item(0)
            folding = Q_mod_value > Q_threshold

        # if len(Q)
        f.write(
            f'{simulation.context.getState().getTime().value_in_unit(nanosecond):10.6f} {simulation.context.getState().getStepCount():20d} {Q:10.4f} {Q_mod_value:10.4f} {folding}\n')
        f.flush()
# save checkpoint for the last state, before simulation is terminated
simulation.saveCheckpoint(cpfile)
