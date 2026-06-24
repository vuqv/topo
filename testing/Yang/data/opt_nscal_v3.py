#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimize domain/interface native-contact scaling factors (nscal) for a CG Go-like
protein model by running multiple trajectories and checking domain stability.

This script iterates over a discrete set of candidate nscal "levels" for each
domain and each domain-domain interface. At each round, it:

- Writes domain/interface scale factors
- Builds a coarse-grained (CG) model (Ca or Ca+SCM)
- Runs `ntraj` independent MD trajectories
- Computes per-domain native-contact fractions (Q / qbb values)
- Marks domains/interfaces as stable/unstable and increases nscal level for
  unstable ones

Key outputs (written under `outpath`):

- `opt_nscal.log`: human-readable record of rounds, chosen nscal values, and
  stability probabilities.
- `round_*/`: one folder per nscal-level iteration.
  - `round_*/setup/`: model files and trajectory outputs (e.g. `*.psf`, `*.cor`,
    `*.top`, `*.prm`/`*.xml`, `*.dcd`, `qbb_*.dat`, etc. depending on scripts).
- `setup/secondary_struc_defs.txt`: secondary-structure elements (SSEs) derived
  from STRIDE and mapped to the cleaned PDB's residue indices.

External dependencies / assumptions:

- Python packages: `parmed`, `numpy`
- External binaries in PATH: `stride`
- External scripts co-located with this script:
  - `create_cg_protein_model.py`
  - `post_trans_single_run_v3.py`
  - `calc_native_contact_fraction.pl` (invoked via `perl`)
  - `parse_cg_prm.py`

Important behavioral notes:

- This script uses `os.chdir(...)` heavily; many paths are implicitly relative
  to the current working directory (CWD).
- Most work is executed via `os.system(...)` / `os.popen(...)`. Errors from
  downstream scripts may not raise Python exceptions; check log files.

Original author: Yang Jiang (PSU)
Refactor: 2024-10-10 by Ian SitariK (PSU)
"""

import os, time, traceback, io, sys, getopt, multiprocessing, random
import parmed as pmd
import numpy as np
import argparse
import subprocess
from pathlib import Path

################### Functions #########################
def parse_nscal_levels(nscal_level_file):
    """Parse a user-supplied nscal "level" file.

    File format (whitespace-delimited):

    - Lines beginning with one of `a`, `b`, `c`, `i` define the candidate levels
      for that class.
    - The remainder of the line are floats (e.g. `a 1.1 1.3 1.5`).

    Where:
    - `a`: alpha-helix domain
    - `b`: beta-sheet domain
    - `c`: mixed alpha/beta domain
    - `i`: interface between two domains

    Returns
    -------
    dict[str, list[float]]
        Mapping of class -> list of candidate nscal levels.

    Notes
    -----
    All four classes must be present and have the same number of levels,
    otherwise the script exits.
    """
    nscal_set = {"a": [],
                 "b": [],
                 "c": [],
                 "i": []}
    f = open(nscal_level_file)
    lines = f.readlines()
    f.close()
    for line in lines:
        words = line.strip().split()
        if words[0] in list(nscal_set.keys()):
            nscal_set[words[0]] = [float(n) for n in words[1:]]
    
    for (k, v) in nscal_set.items():
        if len(v) == 0:
            print('Error: No nscal levels found for class "%s" in the file %s.'%(k, nscal_level_file))
            sys.exit()
        elif len(v) != len(nscal_set["a"]):
            print("Error: Inconsistent number of nscal levels found in the file %s."%(nscal_level_file))
            sys.exit()
    
    return nscal_set

def parse_domain(dom_def_file):
    """Parse a domain definition file and derive interface entries.

    Input file format (one domain per non-comment line):

        <start:end> [<start:end> ...] <class>    # optional comment

    Examples:

        10:55 a
        60:110 120:150 b
        5:40 45:80 c

    Where:
    - `<class>` must be one of `a`, `b`, `c`
    - Each `<start:end>` is an *inclusive* residue-number range from the PDB

    This parser also normalizes residue indices to start at 1 by subtracting
    `(min_start - 1)` across all ranges.

    Returns
    -------
    list[dict]
        A list of entries. The first `N` are domain entries:
        - `{"range": [[start, end], ...], "class": "a"|"b"|"c"}`
        After that, all pairwise interfaces are appended as:
        - `{"range": [i, j], "class": "i"}`
        where `i` and `j` are 0-based indices into the original domain list.
    """
    domain = []
    domian_class = {"a": "Alpha-helix",
                    "b": "Beta-sheet",
                    "c": "Alpha-Beta"}
    f = open(dom_def_file)
    lines = f.readlines()
    f.close()
    
    n = 0 # domain index
    start_min = np.inf
    for line in lines:
        line = line.strip()
        if line == "" or line.startswith("#"):
            continue
        data = line.split("#")[0].strip().split()
        res_range = data[:-1]
        sec_class = data[-1]
        out = "   Domain %d: "%(n+1)
        if sec_class not in list(domian_class.keys()):
            print("Error: Forgot to specify class as a, b or c in domain definition?")
            sys.exit()
        dom = {"range": [],
               "class": sec_class}
        for d in res_range:
            [start, end] = d.split(':')
            out += "%s ~ %s; "%(start, end)
            if int(start) < start_min:
                start_min = int(start)
            dom["range"].append([int(start), int(end)])
        print(out + "Class: %s"%domian_class[sec_class])
        domain.append(dom)
        n += 1
    
    offset = start_min - 1
    for d in domain:
        for r in d["range"]:
            r[0] -= offset
            r[1] -= offset
    
    ndomain = len(domain)
    for i in range(ndomain-1):
        for j in range(i+1, ndomain):
            interface = {"range": [i, j],
                         "class": "i"}
            domain.append(interface)
    return domain

def clean_PDB(input_pdb):
    """Clean and renumber a PDB for downstream CG model generation.

    Operations:
    - Renumber residues sequentially starting at 1 (based on residue order
      inside the loaded structure).
    - Keep only standard amino acids and common histidine protonation states.
    - Save `<basename>_clean.pdb` in the current directory.

    Parameters
    ----------
    input_pdb : str
        Path to the original PDB file.

    Returns
    -------
    str
        Filename of the cleaned PDB (relative to CWD).
    """
    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP'];
    print("-> Cleaning PDB file %s"%input_pdb)
    name = input_pdb.split('/')[-1].split('.pdb')[0]
    struct = pmd.load_file(input_pdb)
    sel_idx = np.zeros(len(struct.atoms))
    for idx, res in enumerate(struct.residues):
        # Force sequential renumbering; later tools assume residue numbers start at 1.
        res.number = idx+1
        if res.name in AA_name_list:
            for atm in res.atoms:
                sel_idx[atm.idx] = 1
    struct[sel_idx].save(name+'_clean.pdb', overwrite=True, altlocs="occupancy")
    print("   PDB file cleaned")
    return name+'_clean.pdb'
    
def creat_CG_model(pdb_file, casm, domain):
    """Create a CG model using `create_cg_protein_model.py`.

    Side effects (in current working directory):
    - Creates and enters a `create_model/` folder.
    - Writes `domain_def.dat` and `go_model.cntrl`.
    - Runs the external CG model builder script and saves output to `go_model.log`.
    - Copies the generated model files (`*.psf`, `*.cor`, `*.top`, `*.prm`)
      back to the parent directory and returns to the parent directory.

    Parameters
    ----------
    pdb_file : str
        Path to PDB file (often relative to where `create_model/` is created).
    casm : int
        0 = Ca model (potential "bt"); 1 = Ca-SCM model (potential "mj").
    domain : list[dict]
        Domain/interface list (from `parse_domain`) with `dom['nscal']` already set.

    Returns
    -------
    (str, str)
        `(prefix, prm_name)` where `prefix` is the model prefix and `prm_name`
        is the CHARMM parameter filename produced by the CG model builder.

    Notes
    -----
    This function relies on the global `full_path` (script directory) to locate
    `create_cg_protein_model.py`.
    """
    if casm == 0:
        print("-> Creating Ca model")
        potential_name = "bt"
    else:
        print("-> Creating Ca-SCM model")
        potential_name = "mj"
    # os.system('mkdir create_model')
    Path("create_model").mkdir(exist_ok=True)
    os.chdir('create_model')
    
    # Count only true domains (exclude interfaces) to size the nscal matrix.
    ndomain = 0
    for d in domain:
        if d['class'] != 'i':
            ndomain += 1
    dom_idx = 0
    dom_idx_map = []
    nscal_array = np.zeros((ndomain,ndomain))
    f = open('domain_def.dat', 'w')
    for d in domain:
        if d['class'] == 'i':
            # Interface entries store indices into the original domain list.
            nscal_array[d['range'][0], d['range'][1]] = d['nscal']
            nscal_array[d['range'][1], d['range'][0]] = d['nscal']
        else:
            nscal_array[dom_idx, dom_idx] = d['nscal']
            for a in d['range']:
                # The CG model builder expects these repeated in order:
                # domain = start-end
                # scale factor = <nscal>
                f.write('domain = %d-%d\n'%(a[0], a[1]))
                f.write('scale factor = %f\n'%d['nscal'])
                dom_idx_map.append(dom_idx)
            dom_idx += 1
    ndomain = len(dom_idx_map)
    for i in range(ndomain-1):
        for j in range(i+1, ndomain):
            # Write the interface scale factors for each pair of domain fragments.
            f.write('scale factor = %f\n'%nscal_array[dom_idx_map[i],dom_idx_map[j]])
    f.close()
    
    ## write control file for making CG model
    f = open('go_model.cntrl', 'w')
    f.write('''pdbfile = ../'''+pdb_file+'''
casm = '''+str(casm)+'''
potential_name = '''+potential_name+'''
domain_file = domain_def.dat
''')
    f.close()

    ## Launch CG script to make model
    CG_script = os.path.join(full_path, 'create_cg_protein_model.py')
    print(f'python {CG_script} -f go_model.cntrl > go_model.log')
    os.system(f"python {CG_script} -f go_model.cntrl > go_model.log 2>&1")
    
    name = pdb_file.split('/')[-1].split('.pdb')[0]
    if casm == 0:
        prefix = name+'_ca'
    else:
        prefix = name+'_ca-cb'
    prm_name = name + '_nscal1_fnn1_go_' + potential_name + '.prm'
    
    if os.path.exists(prefix+'.psf'):
        os.system('cp *.psf ../')
        os.system('cp *.cor ../')
        os.system('cp *.top ../')
        os.system('cp *.prm ../')
        os.chdir('../')
    else:
        print("Error: failed to create CG model from %s\n"%pdb_file)
        sys.exit()
    return (prefix, prm_name)
    
def get_secondary_structure(pdb):
    """Run STRIDE to extract helices/strands and write SSE definitions.

    This function creates `setup/secondary_struc_defs.txt` in the current working
    directory. Each line is:

        <sse_index> <start_resid> <end_resid>

    where residue indices are 1-based indices into the *cleaned PDB structure*
    as loaded by ParmEd.

    Parameters
    ----------
    pdb : str
        Path to a PDB file (typically the cleaned PDB).
    """
    print("-> Getting secondary structure information")
    
    screen_out = os.popen('stride '+pdb).readlines()
    if not screen_out[0].startswith('REM  -------------------- Secondary structure summary -------------------  ~~~~'):
        print(''.join(screen_out))
        sys.exit()

    pdb_struct = pmd.load_file(pdb)
    sec_ele_list = []
    for line in screen_out:
        line = line.strip()
        if line.startswith('LOC '):
            sec_name = line[5:17].strip()
            if 'Helix' in sec_name or 'Strand' in sec_name:
                # STRIDE records residue numbers per chain; we map them to ParmEd residue indices.
                chainid = line[28]
                start_resnum = int(line[21:27].strip())
                end_resnum = int(line[39:45].strip())
                length = end_resnum - start_resnum + 1
                if length >= 4:
                    start_resid = np.nan
                    end_resid = np.nan
                    for res in pdb_struct.residues:
                        if res.chain == chainid and res.number == start_resnum:
                            start_resid = res.idx
                        elif res.chain == chainid and res.number == end_resnum:
                            end_resid = res.idx
                        if not np.isnan(start_resid) and not np.isnan(end_resid):
                            break
                    if np.isnan(start_resid) or np.isnan(end_resid):
                        print('Error: Cannot find residue %d or %d in chain %s in %s.'%(start_resnum, end_resnum, chainid, pdb))
                        sys.exit()
                    sec_ele_list.append([start_resid+1, end_resid+1])
    sec_ele_list.sort(key=lambda x: x[0])
    
    # Downstream scripts expect this file at `setup/secondary_struc_defs.txt`.
    f = open('setup/secondary_struc_defs.txt', 'w')
    for idx, sec_ele in enumerate(sec_ele_list):
        f.write('%d %d %d\n'%(idx+1, sec_ele[0], sec_ele[1]))
    f.close()
    print("   Done.")

def convert_time(t):
    """Convert seconds to `H:M:S` (or pass through NaN).

    Parameters
    ----------
    t : float
        Seconds.

    Returns
    -------
    str
        Formatted time string.
    """
    if np.isnan(t):
        return str(t)
    h = int(t / 3600)
    m = int((t - h * 3600) / 60)
    s = int(t - h * 3600 - m * 60)
    return "%d:%d:%d"%(h,m,s)

# def create_slurm_script(traj_idx, MD_cmd, script_path, job_name_prefix="MD_traj"):
#     """Create a SLURM script for a single MD trajectory.

#     Parameters
#     ----------
#     traj_idx : int
#         Trajectory index (1-based).
#     MD_cmd : str
#         The MD command to execute.
#     script_path : str
#         Path where the SLURM script will be written.
#     job_name_prefix : str
#         Prefix for the SLURM job name.

#     Returns
#     -------
#     str
#         Path to the created SLURM script.
#     """
#     script_content = f"""#!/bin/bash

# #SBATCH -J {job_name_prefix}_{traj_idx}
# #SBATCH -o traj_{traj_idx}.out
# #SBATCH -e traj_{traj_idx}.err

# #SBATCH --partition=standard
# #SBATCH --account=epo2_cr_default
# #SBATCH --gres=gpu:1
# #SBATCH -N 1
# #SBATCH -n 1
# #SBATCH --mem=10G
# #SBATCH -t 7-12:00:00

# cd $SLURM_SUBMIT_DIR

# module load cuda/11.5.0

# # Run MD simulation
# {MD_cmd}
# """
#     with open(script_path, 'w') as f:
#         f.write(script_content)
#     os.chmod(script_path, 0o755)  # Make executable
#     return script_path

def create_slurm_script(
    traj_idx,
    MD_cmd,
    script_path,
    job_name_prefix="MD_traj",
    partition="standard",
):
    """Create a SLURM script for a single MD trajectory.

    Parameters
    ----------
    traj_idx : int
        Trajectory index (1-based).
    MD_cmd : str
        The MD command to execute.
    script_path : str
        Path where the SLURM script will be written.
    job_name_prefix : str
        Prefix for the SLURM job name.
    partition : str
        Partition choice. Supported: "standard", "mgc".

    Returns
    -------
    str
        Path to the created SLURM script.
    """
    partition_map = {
        "standard": {
            "partition": "standard",
            "account": "epo2_cr_default",
        },
        "mgc": {
            "partition": "mgc-nih",
            "account": "epo2_nih",
        },
    }

    if partition not in partition_map:
        raise ValueError(
            f"Unsupported partition '{partition}'. "
            f"Choose from: {list(partition_map.keys())}"
        )

    slurm_cfg = partition_map[partition]

    script_content = f"""#!/bin/bash

#SBATCH -J {job_name_prefix}_{traj_idx}
#SBATCH -o traj_{traj_idx}.out
#SBATCH -e traj_{traj_idx}.err

#SBATCH --partition={slurm_cfg['partition']}
#SBATCH --account={slurm_cfg['account']}
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -t 7-12:00:00

cd $SLURM_SUBMIT_DIR

##module load cuda/11.5.0

## Run MD simulation
{MD_cmd}
"""

    with open(script_path, "w") as f:
        f.write(script_content)

    os.chmod(script_path, 0o755)
    return script_path

def wait_for_slurm_jobs(job_ids, poll_interval=300):
    """Wait for all SLURM jobs to complete by polling squeue.

    Parameters
    ----------
    job_ids : list[str]
        List of SLURM job IDs to wait for.
    poll_interval : int
        Seconds between polling attempts.

    Returns
    -------
    bool
        True if all jobs completed successfully, False otherwise.
    """
    if not job_ids:
        return True
    
    print(f"Waiting for {len(job_ids)} SLURM jobs to complete...")
    
    while True:
        # Check if any jobs are still running by querying each job ID
        running_jobs = []
        for job_id in job_ids:
            result = subprocess.run(
                ['squeue', '-j', job_id, '-h', '--noheader'],
                capture_output=True,
                text=True
            )
            # If squeue returns output, the job is still running
            if result.stdout.strip():
                running_jobs.append(job_id)
        
        if not running_jobs:
            print("All SLURM jobs completed.")
            break
        
        print(f"  {len(running_jobs)} job(s) still running (IDs: {', '.join(running_jobs[:5])}{'...' if len(running_jobs) > 5 else ''})...")
        time.sleep(poll_interval)
    
    return True    

def run_simulation(idx, prefix, prm_name, rand):
    """(Legacy / WIP) Run one trajectory and compute Q.

    Notes
    -----
    - This function is not used by the current `main()` loop (which runs
      trajectories inline).
    - It references variables expected to exist in the global scope
      (`use_gpu`, `worker_idx`, `full_path`, `temperature`, `ppn`, `sim_step`,
      `dom_def_file`).
    - It currently contains an unconditional `quit()` after printing commands,
      which will terminate the process if called.
    """
    if use_gpu == 1:
        gpu = int(multiprocessing.current_process().name.split('-')[-1])-1-worker_idx
    else:
        gpu = -1
    MD_script = os.path.join(full_path, 'post_trans_single_run_v3.py')
    MD_cmd = f'{MD_script} setup/%s.psf setup/%s.cor setup/%s %f %d %d %d %d ../setup/secondary_struc_defs.txt 1.1 setup/%s.cor %d'%(prefix, prefix, prm_name, temperature, ppn, idx+1, rand, sim_step, prefix, gpu)
    print(MD_cmd)
    Q_script = os.path.join(full_path, 'calc_native_contact_fraction.pl')
    Q_cmd = f'{Q_script} -i setup/%s.cor -d ../%s -s ../setup/secondary_struc_defs.txt -t %d.dcd -r 1'%(prefix, dom_def_file, idx+1)
    print(Q_cmd)
    quit()
    os.system(f'{MD_script} setup/%s.psf setup/%s.cor setup/%s %f %d %d %d %d ../setup/secondary_struc_defs.txt 1.1 setup/%s.cor %d'%(prefix, prefix, prm_name, temperature, ppn, idx+1, rand, sim_step, prefix, gpu))
    os.system(f'{Q_script} -i setup/%s.cor -d ../%s -s ../setup/secondary_struc_defs.txt -t %d.dcd -r 1'%(prefix, dom_def_file, idx+1))

########################### MAIN #########################################
def main():
    """Entry point for nscal optimization.

    High-level algorithm:
    - Parse domains and auto-generate all interfaces.
    - Initialize an nscal "level index" for each domain/interface (starting at 0).
    - For each nscal level round:
      - Build a CG model with the current nscal values.
      - Run `ntraj` trajectories and compute Q per domain/interface.
      - If any domain/interface is unstable across trajectories, increment its
        nscal level index and continue to next round; otherwise stop.
    - Append final nscal values to `opt_nscal.log`.
    """

    script_name = f'opt_nscal'
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-i", "--input", type=str, default="", required=True, help="<input.pdb> for CG model creation")
    parser.add_argument("-d", "--domain", type=str, default="", required=True, help="<domain.dat> for domain defination. full path")
    parser.add_argument("-o", "--outpath", type=str, default="./", required=True, help="Path to output directory")
    parser.add_argument("-t", "--temp", type=int, default=310, required=False, help="<Temperature> in Kelvin")
    parser.add_argument("-p", "--tpn", type=int, default=10, required=False, help="<total number of CPUs>. Default 10.")
    parser.add_argument("-j", "--ntraj", type=int, default=10, required=False, help="<number of trajectories>. Default 10. -1 use GPU")
    parser.add_argument("-r", "--restart", type=int, default=0, required=False, help="<0 or 1> restart optimization. Default 0, not restart.")
    parser.add_argument("-s", "--nscal", type=str, default="", required=False, help="<nscal_level.dat> for nscal levels. Default values were obtained from a training set of 18 small single-domain proteins.")
    parser.add_argument("-c", "--casm", type=int, default=0, required=False, help="<0 or 1> CG model type. Default 0, C-alpha model. 1, C-alpha side chain model.")
    parser.add_argument("-n", "--pname", type=str, default="MD", required=False, help="Name of child job")
    parser.add_argument("-q", "--partition", type=str, default="standard", required=False, help="Partition profile for generated SLURM scripts: 'standard' or 'mgc'")
    args = parser.parse_args()

    ## make output folder
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
        print(f'Made directory: {args.outpath}')

    # User inputs / configuration.
    input_pdb = args.input
    dom_def_file = args.domain
    domain = parse_domain(dom_def_file)
    outpath = args.outpath
    temperature = args.temp
    tpn = args.tpn
    ntraj = args.ntraj
    restart = args.restart
    nscal_level_file = args.nscal
    casm = args.casm # use C-alpha model
    prefix_name = args.pname
    partition_name = args.partition
    # Print argument names and values
    for arg_name, arg_value in vars(args).items():
        print(f"{arg_name}: {arg_value}")
        
    # CPU/GPU bookkeeping:
    # - In this codebase, `tpn == -1` is used as a flag to use GPU mode.
    if tpn == -1:
        use_gpu = 1
        ppn = 1
    else:
        use_gpu = 0
        ppn = 1

    nproc = int(tpn/ppn) # Total number of processors
    # Default candidate nscal values (training set derived). Can be overridden by `-s`.
    nscal_set = {"a": [1.1954, 1.4704, 1.7453, 2.0322, 2.5044, 1.7453],
                 "b": [1.4732, 1.8120, 2.1508, 2.5044, 2.5044, 2.1508],
                 "c": [1.1556, 1.4213, 1.6871, 1.9644, 2.5044, 1.6871],
                 "i": [1.2747, 1.5679, 1.8611, 2.1670, 2.5044, 1.8611]}
    if nscal_level_file != "":
        nscal_set = parse_nscal_levels(nscal_level_file)
    print(f'nscal_set:\n{nscal_set}')

    #sim_step = 66666667 #for 1000 ns
    # sim_step = 6666667 #for 100 ns
    # sim_step = 16666667 # 250ns
    sim_step =  33333333 # 500ns
    # sim_step =  66666 # testing-1ns
    # Stability thresholds:
    # - A frame is considered folded for a domain if Q > Q_threshold (or Q == -1).
    # - A trajectory is considered stable for a domain if folded frames fraction >= frame_threshold.
    Q_threshold = 0.6688
    frame_threshold = 0.98
    sleep_time = 10 # s

    os.chdir(outpath)
    print(f'Changed dir to {outpath}')

    if restart == 1 and not os.path.exists('opt_nscal.log'):
        restart = 0
    print(f'Restart: {restart}')

    # Create the "global" setup folder in outpath; SSE file is written here.
    # os.system('mkdir setup')
    Path("setup").mkdir(exist_ok=True)
    # One nscal level index per entry in `domain` (domains + interfaces).
    nscal_index_list = [0 for d in domain]
    start_nscal_level = 0
    # Tracks how many stability-probability blocks have been computed in the log.
    n_calc_Q = -1

    #################################################################
    ## If no restart open log file
    ## else read restart infor from log file
    if restart == 0:
        fo = open('opt_nscal.log', 'w')
        fo.write('## Start at %s\n'%time.asctime(time.localtime(time.time())))
        fo.close()
        
    else:
        f = open('opt_nscal.log')
        lines = f.readlines()
        f.close()
        tag_nscal = False
        dom_idx = 0
        for line in lines:
            if line.startswith('## Round '):
                start_nscal_level = int(line.strip().split(':')[0].strip().split()[-1])-1
                dom_idx = 0
            elif line.startswith('## Final nscal:'):
                print('Finished. No need to restart.')
                sys.exit()
            elif line.startswith('-> Set nscal as:'):
                tag_nscal = True
            elif line.startswith('-> Running simulations...'):
                tag_nscal = False
            elif line.startswith('-> Probability of domain stability:'):
                n_calc_Q += 1
            elif tag_nscal:
                nscal = line.strip().split('=')[-1].strip()
                old_nscal = "%.4f"%(nscal_set[domain[dom_idx]['class']][nscal_index_list[dom_idx]])
                if old_nscal != nscal:
                    nscal_index_list[dom_idx] += 1
                new_nscal = "%.4f"%(nscal_set[domain[dom_idx]['class']][nscal_index_list[dom_idx]])
                if new_nscal != nscal:
                    print('Error: nscal mismatch at round #%d, domain #%d level #%d'%(start_nscal_level+1, dom_idx+1, nscal_index_list[dom_idx]+1))
                    sys.exit()
                dom_idx += 1
        f = open('opt_nscal.log', 'w')
        for idx, line in enumerate(lines):
            if line.startswith('## Round %d:'%(start_nscal_level+1)):
                break
            f.write(line)
        if n_calc_Q == start_nscal_level:
            for line in lines[idx+1:]:
                if line.startswith('-> Probability of domain stability:'):
                    break
                f.write(line)
        f.close()
    #################################################################

    print('Cleaning PDB')
    clean_pdb = clean_PDB(input_pdb)

    print('Get secondary structure elements')
    get_secondary_structure(clean_pdb)

    ### get full path to scripts 
    global full_path
    full_path = os.path.realpath(__file__).split('/')[:-1]
    full_path = '/'+os.path.join(*full_path)
    print(f'full_path: {full_path}')

    worker_idx = 0
    for iteration in range(start_nscal_level, len(nscal_set['a'])):
        # Each "round" corresponds to trying a particular nscal level index (0..K-1),
        # with per-domain/interface indices potentially advanced when unstable.
        round_dir = Path(f"round_{iteration+1}")
        round_dir.mkdir(exist_ok=True)
        os.chdir(round_dir)
        
        # os.system('mkdir round_%d'%(iteration+1))
        # os.chdir('round_%d'%(iteration+1))
        
        if n_calc_Q < start_nscal_level or iteration > start_nscal_level:
            Path("setup").mkdir(exist_ok=True)
            # os.system('mkdir setup')
            fo = open('../opt_nscal.log', 'a')
            fo.write("## Round %d:\n"%(iteration+1))
            fo.write("-> Set nscal as:\n")
            for dom_idx, dom in enumerate(domain):
                # Materialize the current nscal value for this domain/interface.
                dom['nscal'] = nscal_set[dom['class']][nscal_index_list[dom_idx]]
                if dom['class'] == 'i':
                    fo.write('   Interface %d|%d: nscal = %.4f\n'%(dom['range'][0]+1, dom['range'][1]+1, dom['nscal']))
                else:
                    fo.write('   Domain %d: nscal = %.4f\n'%(dom_idx+1, dom['nscal']))
            fo.close()
            
            os.chdir('setup')

            print('Create CG model')
            (prefix, prm_name) = creat_CG_model("../../%s"%clean_pdb, casm, domain)
            print('Model sucessfully made')

            ## parse the CHARMM CG forcefeild file into a .xml file so we can use openmm
            print('Parsing CHARMM CG FF files to make OpenMM .xml FF file')
            os.system('parse_cg_prm.py -p %s -t %s'%(prm_name, prefix+'.top'))
            prm_name = prm_name.split('.prm')[0]+'.xml'
            print(f'File sucessfully made: {prm_name}')
            os.chdir('../')
            
            fo = open('../opt_nscal.log', 'a')
            fo.write("-> Running simulations...\n")
            fo.close()
            
            start_time = [time.time() for i in range(ntraj)]
            start_step = [0 for i in range(ntraj)]
            for i in range(ntraj):
                if os.path.exists('%d.out'%(i+1)):
                    # Best-effort attempt to detect restarts by reading last line.
                    # (Note: uses `tail` via shell.)
                    info = os.popen('tail -n 1 %d.out'%(i+1)).readlines()[0]
                    if not info.strip().startswith('Time') and not info.strip() == '':
                        start_step[i] = int(info.strip().split()[1])
            
            # Generate and submit SLURM scripts for all trajectories
            job_ids = []
            slurm_scripts = []
            round_num = iteration + 1  # Round number (1-based)
            
            print(f"Generating {ntraj} SLURM scripts for parallel MD simulations (Round {round_num})...")
            for i in range(ntraj):
                # Random seed for each trajectory (passed to downstream MD runner).
                rand = int(random.random()*1e7)
                MD_script = os.path.join(full_path, 'post_trans_single_run_v3.py')
                MD_cmd = f'python {MD_script} setup/%s.psf setup/%s.cor setup/%s %f %d %d %d %d ../setup/secondary_struc_defs.txt 1.1 setup/%s.cor %d'%(prefix, prefix, prm_name, temperature, ppn, i+1, rand, sim_step, prefix, 1)
                
                # Create SLURM script for this trajectory (include round number in filename)
                script_path = f'slurm_r{round_num}_traj_{i+1}.sh'
                create_slurm_script(i+1, MD_cmd, script_path, job_name_prefix=f"{prefix_name}_r{round_num}", partition=f"{partition_name}")
                slurm_scripts.append(script_path)
                print(f"  Created: {script_path}")
            
            # Submit all SLURM jobs
            print(f"Submitting {ntraj} SLURM jobs...")
            for script_path in slurm_scripts:
                result = subprocess.run(
                    ['sbatch', script_path],
                    capture_output=True,
                    text=True,
                    cwd=os.getcwd()
                )
                if result.returncode == 0:
                    # Extract job ID from sbatch output (format: "Submitted batch job 12345")
                    job_id = result.stdout.strip().split()[-1]
                    job_ids.append(job_id)
                    print(f"  Submitted {script_path}: job ID {job_id}")
                else:
                    print(f"  ERROR: Failed to submit {script_path}")
                    print(f"  Error output: {result.stderr}")
            
            # Wait for all jobs to complete
            if job_ids:
                wait_for_slurm_jobs(job_ids)
            
            # After all MD simulations complete, run Q calculations for all trajectories
            print(f"Running Q calculations for all {ntraj} trajectories...")
            Q_script = os.path.join(full_path, 'calc_native_contact_fraction.pl')
            for i in range(ntraj):
                Q_cmd = f'perl {Q_script} -i setup/%s.cor -d %s -s ../setup/secondary_struc_defs.txt -t %d.dcd -r 1'%(prefix, dom_def_file, i+1)
                print(f"  Running Q calculation for trajectory {i+1}: {Q_cmd}")
                os.system(Q_cmd)
            
            print(f'Simulations finished and Q calculated')
    

        fo = open('../opt_nscal.log', 'a')
        fo.write("-> Probability of domain stability:\n")
        fo.close()
        if_stable = [[] for d in domain]
        stable_frac = [[] for d in domain]
        for i in range(ntraj):
            # Each `qbb_i.dat` should contain per-domain/interface Q values per frame.
            f = open('qbb_%d.dat'%(i+1))
            lines = f.readlines()
            f.close()
            num_fold_list = np.zeros(len(domain))
            for line in lines[1:]:
                qbb_list = line.strip().split()[:-1]
                for j, qbb in enumerate(qbb_list):
                    # Convention: Q == -1 means "not applicable" (treated as folded/stable).
                    if float(qbb) > Q_threshold or float(qbb) == -1:
                        num_fold_list[j] += 1
            frac = num_fold_list/(len(lines)-1)
            for j, f in enumerate(frac):
                stable_frac[j].append(f)
                if f >= frame_threshold:
                    if_stable[j].append(True)
                else:
                    if_stable[j].append(False)
        fo = open('../opt_nscal.log', 'a')
        tag_break = True
        for dom_idx, dom in enumerate(domain):
            if dom['class'] == 'i':
                fo.write('   Interface %d|%d: '%(dom['range'][0]+1, dom['range'][1]+1))
            else:
                fo.write('   Domain %d: '%(dom_idx+1))
            for frac in stable_frac[dom_idx]:
                fo.write('%.3f '%frac)
            if np.all(if_stable[dom_idx]):
                fo.write('stable.\n')
            else:
                fo.write('instable.\n')
                # update nscal index
                # If any trajectory is unstable for this domain/interface, we increase its nscal
                # level for the next round.
                nscal_index_list[dom_idx] += 1
                tag_break = False
        fo.close()
        os.chdir('../')
        if tag_break:
            # All domains/interfaces are stable at the current nscal settings: stop optimization.
            break

    fo = open('opt_nscal.log', 'a')
    fo.write('## Final nscal:\n')
    for dom_idx, dom in enumerate(domain):
        if dom['class'] == 'i':
            fo.write('   Interface %d|%d: nscal = %.4f\n'%(dom['range'][0]+1, dom['range'][1]+1, dom['nscal']))
        else:
            fo.write('   Domain %d: nscal = %.4f\n'%(dom_idx+1, dom['nscal']))
    fo.close()

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    # Print wallclock runtime when executed as a script.
    print(f'NORMAL TERMINATION: {end_time - start_time}')

