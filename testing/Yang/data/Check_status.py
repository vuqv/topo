import glob
import pandas as pd
import MDAnalysis as mda
import os
import shutil

def list_files_in_highest_numbered_round_folder(directory):
    # Step 1: Get all folder names in the directory
    folders = [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

    # Step 2: Filter folders that start with 'round_' and extract the number after 'round_'
    round_folders = []
    for folder in folders:
        if folder.startswith('round_'):
            try:
                # Extract the number part after 'round_'
                round_number = int(folder.split('_')[1])
                round_folders.append((round_number, folder))  # Store as (number, folder_name)
            except (IndexError, ValueError):
                pass  # Skip if folder name is not in the expected format

    # Step 3: Find the folder with the highest number
    if not round_folders:
        print("No 'round_' folders found.")
        return []

    highest_numbered_folder = max(round_folders, key=lambda x: x[0])[1]
    print(f'highest_numbered_folder: {highest_numbered_folder}')

    # Step 4: List the files in the highest numbered folder
    highest_folder_path = os.path.join(directory, f'{highest_numbered_folder}/setup/')
    print(f'highest_folder_path: {highest_folder_path}')
    files = os.listdir(highest_folder_path)

    return files, highest_folder_path


Finished_count = 0
cmds = [x.split()[-1] for x in open('src/command_files/opt_nscale.cmds').readlines()]
check = {'index':[], 'gene':[], 'pdb':[], 'chain':[], 'Length':[], 'DomDef':[], 'status':[], 'nscales':[]}
for i, c in enumerate(cmds):
    tag = c.split('/')[-2]
    gene, pdb, chain = tag.split('_')
    print(i, tag, gene, pdb, chain)

    # get domain data
    dom_file = f'../Rebuild_AllAtom_structures/data/domains/{tag}.txt'
    print(dom_file)
    dom_data = '\n'.join([x.strip('\n') for x in open(dom_file, 'r').readlines()])
    #print(dom_data)

    # get length
    pdb_file = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/PDBs/{tag}_rebuilt.pdb'
    print(pdb_file)
    u = mda.Universe(pdb_file)
    #print(u)
    CAs = u.select_atoms('name CA')
    #print(CAs, len(CAs))
    Length = len(CAs)

    opt_path = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/opt_nscal.log'
    #print(opt_path)

    # read in the opt_nscal.log and determine if there is a ## Final nscal: line indicating it finished
    opt_lines = [x for x in open(opt_path, 'r').readlines()]
    #print(opt_lines)

    Finished = False
    Finished_idx = 0
    for idx, line in enumerate(opt_lines):
        if '## Final nscal:' in line:
            Finished = True
            Finished_count += 1
            Finished_idx = idx

    ## If finished get the final nscale values from the opt_nscale.log file
    # also move the CG files for that round to the Temp_Quench_Dynamics/{tag}/setup directory
    if Finished:
        print(f'Finished: {gene} {pdb} {chain}')
        final_nscal_lines = '\n'.join([x.replace('   ',' ').strip('\n') for x in opt_lines[Finished_idx + 1 :]])
        
        # make new directory if it doesnt exist
        TQ_dir = os.path.join('../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Temp_Quench_Dynamics/', f'{tag}/setup/')
        if not os.path.exists(TQ_dir):
            os.makedirs(TQ_dir)
            print(f'Made directory: {TQ_dir}')   

        # move datafiles from final round to Temp_Quench_Dynamics directory and make setup directory
        source_directory = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/'
        files, highest_folder_path = list_files_in_highest_numbered_round_folder(source_directory)
        file_pattern = f'{highest_folder_path}/{tag}_*'
        files_to_copy = glob.glob(file_pattern)

        # Copy each file to the destination directory
        for file_path in files_to_copy:
            shutil.copy(file_path, TQ_dir)
            print(f"Copied {file_path} to {TQ_dir}")

        # Find the secondary structure file and copy it
        sec_elems_file = f'../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/setup/secondary_struc_defs.txt'
        shutil.copy(sec_elems_file, TQ_dir)
        print(f'Copied {sec_elems_file} to {TQ_dir}')
       
    else:
        final_nscal_lines = ''
    print(f'{i} {gene} {pdb} {chain} Finished: {Finished} {final_nscal_lines}')
    check['index'] += [i]
    check['gene'] += [gene]
    check['pdb'] += [pdb]
    check['chain'] += [chain]
    check['Length'] += [Length]
    check['DomDef'] += [dom_data]
    check['status'] += [Finished]
    check['nscales'] += [final_nscal_lines]
check = pd.DataFrame(check)
print(check.to_string())
check.to_csv('data/Check.csv', index=False)
print(f'SAVED: data/Check.csv')
print(f'Finished_count: {Finished_count}')