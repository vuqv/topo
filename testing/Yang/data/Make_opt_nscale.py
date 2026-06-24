import glob

pdbs = glob.glob('/storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/post_rebuilt/*')
for f in pdbs:
    tag = f.split('/')[-1].replace('_rebuilt.pdb', '')
    #print(f, tag)

    domain = f'/storage/group/epo2/default/ims86/git_repos/Failure-to-Form_Native_Entanglements/Simulations_of_Native_Entanglement_Misfolding/Rebuild_AllAtom_structures/data/domains/{tag}.txt'
    cmd = f'python src/data/opt_nscal.py -i {f} -d {domain} -o ../../../../git_slugs/Failure-to-Form_Native_Entanglements_slug/Simulations_of_Native_Entanglement_Misfolding/Coarse_Graining_Model/{tag}/'
    print(cmd)
