# Tasks
- [x] change main function in mdrun.py to mdrun
- [x] can the dynamic module change to mdrun module? mdrun module still exposed to the bash, topo-mdrun -f md.ini
- [x] the temperature quenching or equilibrium (conventional mdrun) is handle by option in md.ini, has option: e.g: 
- [x] topo-simulation and dynamics.py seem old, remove if safe
- [x] optimization:
  - [x] config file has optional option: min_contacts [default is 0], this set the contact threshold to run the optimize the strength of domain/interface.
  - [x] for example, two domains has very few contact, meaning their interaction is not strong and hance can be consider as do not folded/stable so do not optimize their strength, just set to the first level of interface and not optimize. Same for intra-domain.
- [x] run_simulation.py is now not needed, is this line from opimization needed:
SCRIPT_DIR = Path(__file__).resolve().parent  # for the run_simulation/split scripts
- [x] Revise the optimization to production code, Turn it into optimize module and expose to system at: topo-optimize
- [ ] For multichain simulations, it is good to always separate chain to add appendix _1..N to the trajectory.dcd
- [ ] docs still need to be revised for comprehensive
- [ ] interacting chains
- [x] Remove shift initial coordinate (enigne), this may affect multichain/translation simulation
- [x] anneal_steps is separated from mdstep. the two process write two dcd file.