# The TOPO runner lives in the package as topo.mdrun. This thin wrapper keeps
# each tutorial self-contained, so `python run_simulation.py -f md.ini` still
# works; it is exactly equivalent to `python -m topo.mdrun -f md.ini`
# (or the installed console command `topo-mdrun -f md.ini`).
from topo.mdrun import mdrun

if __name__ == "__main__":
    mdrun()
