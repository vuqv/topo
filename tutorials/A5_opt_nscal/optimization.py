# The TOPO nscale optimizer lives in the package as topo.optimize. This thin
# wrapper keeps the tutorial self-contained, so `python optimization.py -f
# optimize.ini` still works; it is exactly equivalent to
# `python -m topo.optimize -f optimize.ini` (or the installed console command
# `topo-optimize -f optimize.ini`).
from topo.optimize import optimize

if __name__ == "__main__":
    optimize()
