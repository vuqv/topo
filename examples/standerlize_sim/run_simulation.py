# Import OpenMM library
# Import hpsOpenMM library
import argparse

import topo

# Parse config file:
parser = argparse.ArgumentParser(
    description="\n Usage: python run_simulation.py -f md.ini \n or topo-simulation -f md.ini")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

hps_sim = topo.dynamics.Dynamics(args.input)
hps_sim.run()
