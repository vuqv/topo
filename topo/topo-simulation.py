#!/usr/bin/env python

"""
This script is a command-line program for running a simulation of a protein system using the
TOPO coarse-grained force field. The script uses the argparse library to parse a simulation config file
specified by the user as a command-line argument. The script then imports the topo library,
and creates an instance of the Dynamics class, passing in the parsed config file.
Finally, the script runs the simulation by calling the run() method on the Dynamics class instance.
"""

import argparse

import topo

# Parse config file:
parser = argparse.ArgumentParser(
    description="\n Usage: python run_simulation.py -f md.ini \n or topo-simulation -f md.ini")
parser.add_argument('-input', '-f', type=str, help='simulation config file')
args = parser.parse_args()

topo_sim = topo.dynamics.Dynamics(args.input)
topo_sim.run()
