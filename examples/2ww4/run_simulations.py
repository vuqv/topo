import argparse
import configparser
import time
import warnings
from distutils.util import strtobool
from json import loads

import numpy as np
import openmm as mm
from openmm import unit

from parmed.exceptions import OpenMMWarning

import topo


cg_model = topo.models.buildCoarseGrainModel('2ww4.pdb', minimize=False, model='topo', domain_def='domain.yaml', stride_output_file='stride.dat', box_dimension=False)
