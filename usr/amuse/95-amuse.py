# ~/.ipython/profile_amuse/startup
import math
import numpy as np
import matplotlib.pyplot as plt
#
import amuse
from amuse.io import read_set_from_file
from amuse.io import write_set_to_file
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
#
from amuse.lab import *
#
_aversion = "TBD"
print(f"amuse Version {_aversion} loaded from {amuse.__file__}")
