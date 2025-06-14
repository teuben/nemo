#! /usr/bin/env python
#
#     convert snapshots using AMUSE, it can also create a Plummer sphere from scratch
#
#     benchmark creation of Plummer (remember, premature optimization is the root of all evil)
#     - this code for 10M particles:   36.52user 3.94system 0:40.52elapsed 99%CPU txt
#                                     299.94user 6.16system 5:06.33elapsed 99%CPU dyn
#                                      10.07user 3.18system 0:13.27elapsed 99%CPU hdf5 / amuse
#                                       4.07user 1.96system 0:06.04elapsed 99%CPU -mem-
#     - NEMO's mkplummer:               2.61user 0.66system 0:03.28elapsed 99%CPU

import os
import sys
import time
import argparse
import numpy as np

my_version = "14-jun-2025"

my_help = [f"Version: {my_version}",
           "",
           "Convert snapshots using AMUSE, ",
           "optionally create a Plummer sphere if only output is requested",
           "",
           "The following formats are currently supported:",
           "",           
           "amuse       HDF5 file",
           "amuse-txt   text files with AMUSE header format",
           "csv         comma separated files",
           "dyn         Starlab binary structured file",
           "gadget      Gadget binary data file",
           "hdf5        HDF5 file",
           "nemo        NEMO binary structured file",
           "starlab     Starlab binary structured file",
           "tsf         NEMO binary structured file",
           "txt         text file containing a table of values separated by a predefined",
           "vts         text file containing a table of values separated by a predefined character",
           "vtu         text file containing a table of values separated by a predefined character",
           "",
           "Each format option does have it's own options, YET TO BE IMPLEMENTED",
           "Some formats allow overwrite mode, and don't need the -w flag",
           ]


def commandLine():
    parser = argparse.ArgumentParser(description="\n".join(my_help),
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--ifile',     '-i',  help="input snapshot name",    default=None)
    parser.add_argument('--ofile',     '-o',  help="output snapshot name",   default=None)
    parser.add_argument('--ifmt',      '-I',  help="input file format",      default='amuse')
    parser.add_argument('--ofmt',      '-O',  help='output file format',     default='txt')
    parser.add_argument('--nbody',     '-n',  help="nbody if creating",      default=128)
    parser.add_argument('--seed',      '-s',  help='seed if creating',       default=123)
    parser.add_argument('--overwrite', '-w',  help='force overwrite output', action="store_true")
    parser.add_argument('--time',      '-t',  help='extract time (-1=last)', default=-1)
    
    parser.add_argument('--debug',     '-d',  help='add some debugging',     action="store_true")
    # @todo plummer creation has some options we don't expose here
    # :argument radius_cutoff: Cutoff value for the radius (defaults to 22.8042468)
    # :argument mass_cutoff: Mass percentage inside radius of 1
    args = parser.parse_args()
    if args.debug:
        print(args)
    return args
   
args = commandLine()

nbody = int(args.nbody)
seed  = int(args.seed)
ifile = args.ifile
ofile = args.ofile

ifmt  = args.ifmt
ofmt  = args.ofmt
wmode = args.overwrite
ntime = int(args.time)    # -1=last   0=first


try:
    from amuse.units import nbody_system
    from amuse.ic.plummer import new_plummer_model
    from amuse.io import read_set_from_file
    from amuse.io import write_set_to_file
except:
    print("no AMUSE found in your python environment")
    sys.exit(1)



if ifile is None and ofile is not None:
    # write new
    np.random.seed(seed)
    stars = new_plummer_model(nbody)
    write_set_to_file(stars, ofile, format=ofmt, overwrite=wmode)
    print(f"Wrote Plummer model w/ nbody={nbody} to {ofile} in {ofmt} format")
elif ifile is not None:
    snapshots = read_set_from_file(ifile, format=ifmt)
    itime = 0
    for bodies in list(snapshots.history):
        model_time = bodies.get_timestamp()
        nbody = len(bodies)
        if ofile is not None:
            # re-write old
            write_set_to_file(bodies, ofile, format=ofmt, overwrite=wmode)
            if ntime >= 0 and itime == ntime:
                print(f"Wrote selected snapshot {itime} at {model_time} with {nbody} bodies to {ofile} in {ofmt} format {wmode}")
                break
            else:
                print(f"Processing snapshot number {itime} at {model_time} with {nbody} bodies")
        itime = itime + 1
    if ntime < 0:
        print(f"Last snapshot {model_time} written, use -t to select a number < {itime-1}")
else:
    # no input or output, just make a model in memory
    np.random.seed(seed)
    stars = new_plummer_model(nbody)
    print(f"Created Plummer model w/ nbody={nbody} in memory, no I/O.")
    print("Use -h or --help to get help on various options.")



# close shop but don't report provenance

from amuse.support import literature
literature.TrackLiteratureReferences.suppress_output()
