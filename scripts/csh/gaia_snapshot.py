#! /usr/bin/env python
#
#  Turn gaia data into an NEMO snapshot sphere , with straight masses as well as
#  densities based on nearby neighbors
#
#  @todo   use velocities as well
#
#  GAIA archive (for manual downloads)      https://gea.esac.esa.int/archive/
#  astroquery (for programmatic downloads)  https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
#
#  7-may-2022 - created on a rainy cold windy Greenman Festival in Greenbelt while playing Mahjong - PJT
#

import os
import sys
import numpy as np
from astroquery.gaia import Gaia

db   = "gaiadr2.gaia_source"
db   = "gaiaedr3.gaia_source"
pmax = 4
dump = False

print("Using %s with parallax > %g mas" % (db,pmax))

if dump:
    job = Gaia.launch_job_async("SELECT l,b,parallax "
                                "FROM %s "
                                "WHERE parallax >= %g" % (db,pmax),
                                dump_to_file=True, output_format="ecsv")
    f1 = job.outputFile
else:
    job = Gaia.launch_job_async("SELECT l,b,parallax "
                                "FROM %s "
                                "WHERE parallax >= %g" % (db,pmax))
    f1 = None

r = job.get_results()
n = len(r)
print('Found',n,'stars')

if 'NEMO' in os.environ:
    print("Converting using NEMO")
    if f1 != None:
        expr = "%3*cosd(%1)*cosd(%2),%3*sind(%1)*cosd(%2),%3*sind(%2)"
        f2 = 'gaia%g.snap'  % pmax        
        cmd = "zcat %s | tabcomment - delete=t | awk -F, '{print $1,$2,1000/$3}' | tabmath - - '%s' all | tabtos - %s block1=x,y,z nbody=%d" % (f1,expr,f2,n)
        print("Dumping from %s to %s" % (f1,f2))
    else:
        l = r['l']
        b = r['b']
        d = 1000.0/r['parallax']
        rpd = np.pi/180.0
        x = d * np.cos(l*rpd) * np.cos(b*rpd)
        y = d * np.sin(l*rpd) * np.cos(b*rpd)
        z = d * np.sin(b*rpd)
        f1 = 'gaia%g.xyz'   % pmax
        f2 = 'gaia%g.snap'  % pmax
        np.savetxt(f1,np.c_[x,y,z], fmt="%f")
        cmd = "tabtos %s %s block1=x,y,z nbody=%d debug=-1" % (f1,f2,n)
        print("Converting from %s to %s" % (f1,f2))
    print(cmd)
    os.system(cmd)
    print("Written",f2)
    f3 = 'gaia%g_rho.snap' % pmax
    cmd = 'uns_addmass %s - | uns_density - %s all' % (f2,f3)
    os.system(cmd)
    print("Written",f3)
else:
    print('No NEMO in your environment. Check out %s' % f)

