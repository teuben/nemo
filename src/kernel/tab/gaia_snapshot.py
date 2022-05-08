#! /usr/bin/env python
#
#  Turn gaia data into an XYZ sphere.
#  @todo   use velocities as well
#
#  https://gea.esac.esa.int/archive/
#  See also examples on https://astroquery.readthedocs.io/en/latest/gaia/gaia.html

import os
import sys
import numpy as np
from astroquery.gaia import Gaia

db   = "gaiadr2.gaia_source"
db   = "gaiaedr3.gaia_source"
pmax = 10
dump = False

print("Using %s with parallax > %g mas" % (db,pmax))

if dump:
    job = Gaia.launch_job_async("SELECT l,b,parallax "
                                "FROM %s "
                                "WHERE parallax >= %g" % (db,pmax),
                                dump_to_file=True, output_format="ecsv")
    f = job.outputFile
    print(f)
else:
    job = Gaia.launch_job_async("SELECT l,b,parallax "
                                "FROM %s "
                                "WHERE parallax >= %g" % (db,pmax))
    f = None

r = job.get_results()
n = len(r)
print('Found',n,'stars')

if 'NEMO' in os.environ:
    if dump:
        expr = "%3*cosd(%1)*cosd(%2),%3*sind(%1)*cosd(%2),%3*sind(%2)"
        f2 = 'gaia%g.snap'  % pmax        
        cmd = "zcat %s | tabcomment - delete=t | awk -F, '{print $1,$2,1000/$3}' | tabmath - - '%s' all | tabtos - %s block1=x,y,z nbody=%d" % (f,expr,f2,n)
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
    print(cmd)
    os.system(cmd)
    print("Written",f2)
else:
    print('No NEMO in your environment. Check out %s' % f)

