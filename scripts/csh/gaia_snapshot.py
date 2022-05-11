#! /usr/bin/env python
#
#  Turn gaia data into an NEMO snapshot sphere , with straight masses as well as
#  densities based on nearby neighbors
#
#  Should work with python-unsiotools, or if not available, uns_tools
#
#  @todo   use velocities as well
#
#  GAIA archive (for manual downloads)      https://gea.esac.esa.int/archive/
#  astroquery (for programmatic downloads)  https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
#
#  7-may-2022 - created on a rainy cold windy Greenman Festival in Greenbelt while playing Mahjong - PJT
#  10-may-2022 - merged in JCL's unsio method - if it fails, the old method will be attempted
#

import os
import sys
import time
import argparse
import numpy as np
from astroquery.gaia import Gaia
try:
    import unsio.output as uns_out 
    import unsiotools.simulations.cfalcon as falcon
    Qtools = False
except:
    print("Warning: unsio not available, trying uns_* tools", file=sys.stderr)
    print("This also forces NEMO tools")
    print("pip install python-unsiotools", file=sys.stderr)
    Qtools = True


# defaults
def_db    = "gaiaedr3.gaia_source"    #   others:  "gaiadr2.gaia_source"
def_pmax  = 25
def_dump  = False
def_dens  = False
def_nemo  = Qtools


def computeDensity(pos):
    # generate fake masses
    mass=np.ones(int(len(pos)/3),dtype=np.float32)
    tic=time.perf_counter()
    cf=falcon.CFalcon()
    ok,rho,hsml=cf.getDensity(pos,mass)
    toc=time.perf_counter()
    print(f"Density computed in {toc-tic:0.4f} seconds",file=sys.stderr)
    return rho,hsml

def saveSnapshot(myoutfile,pos,rho,hsml):
    # we instantiate a CUNS_OUT object
    my_out=uns_out.CUNS_OUT(myoutfile,"nemo") #
    timex=0
    status=my_out.setData(timex,"time")         # set time    
    status=my_out.setData(pos  ,"all","pos" )   # set particles position    
    if rho is not None:
        status=my_out.setData(rho  ,"all","rho" )   # set particles density    
        status=my_out.setData(hsml ,"all","hsml")   # set particles hsml    
    my_out.save()  # write on file system    
    my_out.close() # close


def process(args):
    db    = args.db
    pmax  = args.pmax
    dump  = args.dump
    nemo  = args.nemo

    print("Using %s with parallax > %g mas" % (db,pmax), file=sys.stderr)

    if os.path.exists(args.filename):
        print(f'File {args.filename} already exists', file=sys.stderr)
        sys.exit(0)

    tic=time.perf_counter()
    
    if dump:
        job = Gaia.launch_job_async("SELECT l,b,parallax "
                                    "FROM %s "
                                    "WHERE parallax >= %g" % (db,pmax),
                                    dump_to_file=True, output_format="ecsv")
        f1 = job.outputFile
        print(f"Dumpfile {f1}", file=sys.stderr)
    else:
        job = Gaia.launch_job_async("SELECT l,b,parallax "
                                    "FROM %s "
                                    "WHERE parallax >= %g" % (db,pmax))
        f1 = None

    r = job.get_results()
    toc=time.perf_counter()
    n = len(r)
    print(f'Found {n} stars in {toc-tic:0.4f} seconds',file=sys.stderr)

# nemo dump    Converted 574531 stars in 1.3434 seconds
#      dump    Converted 574531 stars in 0.9047 seconds

    tic=time.perf_counter()
    if args.nemo:
        using = "NEMO"
        if not 'NEMO' in os.environ:
            print("NEMO environment not loaded", file=sys.stderr)
            sys.exit(0)
        print("Converting using NEMO tools")
        if dump:
            expr = "%3*cosd(%1)*cosd(%2),%3*sind(%1)*cosd(%2),%3*sind(%2)"
            # f2 = 'gaia%g.snap'  % pmax
            f2 = args.filename
            cmd = "zcat %s | tabcomment - delete=t | awk -F, '{print $1,$2,1000/$3}' | tabmath - - '%s' all | tabtos - %s block1=x,y,z nbody=%d" % (f1,expr,f2,n)
            print("Dumping from %s to %s" % (f1,f2))
            print(cmd)
            os.system(cmd)
            print("Written",f2)
        else:
            l = r['l']
            b = r['b']
            d = 1000.0/r['parallax']
            rpd = np.pi/180.0
            x = d * np.cos(l*rpd) * np.cos(b*rpd)
            y = d * np.sin(l*rpd) * np.cos(b*rpd)
            z = d * np.sin(b*rpd)
            f1 = 'gaia%g.xyz'   % pmax
            f2 = args.filename            
            np.savetxt(f1,np.c_[x,y,z], fmt="%f")
            if args.dens:
                cmd = "tabtos %s - block1=x,y,z nbody=%d debug=-1 |uns_addmass - - | uns_density - %s all" % (f1,n,f2)
            else:
                cmd = "tabtos %s %s block1=x,y,z nbody=%d debug=-1" % (f1,f2,n)
            print("Converting from %s to %s" % (f1,f2))
            print(cmd)
            os.system(cmd)
            print("Written",f2)
    else:
        using = "UNSIO"
        l = r['l']
        b = r['b']
        d = 1000.0/r['parallax']
        rpd = np.pi/180.0
        x = d * np.cos(l*rpd) * np.cos(b*rpd)
        y = d * np.sin(l*rpd) * np.cos(b*rpd)
        z = d * np.sin(b*rpd)
        # flatten positions to 1D array [x,y,z,x,y,z,x,y,z....]
        pos=np.concatenate((x.data,y.data,z.data))
        pos.shape=(3,len(x.data))
        pos=np.float32(pos.T.reshape(-1))        
        # compute density
        if args.dens:
            rho,hsml=computeDensity(pos)
        else:
            rho,hsml=(None,None)
        # save data using unsio
        saveSnapshot(args.filename,pos,rho,hsml)

    toc=time.perf_counter()
    print(f'Converted {n} stars in {toc-tic:0.4f} seconds using {using}',file=sys.stderr)

def commandLine():
    # help
    parser = argparse.ArgumentParser(description="Download and convert GAIA Data to NEMO Snapshot",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    

    parser.add_argument('filename',    help="output snapshot name",    default=None)
    parser.add_argument('--pmax',      help="pmax value",              default=def_pmax,type=float)
    parser.add_argument('--db',        help="Gaia DB",                 default="gaiaedr3.gaia_source")
    parser.add_argument('--dump',      help='force ECSV dump',         dest="dump",  action="store_true" , default=def_dump)
    parser.add_argument('--density',   help='compute local density',   dest="dens",  action="store_true",  default=def_dens)
    parser.add_argument('--nemo',      help='force NEMO tools [slow]', dest="nemo",  action="store_true",  default=def_nemo)    
    # parse
    args = parser.parse_args()
    # start main function
    process(args)

# -----------------------------------------------------
# main program
if __name__ == '__main__':
  commandLine()
