/*
 *  preprocessor to run (a series of) GALAXY:
 *
 *  The first argument must be a (non-existent) directory name, in
 *  which a new 'galaxy.dat' file will be written, and in which
 *  galaxy will be run. The 'galaxy' program is assumed to be in the
 *  PATH
 *
 *	24-jun-1997		written
 *	18-jul-2001		default is galaxy.exe
 *      18-mar-2004             finally made it work
 *       7-mar-2006     V2.0    no need for snap2ini
 *                      V2.1    also write output snap
 *      17-sep-2013     V2.2    new run interface
 *      17-dec-2019     V3.0    force this only runs the old galaxy13
 *
 * 
 *  @TODO: this program should run in SinglePrecision, since galaxy does
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <history.h>
#include <unfio.h>
#include <run.h>


#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>

string defv[] = {
    "in=???\n         input snapshot (nemo format)",
    "outdir=???\n     output directory in which GALAXY will be run",
    "scale=15.0\n     number of grid cells per unit length",
    "dt=0.05\n        integration time step",
    "dtout=0.5\n      time between particle outputs",
    "dtlog=0.1\n      time between integral checks",
    "tstop=1.0\n      end time",
    "grid=33,33,33\n  number of grid cells in (x,y,z,)",
    "format=%g\n      format for pos,vel for galaxy.ini",
    "header=\n        If given, use this for unfio header/trailer size",
    "exe=galaxy13\n   name of GALAXY13 executable",
    "VERSION=3.0\n    17-dec-2019 PJT",
    NULL,
};

string usage = "frontend to run Sellwood's galaxy13 code";

string cvsid="$Id$";

int nemo_main()
{
    int i, n, nbody, bits, ndata, count, ngrid[3];
    real scale, dt, dtout, dtlog, tstop, tsnap, mass;
    string exefile = getparam("exe");
    string rundir = getparam("outdir");
    string fmt = getparam("format");
    stream datstr, instr, outstr;
    char dname[256];
    char command[256];
    char fmt6[256];
    float *gdata, *gd;
    Body *bp, *btab = NULL;

    if (hasvalue("header"))
      unfsize(getiparam("header"));

    n = nemoinpi(getparam("grid"),ngrid,3);
    if (n>0 && n<=3) {
        for (i=n; i<3; i++)
            ngrid[i] = ngrid[i-1];
    } else
        error("%d Syntax error: %s (need 1,2 or 3 integers)",
                n,getparam("grid"));
    scale = getdparam("scale");
    dt = getdparam("dt");
    dtout = getdparam("dtout");
    dtlog = getdparam("dtlog");
    tstop = getdparam("tstop");
    sprintf(fmt6,"%s %s %s %s %s %s\n",fmt,fmt,fmt,fmt,fmt,fmt);

    run_mkdir(rundir);


    /* prepare the parameter file for galaxy */

    sprintf(dname,"%s/%s",rundir,"galaxy.dat");
    datstr = stropen(dname,"w");    
    fprintf(datstr,"%d %d %d\n",ngrid[0],ngrid[1],ngrid[2]);
    fprintf(datstr,"%g\n",scale);
    fprintf(datstr,"%g\n",dt);
    fprintf(datstr,"%g\n",dtout);
    fprintf(datstr,"%g\n",dtlog);
    fprintf(datstr,"%g\n",tstop);
    strclose(datstr);

    /* read snapshot and output it in something galaxy understands */

    instr = stropen(getparam("in"),"r");
    get_history(instr);
    if (!get_tag_ok(instr, SnapShotTag)) error("no snapshot");
    get_snap(instr, &btab, &nbody, &tsnap, &bits);
    strclose(instr);

    if ( (bits & PhaseSpaceBit) == 0) error("no phasespace");
    for (bp=btab, mass=0.0;bp<btab+nbody; bp++)
      mass += Mass(bp);
    
    sprintf(dname,"%s/%s",rundir,"galaxy.ini");
    datstr = stropen(dname,"w");
    fprintf(datstr,"%g %g %d\n",tsnap,mass,nbody);
    for (bp=btab;bp<btab+nbody; bp++)
      fprintf(datstr,fmt6,
	      Pos(bp)[0],Pos(bp)[1],Pos(bp)[2],
	      Vel(bp)[0],Vel(bp)[1],Vel(bp)[2]);
    strclose(datstr);

    /* run the fortran program */

    run_cd(rundir);
    if (run_sh(exefile))
      error("Problem executing %s",exefile);

    /* Output data from native galaxy (.res) format to snapshot */
    
    sprintf(dname,"%s","galaxy.res");
    instr = stropen(dname,"r");

    sprintf(dname,"%s","galaxy.snap");
    outstr = stropen(dname,"w");
    put_history(outstr);

    ndata = 7*sizeof(float)*nbody;
    gdata = (float *) allocate(ndata);

    bits = (TimeBit | PhaseSpaceBit | PotentialBit);
    while (1) {
      count = unfread(instr,gdata,ndata);
      if (count <= 0) break;
      printf("Time=%g\n",gdata[0]);
      tsnap = gdata[0];
      count = unfread(instr,gdata,ndata);
      if (count <= 0) error("Problem reading posvelphi from galaxy.res");
      for (bp=btab, gd=gdata;bp<btab+nbody; bp++) {
	Pos(bp)[0] = *gd++;
	Pos(bp)[1] = *gd++;
	Pos(bp)[2] = *gd++;
	Vel(bp)[0] = *gd++;
	Vel(bp)[1] = *gd++;
	Vel(bp)[2] = *gd++;
	Phi(bp)    = *gd++;
      }
      put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }
}


