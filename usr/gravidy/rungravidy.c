/*
 *  rungravidy:   NEMO frontend for gravidy
 *	
 *  15-feb-2019  V0.1  initial draft                  PJT
 *
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <history.h>
#include <filefn.h>
#include <run.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/get_snap.c>


string defv[] = {
    "in=???\n            input file (optional - see nbody= ) ",
    "outdir=???\n     output run directory (required)",

    "eta=0.02\n       Time-step parameter for irregular force polynomial",
    "deltat=0.25\n    Output time interval in units of the crossing time",
    "tcrit=2.0\n      Termination time in units of the crossing time",
    "eps=0.05\n       Potential softening parameter",

    "exe=gravidy\n    Name of the executable",

    "VERSION=0.2\n    18-feb-2019 PJT",
    NULL,
};

string usage="front end for gravidy N-body code";


void nemo_main(void)
{
  int nbody;
    real eta, deltat, tcrit, qe, eps, tcomp;
    string exefile = getparam("exe");
    string tabfile = "gravidy.in";
    string rundir = getparam("outdir");
    string infile, fname;
    char dname[256], runcmd[256];
    stream datstr, histr;

    eta = getdparam("eta");
    deltat = getdparam("deltat");
    tcrit = getdparam("tcrit");
    eps = getdparam("eps");

    infile = getparam("in");
    fname = fullname(infile);
    dprintf(0,"fullname: %s\n",fname);

    run_mkdir(rundir);

    stream instr = stropen(infile,"r");
    nbody = get_snap_nbody(instr);
    strclose(instr);
    dprintf(0,"Grabbing nbody=%d\n",nbody);
    sprintf(runcmd,"snapprint %s i,m,x,y,z,vx,vy,vz format=%%.15f > %s/%s", infile,rundir,tabfile);
    run_sh(runcmd);
    

    run_cd(rundir);
    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);


    sprintf(runcmd,"%s -i %s -t %g -o %s -s %g -e %g -z %g -l",exefile,tabfile,tcrit,rundir,eps,eta,deltat);
    run_sh(runcmd);

    sprintf(runcmd,"cat %s.out.snapshot_* | sed 's/# Time://' | tabtos - - time skip,m,pos,vel nbody=%d | csf - OUT3.snap SnapShot",
	    rundir, nbody);
    run_sh(runcmd);    

}
