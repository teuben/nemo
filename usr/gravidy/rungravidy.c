/*
 *  rungravidy:   NEMO frontend for gravidy
 *	
 *  15-feb-2019  V0.1  initial draft                                    PJT
 *  14-jun-2023  V0.3  reminders on original flags, added some more     PJT
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

// see also https://gravidy.xyz/include/use.html#command-line-options

string defv[] = {
    "in=???\n         Input file, required   (-i --input) ",
    "outdir=???\n     Output run directory, required (-o --output)",

    "eta=0.01\n       Time-step parameter for irregular force polynomial (-e --eta)",
    "deltat=0.125\n   Output time interval in units of the crossing time) (-z --interval)",
    "tcrit=1\n        Termination time in units of the crossing time (-t --time)",
    "eps=0.0001\n     Gravitational softening parameter (-s --softening)",

    "lagrange=t\n     Print information of the Lagrange Radii in every integration time (-l --lagrange)",
    "all=f\n          Print all the information of N-particles in every integration time (-a --all)",
    "gpu=0\n          GPUs to use, by default is the maximum available devices (use even numbers) (-g --gpu)",

    "dryrun=f\n       Dryrun?",
    "exe=gravidy\n    Name of the executable",

    "VERSION=0.4\n    16-jun-2023 PJT",
    NULL,
};

string usage="front end for gravidy N-body code";


void nemo_main(void)
{
    int nbody;
    real eta, deltat, tcrit, eps;
    string exefile = getparam("exe");
    string tabfile = "gravidy.in";
    string rundir = getparam("outdir");
    string infile, fname;
    bool Qdry = getbparam("dryrun");
    bool Qlag = getbparam("lagrange");
    bool Qall = getbparam("all");
    int ngpu = getiparam("gpu");
    char runcmd[256];
    stream histr;

    eta = getdparam("eta");
    deltat = getdparam("deltat");
    tcrit = getdparam("tcrit");
    eps = getdparam("eps");

    infile = getparam("in");
    fname = fullname(infile);
    dprintf(0,"fullname: %s\n",fname);

    if (!Qdry)
      run_mkdir(rundir);

    stream instr = stropen(infile,"r");
    nbody = get_snap_nbody(instr);
    strclose(instr);
    dprintf(0,"Grabbing nbody=%d\n",nbody);
    
    sprintf(runcmd,"snapprint %s i,m,x,y,z,vx,vy,vz format=%%.15f > %s/%s", infile,rundir,tabfile);
    if (Qdry)
      printf("%s\n",runcmd);
    else {
      dprintf(0,"%s\n",runcmd);
      run_sh(runcmd);
    }

    if (!Qdry) {
      run_cd(rundir);
      histr = stropen("history","w");
      put_history(histr);
      strclose(histr);
    }

    if (ngpu > 0) warning("gpu=%d not used yet", ngpu);


    sprintf(runcmd,"%s -i %s -t %g -o %s -s %g -e %g -z %g %s %s",
	    exefile,tabfile,tcrit,rundir,eps,eta,deltat,
	    Qlag ? "-l" : "",
	    Qall ? "-a" : "");
    if (Qdry)
      printf("%s\n",runcmd);
    else {
      dprintf(0,"%s\n",runcmd);      
      run_sh(runcmd);
    }

    sprintf(runcmd,"cat %s.out.snapshot_* | sed 's/# Time://' | tabtos - - time skip,m,pos,vel nbody=%d | csf - OUT3.snap SnapShot",
	    rundir, nbody);
    
    if (Qdry)
      printf("%s\n",runcmd);
    else {
      dprintf(0,"%s\n",runcmd);      
      run_sh(runcmd);
    }

}
