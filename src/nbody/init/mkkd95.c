/*
 *  mkkd95:   make a Kuijken-Dubinski (1995) compostite B-D-H model
 *       This is simply a driver that runs the GalactICS routines
 *	
 *  6-mar-04    Created              PJT
 *
 *
 */

#include <stdinc.h>
#include <getparam.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
  "out=???\n        output snapshot (a rundirectory $out.tmpdir is also created)",
  "nbulge=4000\n    Number of particles in bulge",
  "ndisk=8000\n     Number of particles in disk",
  "nhalo=6000\n     Number of particles in halo",

  "psi0=-4.6\n       in.dbh: Psi0    (HALO)",
  "v0=1.42\n         in.dbh: v0",
  "q=1\n             in.dbh: q",
  "rck2=0.1\n        in.dbh: (rc/rk)^2",
  "ra=0.8\n          in.dbh: ra",

  "md=0.867\n        in.dbh: M_d      (DISK)",
  "rd=1\n            in.dbh: R_d",
  "router=5\n        in.dbh: R_outer",
  "zd=0.1\n          in.dbh: z_d",
  "drtrunc=0.5\n     in.dbh: dR_trunc",

  "rhob=14.45\n      in.dbh: rho_b    (BULGE)",
  "psicut=-2.3\n     in.dbh: psi_cut",
  "sigb=0.714\n      in.dbh: sib_b",

  "deltar=0.01\n     in.dbh: delta_r",
  "nr=2400\n         in.dbh: nr",
  "nharm=10\n        in.dbh: number of harmonics",

  "sigvr0=0.47\n     in.diskdf: central radial velocity dispersion",
  "sigr0=1.0\n       in.diskdf: scalelength of sig_r^2",
  
  "ncorr=50\n        in.diskdf: number of intervals for correction function",
  "niter=10\n        in.diskdf: number of iterations",
  
  "fstreamb=0.75\n   in.bulge:",
  "fstreamh=0.5\n    in.halo:",
  
  "seed=0\n          in.- Random Seed ",
  "zerocm=t\n        in.- Center the snapshot?",

  "bin=.\n           directory in which KD95 binaries live",

  "VERSION=1.0\n     6-mar-04 PJT",
  NULL,
};

string usage="Kuijken-Dubinsky-95 composite bulge-disk-halo model";


nemo_main()
{
    int nbody, nfix, nrand, nrun, kstart;
    real eta, deltat, tcrit, qe, eps, tcomp;
    real alphas, body1, bodyn;
    real q, vxrot, vzrot, rbar, zmbar;
    string out=getparam("out");
    string infile, kd_bindir = getparam("bin");
    char fullname[256], runcmd[256], rundir[256];
    stream datstr, histr;
    int iseed = getiparam("seed");
    int zerocm = getbparam("zerocm") ? 1 : 0;

    int nbulge = getiparam("nbulge");
    int ndisk  = getiparam("ndisk");
    int nhalo  = getiparam("nhalo");

    sprintf(rundir,"%s.tmpdir",out);

    make_rundir(rundir);
    goto_rundir(rundir);


    datstr = stropen("in.dbh","w"); 
    fprintf(datstr,"%s\n",  
	    nhalo > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g %g %g\n",
	    getdparam("psi0"),
	    getdparam("v0"),
	    getdparam("q"),
	    getdparam("rck2"),
	    getdparam("ra"));

    fprintf(datstr,"%s\n",  
	    ndisk > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g %g %g\n",
	    getdparam("md"),
	    getdparam("rd"),
	    getdparam("router"),
	    getdparam("zd"),
	    getdparam("drtrunc"));

    fprintf(datstr,"%s\n",  
	    nbulge > 0 ? "y" : "n");
    fprintf(datstr,"%g %g %g\n",
	    getdparam("rhob"),
	    getdparam("psicut"),
	    getdparam("sigb"));
    fprintf(datstr,"%g %d\n",
	    getdparam("deltar"),
	    getiparam("nr"));
    fprintf(datstr,"%d\n",
	    getiparam("nharm"));
    fprintf(datstr,"%s\n",
	    "dbh.ps/ps");
    strclose(datstr);


    datstr = stropen("in.diskdf","w"); 
    fprintf(datstr,"%g %g\n",
	    getdparam("sigvr0"),
	    getdparam("sigr0"));
    fprintf(datstr,"%d\n",
	    getiparam("ncorr"));
    fprintf(datstr,"%d\n",
	    getiparam("niter"));
    fprintf(datstr,"%s\n",
	    "diskdf.ps/ps");
    strclose(datstr);

    datstr = stropen("in.bulge","w"); 
    fprintf(datstr,"%g\n",
	    getdparam("fstreamb"));
    fprintf(datstr,"%d\n",
	    getiparam("nbulge"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);

    datstr = stropen("in.disk","w"); 
    fprintf(datstr,"%d\n",
	    getiparam("ndisk"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);

    datstr = stropen("in.halo","w"); 
    fprintf(datstr,"%g\n",
	    getdparam("fstreamh"));
    fprintf(datstr,"%d\n",
	    getiparam("nhalo"));
    fprintf(datstr,"%d\n",
	    iseed);
    fprintf(datstr,"%d\n",
	    zerocm);
    fprintf(datstr,"%s\n",
	    "dbh.dat");
    strclose(datstr);


    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);

    datstr = stropen("make-it","w"); 
    fprintf(datstr,"#! /bin/sh\n");
    fprintf(datstr,"# created by NEMO\n");
    fprintf(datstr,"%s/dbh < in.dbh\n",kd_bindir);
    fprintf(datstr,"%s/getfreqs\n",kd_bindir);
    fprintf(datstr,"%s/diskdf < in.diskdf\n",kd_bindir);
    fprintf(datstr,"%s/genbulge < in.bulge  > bulge\n",kd_bindir);
    fprintf(datstr,"%s/gendisk  < in.disk   > disk\n",kd_bindir);
    fprintf(datstr,"%s/genhalo  < in.halo   > halo\n",kd_bindir);
    fprintf(datstr,"%s/mergerv disk bulge halo > galaxy\n");
    fprintf(datstr,"tabtos galaxy ../%s nbody,time mass,pos,vel\n",out);
    strclose(datstr);

    run_program("chmod +x make-it; ./make-it");

}


goto_rundir(string name)
{
    if (chdir(name))
        error("Cannot change directory to %s",name);
}

make_rundir(string name)
{
    if (mkdir(name, 0755))
        warning("Run directory %s already exists",name);
}

run_program(string cmd)
{
    system(cmd);
}


/*
 *	Order of input lines in "nbody1.in" for a new run (KSTART=1)
 *
 *          variables                   condition               where
 *  ----------------------------    -------------------
 *  KSTART TCOMP                                            nbody1.f    MAIN
 *  Nbody NFIX NRAND NRUN                                   input.f     INPUT
 *  ETA DELTAT TCRIT QE EPS                                 input.f     INPUT
 *  KZ(1..15)                                               input.f     INPUT
 *      ALPHAS BODY1 BODYN          KZ(4).NE.2              data.f      DATA
 *      SEMI ECC                    KZ(12).NE.0             data.f      DATA
 *  Q VXROT VZROT RBAR ZMBAR                                scale.f     SCALE
 *      NFRAME DELTAF               KZ(7).GT.0              scale.f     SCALE
 *      XCM ECC                     KZ(8).GT.0              subsys.f    SUBSYS
 *
 */
 
