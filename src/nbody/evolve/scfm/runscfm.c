/*
 *  MAIN.C: NEMO driver for SCFM
 *		this program also needs 'tabtos' in $NEMOBIN
 *	
 *      22-jan-98   V1.0    cloned off runbody2                 pjt
 *	16-dec-99   V1.1    exe is now scfm.exe within NEMO     pjt
 *			    back to non-exe
 *	 3-apr-01   V1.2    back to exe .... why on earth....   pjt
 *      17-mar-06   V1.3    using fullname()  for in=           pjt
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filefn.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

string defv[] = {
    "in=???\n         input file (standard snapshot)",
    "outdir=???\n     output run directory (required)",
    "nsteps=100\n     Number of timesteps  ",
    "noutbod=100\n    Output system state (on SNAPxxx)",
    "noutlog=10\n     Output logfile data (on SCFLOG)",
    "dtime=0.01\n     Integration timestep ",
    "gravconst=1.0\n  Gravitional Constant ",
    "selfgrav=t\n     include selfgravity?",
    "inptcoef=f\n     ",
    "outpcoef=f\n     ",
    "zeroodd=f\n      ",
    "zeroeven=f\n     ",
    "fixacc=f\n       ",
    "headline=\n      run comment",
    "VERSION=1.3\n    17-mar-06 PJT",
    NULL,
};

string usage="Self Consistent Field Method N-body code";

string cvsid="$Id$";

nemo_main()
{
    string exefile = "scfm.exe";
    string parfile = "SCFPAR";
    string rundir = getparam("outdir");
    string infile, fname;
    char dname[256], runcmd[256];
    string headline = getparam("headline");
    stream datstr, histr;
    int nsteps = getiparam("nsteps");
    int noutbod = getiparam("noutbod");
    int noutlog = getiparam("noutlog");
    real dtime = getrparam("dtime");
    real gravconst = getrparam("gravconst");
    bool selfgrav = getbparam("selfgrav");
    bool inptcoef = getbparam("inptcoef");
    bool outpcoef = getbparam("outpcoef");
    bool zeroodd  = getbparam("zeroodd");
    bool zeroeven  = getbparam("zeroeven");
    bool fixacc = getbparam("fixacc");

    make_rundir(rundir);

    sprintf(dname,"%s/%s",rundir,parfile);
    datstr = stropen(dname,"w");    


    fprintf(datstr,"C Basic input parameters\n");
    fprintf(datstr,"%s\n",headline);
    fprintf(datstr," %d\n",nsteps);
    fprintf(datstr," %d\n",noutbod);
    fprintf(datstr," %d\n",noutlog);
    fprintf(datstr," %g\n",dtime);
    fprintf(datstr," %g\n",gravconst);
    fprintf(datstr,"%s\n", selfgrav ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"%s\n", inptcoef ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"%s\n", outpcoef ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"%s\n", zeroodd  ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"%s\n", zeroeven ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"%s\n", fixacc   ? ".TRUE." : ".FALSE.");
    fprintf(datstr,"C end of parameter input\n");
    fclose(datstr);

    goto_rundir(rundir);
    histr = stropen("history","w");
    put_history(histr);
    strclose(histr);

    if (hasvalue("in")) {
	infile = getparam("in");
	if (*infile == '-') {		/* do something special for pipes */
	  sprintf(runcmd,"snapprint - m,x,y,z,vx,vy,vz header=t > SCFBI");
	} else if (*infile == '/') {	/* regular file */
	  fname = fullname(infile);
	  sprintf(runcmd,"snapprint %s m,x,y,z,vx,vy,vz header=t > SCFBI",fname);
        } 
        dprintf(0,"%s\n",runcmd);
        if (system(runcmd)) error("Error converting input data");
    }
#if 0
    sprintf(runcmd,"cp ../SCFPAR .; %s",exefile);
#else
    sprintf(runcmd,"%s",exefile);
#endif
    run_program(runcmd);

    sprintf(runcmd,
        "cat SNAP??? | tabtos - scfm.dat nbody,time mass,pos,vel,phi");
    run_program(runcmd);

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
