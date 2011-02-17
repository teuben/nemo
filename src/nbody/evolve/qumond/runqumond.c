/*
 *  RUNQUMOND:   run the QuMond program
 */

#include <nemo.h>

string defv[] = {
  "in=???\n              Input snapshot",
  "outdir=???\n          Run directory",
  "out=\n                Output snapshot(s) if to convert back to NEMO?",
  "numbs=16\n            Number of integration steps",
  "aexpn=0.00544081442\n Initial scale factor (z=1/aexpn-1)",
  "adiv=0.002\n          Normalisation of time steps",
  "om0=0.2623\n          Omega_0 (Omega_0=Omega_cdm+Omega_nu+Omega_baryon)",
  "vsca=512.0\n          length of box (Mpc/h)",
  "hubble=73.20\n        Hubble param used in cosmics",
  "mond=2\n              Desired nu-function (0:no MOND, 1:some weird fcn, 2:simple nu fcn",
  "freq=8\n              frequency to output a binary file",
  "brand=0\n             Restarting option (0=new and ascii, 2=old and binary)",
  "au0=1.0\n             a_0 is empirically ~ 1.2e-8m/s^2. Want to rescale it by a factor s.t. g_0=factor*a_o",
  "VERSION=0.1\n         17-feb-2011 PJT",
  NULL,
};

string usage="Frontend to the QuMond program";

string cvsid="$Id$";

void nemo_main()
{
  int nr, nth, nph, lmax ;
  int model,mond_ind;
  real spl_ord, rmap  ;
  real a0,  dt_iter,scale;
  int iter_max;
  real tol;
  int new, id_new, nout,iene, mrate;
  real tmax, dt_min, cfl;
  int lp_ord;
  stream parstr;
  string infile = getparam("in");
  string outfile = getparam("out");
  string outdir = getparam("outdir");
  string exefile = "QuMond";
  string parfile = "start.txt";
  char dname[256], cmd[256];

  real  aexpn = getdparam("aexpn");
  real  adiv = getdparam("adiv");
  real  om0 = getdparam("om0");
  real  vsca = getdparam("vsca");
  real  hubble = getdparam("hubble");
  int   mond = getiparam("mond");
  int   freq = getiparam("freq");
  int   brand = getiparam("brand");
  real  au0 = getdparam("au0");

  printf("Initial a=%g\n",aexpn);
  printf("Initial z=%g\n",1/aexpn-1);
  printf("Initial timestep = %g\n",aexpn*adiv);


  error("code not ready to run yet");
  
  make_rundir(outdir);
  sprintf(dname,"%s/%s",outdir,parfile);

  parstr = stropen(dname,"w");

  fprintf(parstr,"nr, nth, nph, lmax \n");
  fprintf(parstr,"%d,  %d,  %d,  %d \n",nr, nth, nph, lmax);

  fprintf(parstr,"model,mond_ind,spl_ord, rmap\n");
  fprintf(parstr,"%d, %d, %g, %g\n",model,mond_ind,spl_ord, rmap);

  fprintf(parstr,"a0,  dt_iter,scale\n");
  fprintf(parstr,"%g, %g, %g\n",a0,  dt_iter,scale);

  fprintf(parstr,"iter_max, tol\n");
  fprintf(parstr,"%d, %g\n",iter_max, tol);

  fprintf(parstr,"new, id_new, nout,iene,mrate\n");
  fprintf(parstr,"%d, %d, %d, %d, %d\n",new, id_new, nout,iene,mrate);
  fprintf(parstr,"tmax, dt_min, cfl, lp_ord\n");
  fprintf(parstr,"%g, %g, %g, %d\n",tmax, dt_min, cfl, lp_ord);
  fprintf(parstr,"! --------------\n");
  fprintf(parstr,"! this file was produced automatically by NEMO::runmond\n");
  strclose(parstr);

  sprintf(cmd,"snapmody in=%s out=%s/%s",infile,outdir,"mout00.bin");
  run_program(cmd);

  goto_rundir(outdir);
  sprintf(cmd,"rmond");
  run_program(cmd);

  if (hasvalue("out")) {
    sprintf(cmd,"cat mout??.bin | modysnap - %s",outfile);
    run_program(cmd);
  }

}

goto_rundir(string name)
{
  dprintf(0,"CHDIR: %s\n",name);
  if (chdir(name))
    error("Cannot change directory to %s",name);
}

make_rundir(string name)
{
  dprintf(0,"MKDIR: %s\n",name);
  if (mkdir(name, 0755))
    warning("Run directory %s already exists",name);
}

run_program(string cmd)
{
  dprintf(0,"EXEC: %s\n",cmd);
  system(cmd);
}


