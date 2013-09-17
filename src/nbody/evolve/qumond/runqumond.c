/*
 *  RUNQUMOND:   run the QuMond program
 *
 *    assemble parameters
 *    needs qmics.dat
 *    create directory, run program
 *
 *    28-apr-2011   0.3   ?
 *    17-sep-2013   0.4  new run interface
 */

#include <nemo.h>

string defv[] = {
  "in=???\n              Input snapshot",
  "outdir=???\n          Run directory",
  "out=\n                Output snapshot(s) if to convert back to NEMO?",
  "numbs=8\n             Number of integration steps",
  "aexpn=0.00544081442\n Initial scale factor (z=1/aexpn-1)",
  "adiv=0.002\n          Normalisation of time steps",
  "om0=0.2623\n          Omega_0 (Omega_0=Omega_cdm+Omega_nu+Omega_baryon)",
  "vsca=512.0\n          length of box (Mpc/h)",
  "hubble=73.20\n        Hubble param used in cosmics",
  "mond=2\n              Desired nu-function (0:no MOND, 1:some weird fcn, 2:simple nu fcn",
  "freq=8\n              frequency to output a binary file",
  "brand=0\n             Restarting option (0=new and ascii, 2=old and binary)",
  "au0=1.0\n             a_0 is empirically ~ 1.2e-8m/s^2. Want to rescale it by a factor s.t. g_0=factor*a_o",
  "VERSION=0.4\n         17-sep-2013 PJT",
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
  stream parstr, outstr;
  string infile = getparam("in");
  string outfile = getparam("out");
  string outdir = getparam("outdir");
  string exefile = "QMcode";
  string datfile = "qmics.dat";
  string parfile = "start.txt";
  string logfile = "QMcode.log";
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

  
  run_mkdir(outdir);
  sprintf(dname,"%s/%s",outdir,parfile);

  parstr = stropen(dname,"w");

  fprintf(parstr,"%d\n",getiparam("numbs"));
  fprintf(parstr,"%s\n",getparam("aexpn"));
  fprintf(parstr,"%s\n",getparam("adiv"));
  fprintf(parstr,"%s\n",getparam("om0"));
  fprintf(parstr,"%s\n",getparam("vsca"));
  fprintf(parstr,"%s\n",getparam("hubble"));
  fprintf(parstr,"%d\n",getiparam("mond"));
  fprintf(parstr,"%d\n",getiparam("freq"));
  fprintf(parstr,"%d\n",getiparam("brand"));
  fprintf(parstr,"%s\n",getparam("au0"));
  strclose(parstr);

  sprintf(cmd,"cp %s %s/%s", infile, outdir, datfile);
  run_sh(cmd);

  run_cd(outdir);
  sprintf(cmd,"%s < %s > %s ", exefile, parfile, logfile);
  run_sh(cmd);

  if (hasvalue("out")) {
    warning("no postprocessing yet - only history written");
  } else {
    outstr = stropen("history.bin","w");
    put_history(outstr);
    strclose(outstr);
  }

}

