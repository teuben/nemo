/*
 *  RUNMOND:   run the mond program
 * 
 *  2009 ?
 *  17-sep-2013  0.3   new run interface      PJT
 */

#include <nemo.h>

string defv[] = {
  "in=???\n         Input snapshot",
  "outdir=???\n     Run directory",
  "out=\n           Output snapshot if to convert back to NEMO?",
  "nr=128\n         some help",
  "nth=64\n         some help",
  "nph=128\n        some help",
  "lmax=32\n        some help",
  "model=0\n        reads particles on mout00.bin file (where n=id_new) **do not change**",
  "mond_ind=1\n     0=Newtonian, 1=MOND, 2=deepMOND",
  "spl_ord=2\n      some help",
  "rmap=1\n         some help",
  "a0=1.0\n         some help",
  "dt_iter=0.3\n    some help",
  "scale=1.0\n      some help",
  "iter_max=30\n    some help",
  "tol=5\n          some help",
  "new=0\n          some help",
  "id_new=0\n       index to read from for initial conditions",
  "nout=10\n        ** number of snapshots to output",
  "iene=5\n         some help",
  "mrate=10\n       some help",
  "tmax=10.0\n      ** end of integration",
  "dt_min=1.e-5\n   some help",
  "cfl=0.3\n        some help",
  "lp_ord=2\n       some help",
  "VERSION=0.3\n    17-sep-2013 PJT",
  NULL,
};

string usage="Frontend to the rmond program";

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
  string exefile = "rmond";
  string parfile = "input.data";
  char dname[256], cmd[256];

  nr = getiparam("nr");
  nth = getiparam("nth");
  nph = getiparam("nph");
  lmax = getiparam("lmax");
  model = getiparam("model");
  mond_ind = getiparam("mond_ind");
  spl_ord = getrparam("spl_ord");
  rmap = getrparam("rmap");
  a0 = getrparam("a0");
  dt_iter = getrparam("dt_iter");
  scale = getrparam("scale");
  iter_max = getiparam("iter_max");
  tol = getrparam("tol");
  new = getiparam("new");
  id_new = getiparam("id_new");
  nout = getiparam("nout");
  iene = getiparam("iene");
  mrate = getiparam("mrate");
  tmax = getrparam("tmax");
  dt_min = getrparam("dt_min");
  cfl = getrparam("cfl");
  lp_ord = getiparam("lp_ord");




  run_mkdir(outdir);
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
  run_sh(cmd);

  run_cd(outdir);
  sprintf(cmd,"rmond");
  run_sh(cmd);

  if (hasvalue("out")) {
    sprintf(cmd,"cat mout??.bin | modysnap - %s",outfile);
    run_sh(cmd);
  }

}
