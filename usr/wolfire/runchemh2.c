/*
 *  RUNCHEMH2:   run the CHEMH2 program
 *
 *    assemble parameters
 *    needs chemie6.dat and optionally npdr.in from the Meudon code
 *    create directory, run program
 *
 *     5-feb-2018   0.1  quick hack
 */

#include <nemo.h>
#include <run.h>

string defv[] = {
  "rundir=???\n          Run directory",
  "dat=chemie6.dat\n     Example Input datafile",
  "DENS0=\n              2nd line...",
  "ZETACR=\n             ...",
  "G0=\n                 ...",
  "ABUNC=\n              ...",
  "ABUNO=\n              ...",
  "XPRES=\n              3rd line...",
  "ABUNMG=\n             ...",
  "ABUNSI=\n             ...",
  "ABUNFE=\n             ...",
  "ABUNS=\n              ...",
  "ABUNF=\n              ...",
  "ABUNCL=\n             ...",
  "FGPUMP=\n             5th line...",
  "IBRLO=\n              ...",
  "ISO=\n                ...",
  "ITURB=\n              ...",
  "ITHP=\n               ...",
  "VERSION=0.1\n         5-feb-2018 PJT",
  NULL,
};

string usage="Frontend to the CHEMH2 program";

string cvsid="$Id$";






local string fn1[] = {
  "chemie6.dat",
  "euv633.dat",
  "marks_spect.snreuv.dat4",
  NULL,
};

local string fn2[] = {
  "npdr.in",
  NULL,
};

void nemo_main()
{
  stream parstr, outstr, datstr;
  string datfile = getparam("dat");
  string outdir = getparam("rundir");
  string exefile = "chemh2";
  string parfile = "start.txt";
  string logfile = "chemh2.log";
  char dname[256], cmd[256];
  char line[256];

  /* 1st line:  no parameters ever change */
  /* 1000 format(18I3)
  /* NR,NMOL,IEL,NREACT,NPHOT1,NPHOT2,NCR1,NCR2,IH2,IH,IC,IE,IO,ICP,IH2O,IOH,ICH */

  /* 2nd line:  underlines parameters can change */
  /* 1005 format(2x,f6.3,2x,e8.2,2x,f6.3,6(1PE8.1))   */
  /* DENS0,ZETACR,G0,ABUNC,ABUNO,DVDOP,YPE,XD,EB */
  /* ----- ------ -- ----- -----                 */

  /* 3rd line: */
  /* 1006 format(2x,f6.3,10e8.1) */
  /* XPRES,ALUM,TD0,TEM,ABUNMG,ABUNSI,ABUNFE,ABUNS,ABUNF,ABUNCL */
  /* -----              ------ ------ ------ ----- ----- ------ */

  /* 4th line: no parameters ever change */
  /* T1000,T2000 */

  /* 5th line: */
  /* format(*)
  /* FGPUMP,IBRLO,ISO,ITURB,ITHP */
  /* -----  ----- --- ----- ---- */

  /* 6th line and below:  they never change */

  
  run_mkdir(outdir);

  sprintf(dname,"%s/%s",outdir,fn1[0]);      /* chemie6.dat */
  parstr = stropen(dname,"w");

  datstr = stropen(datfile,"r");

  for(;;) {
    if (get_line(datstr, line) < 0) break;
    fprintf(parstr,"%s\n",line);
  }
  strclose(parstr);


  // run_sh(cmd);     need to get symlinks to the other files

  run_cd(outdir);
  sprintf(cmd,"%s > %s ", exefile, logfile);
  run_sh(cmd);

  outstr = stropen("history.bin","w");
  put_history(outstr);
  strclose(outstr);
}

