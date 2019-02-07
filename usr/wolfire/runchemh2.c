/*
 *  RUNCHEMH2:   run the CHEMH2 program
 *
 *    assemble parameters
 *    patch chemie6.dat and optionally npdr.in from the Meudon code
 *    create directory, run program
 *
 *     5-feb-2018   0.1  quick initial hacks
 */

#include <nemo.h>
#include <filefn.h>
#include <extstring.h>
#include <run.h>

string defv[] = {
  "rundir=???\n          Run directory",
  "exe=chemh2\n          Name of executable (in $PATH)",
  "dat=chemie6.dat\n     Example Input datafile",
  "DENS0=\n              2nd line...",
  "ZETACR=\n             ...",
  "G0=\n                 ...",
  "ABUNC=\n              ...C abundance",
  "ABUNO=\n              ...O abundance",
  "DVDOP=\n              ...line width (km/s??)",
  "XPRES=\n              3rd line...",
  "ABUNMG=\n             ...Mg abundance",
  "ABUNSI=\n             ...Si abundance",
  "ABUNFE=\n             ...Fe abundance",
  "ABUNS=\n              ...S abundance",
  "ABUNF=\n              ...F abundance",
  "ABUNCL=\n             ...Cl abundance",
  "FGPUMP=\n             5th line...",
  "IBRLO=\n              ...meudon flag",
  "ISO=\n                ...",
  "ITURB=\n              ...",
  "ITHP=\n               ...",
  "GRIDSTEP=25,29,0,6\n  Written to the GRIDOUTPUT/GRIDSTEP file",
  "VERSION=0.6\n         7-feb-2018 PJT",
  NULL,
};

string usage="Frontend to the CHEMH2 program";

string cvsid="$Id$";




/* TODO's
 *
 *     - the dat= file should be able to be read , edited, and written. now they both need to co-exist.
 */


/* files needed in runtime,
 * perhaps these should be listed in an acii file given by tab= or so?
 *
 */

local string fn1[] = {          // unit
  "chemie6.dat",                // 10
  // "chemie7.dat",                // 10   fatal error? or soon the new format?
  "euv633.dat",                 // 20 
  "marks_spect.snreuv.dat4",    // 21
  "avrout.dat",                 // 22
  "contable.dat",               // 29
  "twophasepg10fh20.dat",       // 59
  "twophasepg10fh21em4.dat",    // 60
  "twophasepg10fh21em3.dat",    // 61
  "twophasepg10fh23em3.dat",    // 62
  "twophasepg10fh21em2.dat",    // 63
  "twophasepg10fh23em2.dat",    // 64
  "twophasepg10fh21em1.dat",    // 65
  "twophasepg10fh23em1.dat",    // 66
  "twophasepg10fh25em1.dat",    // 67
  "twophasepg10fh27em1.dat",    // 68
  "twophasepg10fh29em1.dat",    // 69
  "q_oh2_12c16o.dat",           // 23
  "q_ph2_12c16o.dat",           // 23
  // "GRIDOUTPUT",                 // this is a whole complex tree of read and write
  NULL,
};

local string fn2[] = {
  "npdr.in",
  NULL,
};

local string gridoutput_files =
  "CI CO CP FeI FeII H2 HCOP OI SI SiI SiII TS TSav0p01 TSav0p03 TSav0p1 TSav0p3 TSav1p0 TSav3p0";


local string patch_f(string key,string line,int start,int len, string fmt);
local string patch_e(string key,string line,int start,int len, string fmt);
local string patch_i(string key,string line,int start,int len, string fmt);

#define MAXGRIDSTEP  4



void nemo_main()
{
  stream parstr, outstr, datstr, tmpstr;
  string datfile = getparam("dat");
  string outdir  = getparam("rundir");
  string exefile = "chemh2";
  string parfile = "start.txt";
  string logfile = "chemh2.log";
  string l5[5];
  char dname[256], cmd[256];
  char line[256];
  int i, nf1, nf2, ibrlo;
  int gridstep[MAXGRIDSTEP];

  /* fetch/check some variables */
  if (nemoinpi(getparam("GRIDSTEP"),gridstep,MAXGRIDSTEP) != MAXGRIDSTEP)
    error("GRIDSTEP= needs 4 integers");

  /* (debug) report data files needed */
  
  nf1 = xstrlen(fn1, sizeof(string))-1;
  nf2 = xstrlen(fn2, sizeof(string))-1;
  for (i=0; i<nf1; i++)
    dprintf(1,"CHEMH2 data1: %s\n",fn1[i]);
  for (i=0; i<nf2; i++)
    dprintf(1,"CHEMH2 data2: %s\n",fn2[i]);

  /* create run directory */

  run_mkdir(outdir);

  /*************************************************************************************/  
  /* read input par file, copy to output par in run directory while editing parameters */
  
  datstr = stropen(datfile,"r");             /* input par, normally chemie6.dat, but can be anything */
  sprintf(dname,"%s/%s",outdir,fn1[0]);      /* output par name should chemie6.dat */
  parstr = stropen(dname,"w");

		
  /* 1st line:  no parameters ever change */
  /* 1000 format(18I3) */
  /* 200 76 12322 16 33 92 95  5  1  3 10  4 13  8  9  7 25  */
  /* NR,NMOL,IEL,NREACT,NPHOT1,NPHOT2,NCR1,NCR2,IH2,IH,IC,IE,IO,ICP,IH2O,IOH,ICH */
  
  get_line(datstr, line);
  fprintf(parstr,"%s\n",line);

  /* 2nd line:  underlined parameters can change */
  /* 1005 format(2x,f6.3,2x,e8.2,2x,f6.3,6(1PE8.1))   */
  /*    1.000  0.20E-15   0.750 1.6E-04 3.2E-04 1.5E+00 1.0E-01 4.4E-01 2.7E+00 */
  /* DENS0,ZETACR,G0,ABUNC,ABUNO,DVDOP,YPE,XD,EB */
  /* ----- ------ -- ----- -----                 */

  get_line(datstr, line);
  patch_f("DENS0",  line, 2,6, "%6.3f");
  patch_e("ZETACR", line,10,8, "%8.2e");
  patch_f("G0",     line,20,6, "%6.3f");
  patch_e("ABUNC",  line,26,8, "%8.1e");
  patch_e("ABUNO",  line,34,8, "%8.1e");
  patch_e("DVDOP",  line,42,8, "%8.1e");  
  fprintf(parstr,"%s\n",line);

  /* 3rd line: */
  /* 1006 format(2x,f6.3,10e8.1) */
  /*    0.000 0.0E+00 0.0E+00 0.0E+00 1.1E-06 1.7E-06 1.7E-07 2.8E-05 1.8E-08 1.8E-07 */
  /* XPRES,ALUM,TD0,TEM,ABUNMG,ABUNSI,ABUNFE,ABUNS,ABUNF,ABUNCL */
  /* -----              ------ ------ ------ ----- ----- ------ */

  get_line(datstr, line);
  patch_f("XPRES",  line, 2,6, "%6.3f");    /* TODO: these needs a classic fortran E8.1 */
  patch_f("ABUNMG", line,32,8, "%8.1e");  
  patch_f("ABUNSI", line,40,8, "%8.1e");  
  patch_f("ABUNFE", line,48,8, "%8.1e");  
  patch_f("ABUNS",  line,56,8, "%8.1e");  
  patch_f("ABUNF" , line,64,8, "%8.1e");  
  patch_f("ABUNCL", line,72,8, "%8.1e");  
  fprintf(parstr,"%s\n",line);

  /* 4th line: no parameters ever change */
  /*  2.5E+00 1.8E+00 */
  /* 1001 format(9E8.1)
  /* T1000,T2000 */

  get_line(datstr, line);
  fprintf(parstr,"%s\n",line);

  /* 5th line: */
  /*  1.0 1 0 0 0 */
  /* format(*)
  /* FGPUMP,IBRLO,ISO,ITURB,ITHP */
  /* -----  ----- --- ----- ---- */

  get_line(datstr, line);
  l5[0] = patch_f("FGPUMP",  line,0,0, "%g");    
  l5[1] = patch_i("IBRLO",   line,1,0, "%d");  
  l5[2] = patch_i("ISO",     line,2,0, "%d");  
  l5[3] = patch_i("ITURB",   line,3,0, "%d");  
  l5[4] = patch_i("ITHP",    line,4,0, "%d");
  sprintf(line," %s %s %s %s %s",l5[0],l5[1],l5[2],l5[3],l5[4]);
  fprintf(parstr,"%s\n",line);

  /* 6th line and below:  they never change */

  for(;;) {
    if (get_line(datstr, line) < 0) break;
    fprintf(parstr,"%s\n",line);
  }
  strclose(parstr);

  ibrlo = atoi(l5[1]);
  if (ibrlo)
    warning("CHEMH2 running with the Merged Meudon code; not checked yet");

  /*************************************************************************/

  run_cd(outdir);

  //   get symlinks to the datafiles
  string fn, chempath = getenv("CHEMPATH");
  if (chempath == NULL) chempath = strdup("..");
  dprintf(1,"CHEMPATH: %s\n",chempath);
    
  for (i=1; i<nf1; i++) {  // skip the first one, chemie6.dat, we wrote it manually
    fn = pathfind(chempath,fn1[i]);
    if (fn != NULL) {
      // dprintf(1,"SYMLINK %s %s\n",fn1[i],fn);
      // symlink(fn,fn1[i]);
      sprintf(cmd,"cp -a %s .", fn);
      dprintf(1,"%s\n",cmd);      
      run_sh(cmd);
    } else
      warning("Missing file %s along $CHEMPATH\n",fn1[i]);
  }
  for (i=0; i<nf2; i++) {  // skip the first one, chemie6.dat, we wrote it manually
    fn = pathfind(chempath,fn2[i]);
    if (fn != NULL) {
      // dprintf(1,"SYMLINK %s %s\n",fn2[i],fn);
      // symlink(fn,fn2[i]);
      sprintf(cmd,"cp -a %s .", fn);
      dprintf(1,"%s\n",cmd);
      run_sh(cmd);
    } else
      warning("Missing file %s along $CHEMPATH\n",fn2[i]);
  }

  /* deal with GRIDOUTPUT */
  mkdir("GRIDOUTPUT",0775);

  tmpstr = stropen("GRIDOUTPUT/GRIDSTEP","w");
  fprintf(tmpstr,"%12d%12d%12d%12d\n",gridstep[0],gridstep[1],gridstep[2],gridstep[3]);

  sprintf(cmd,"cd GRIDOUTPUT; touch %s", gridoutput_files);
  run_sh(cmd);

  sprintf(cmd,"%s > %s ", exefile, logfile);
  run_sh(cmd);

  outstr = stropen("history.bin","w");
  put_history(outstr);
  strclose(outstr);
}


/*
 * The following 3 patch routines either patch a line with parameters in place (len>0)
 * or return the string representing the N'th (start) parameter on the line if len==0
 * The latter is needed if free format(*) was used.
 * Since line is overwritten in case len>0 this routine cannot be mixed between these
 * two types of calls .  After all, this is a hack.
 * To make matters worse, there may also be some memory leaks.
 *
 */

string patch_f(string key,string line,int start,int len, string fmt)
{
  int i;
  real val;
  char fval[64];
  string *sp, rv;

  if (len==0) {
    sp = burststring(line," ");
    if (hasvalue(key))
      rv = getparam(key);
    else {
      rv = strdup(sp[start]);
      freestrings(sp);
    }
    dprintf(1,"KEY: %s[%d] = %s\n",key,start,rv);    
    return rv;
  }
  
  if (!hasvalue(key)) return;
  val = getrparam(key);
  sprintf(fval,fmt,val);
  dprintf(1,"KEY: %s=%f  %s->%s\n",key,val,fmt,fval);
  for (i=0; i<len; i++)
    line[start+i] = fval[i];
  return fval;
}

string patch_e(string key,string line,int start,int len, string fmt)
{
  int i;
  real val;
  char fval[64];

  // no len=0 needed here, since we never use that option
  
  if (!hasvalue(key)) return;
  val = getrparam(key);
  sprintf(fval,fmt,val);
  dprintf(1,"KEY: %s=%f  %s->%s\n",key,val,fmt,fval);
  for (i=0; i<len; i++)
    line[start+i] = fval[i];
  return fval;
}

string patch_i(string key,string line,int start,int len, string fmt)
{
  int i;
  int val;
  char fval[64];
  string *sp,rv;
  
  if (len==0) {
    sp = burststring(line," ");
    dprintf(1,"KEY: %s len=%d\n",key,xstrlen(sp,sizeof(string)));
    if (hasvalue(key))
      rv = getparam(key);
    else {
      rv = strdup(sp[start]);
      freestrings(sp);
    }
    dprintf(1,"KEY: %s[%d] = %s\n",key,start,rv);
    /* memleak */
    return rv;
  }
  
  if (!hasvalue(key)) return;
  val = getiparam(key);
  sprintf(fval,fmt,val);
  dprintf(1,"KEY: %s=%d  %s->%s\n",key,val,fmt,fval);
  for (i=0; i<len; i++)
    line[start+i] = fval[i];
  return fval;
}
