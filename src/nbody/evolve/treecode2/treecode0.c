/*
 *  TREECODE0:   Nemo driver to the TREECODE family
 *
 *	Writes a script that:
 *	Converts input snapshot(5) to the "205" ascii format
 *	runs the (fortran) treecode program, and converts
 *	the output data (in 205 format) back to snapshot(5)
 *	The script is then executed via execl().
 *
 *	The following executables need to be present in the
 *	searchpath:
 *	   stoa, treecode*, and atos
 *	   
 *
 * History:
 *	26-jul-91	Created		PJT
 *	23-oct-91       extra tr pipe to delete , in transfer	PJT
 */
#include <stdinc.h>
#include <getparam.h>

string defv[] = {
                    /* NEMO parameters */
  "prog=treecode2\n             Name of TREECODE program to run",
  "in=???\n                     Snapshot input filename",
  "out=???\n                    Snapshot output filename",
  "log=???\n                    Logfile from TREECODE program",
  "clean=t\n                    Clean up temp files?",
  "exec=t\n			Execute the script?",
                    /* TREECODE parameters */
  "text=Test Run\n              Random verbiage for this run",
  "nsteps=64\n                  total # steps for this run",
  "noutbod=8\n                  # steps before major output of mas,phase,phi",
  "noutlog=1\n                  # steps before output in logfile",
  "dt=0.03125\n                 timestep",
  "tol=1.0\n                    tolerance",
  "eps=0.05\n                   softening length",
  "usequad=f\n                  Use quadrupole corrections?",
  "VERSION=1.1\n                23-oct-91 PJT",
  NULL,
};


nemo_main()
{
    stream str;
    char tempname[100], *iname, *oname, *pname, *olog;
    int i, pid;
    bool Qclean, Qexec;

    iname = getparam("in");
    str = stropen(iname,"r");           /* open - to test if present */
    strclose(str);                      /* and close */
    pname = getparam("prog");
    if (*pname == NULL) error("Need prog= name");
    oname = getparam("out");		/* use oname for temp dir */
    olog = getparam("log");

    /* make script */
    pid = getpid();
    sprintf(tempname,"tmpscript.%d",pid);
    dprintf(1,"[Writing temporary shell script %s\n",tempname);
    str = stropen(tempname,"w");
    fprintf(str,"#! /bin/csh -f\n");
    fprintf(str,"# File created by %s - do not edit\n",getargv0());
    fprintf(str,"set odir=$cwd;mkdir %s_dir;if($status)goto error;cd %s_dir\n",
                 oname,oname);
    fprintf(str,"stoa $odir/%s TREEBI\n",
                 iname);

    fprintf(str,"echo %s > TREEPAR\n",getparam("text"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("nsteps"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("noutbod"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("noutlog"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("dt"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("tol"));
    fprintf(str,"echo %s>> TREEPAR\n",getparam("eps"));
    if (getbparam("usequad"))
        fprintf(str,"echo .TRUE.>> TREEPAR\n");
    else
        fprintf(str,"echo .FALSE.>> TREEPAR\n");
    fprintf(str,"%s ; cd ..\n",
                 pname);
    fprintf(str,"tr D e < %s_dir/TREEBOA|tr -d ,|atos - %s mass,phase,phi\n",
                 oname, oname);
    fprintf(str,"cat %s_dir/TREELOG > %s\n",
                 oname, olog);
    fprintf(str,"error:\ncd $odir\n");
    if (getbparam("clean"))
        fprintf(str,"rm -rf %s_dir; rm -f %s\n",
                 oname, tempname);
    strclose(str);

    if (getbparam("exec")) {
        dprintf(0,"[Created and executing script %s]\n",tempname);
        execl("/bin/csh","csh",tempname,NULL);
        error("execl: %s\n",tempname);
    } else 
        dprintf(0,"[Created script %s]\n",tempname); 
}
