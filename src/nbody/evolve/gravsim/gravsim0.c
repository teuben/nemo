/*
 *  GRAVSIM0:   Nemo driver to the GRAVSIM family
 *
 *	Writes a script that:
 *	Converts input snapshot(5) via the "205" ascii format
 *	to gdff, adds the appropriate sections to the gdff
 *	file via standard NEMO keywords, runs the gravsim
 *	program, and converts the output datafile back to
 *	205, which is turned into output snapshot(5) again....
 *	The script is then executed via execl().
 *
 *	The following executables need to be present in the
 *	searchpath:
 *	   stoa, 205-to-gdff, gdff-editor, gravsim*, gdff-to-205 and atos
 *	   
 *
 * History:
 *	28-may-91	Created		PJT
 */
#include <stdinc.h>
#include <getparam.h>

string defv[] = {
/* NEMO parameters */
  "in=???\n                         Snapshot input filename",
  "out=???\n                        Snapshot output filename",
  "prog=gravsim\n                   GRAVSIM Program to run",
  "clean=t\n			    Clean up temp files?\n\
            Remaining parameters are for GRAVSIM family",
/* GRAVSIM prameters */
  "title=\n	                    Title for the run",
  "processors=1\n                   Number of processors that can be used",
  "cells=\n                         Number of cells to be allocated",
/*  "phi-and-acceleration-output=f\n  Extra output for phi/acc?",   */
  "split-number=\n                  ...",
  "start-time=0.0\n                 Start time for run",
  "end-time=2.0\n                   Stop time for run",
  "eps=0.05\n                       Softening length",
  "tol=1\n                          Tolerance",
  "data-step=0.25\n                 Major Output of data",
  "report-step=0.03125\n            Minor Output of energy etc.",
  "g-value=1\n                      Gravitational constant",
  "time-step=0.03125\n              Integration time step",
  "VERSION=1\n                      28-may-91 PJT",
  NULL,
};

string usage = "NEMO driver for gravsim N-body integrator";

#define NEMOPAR 3       /* define how many are NEMO pars, and GRAVSIM pars */

void
nemo_main()
{
    stream str;
    char tempname[100], *iname, *oname, *pname, *cp, key[20];
    int i, pid;

    iname = getparam("in");
    str = stropen(iname,"r");           /* open - to test if present */
    strclose(str);                      /* and close */
    pname = getparam("prog");
    if (*pname == NULL) error("Need prog= name");
    oname = getparam("out");

    /* make script */
    pid = getpid();
    sprintf(tempname,"gravsim%d",pid);
    dprintf(1,"[Writing temporary shell script %s\n",tempname);
    str = stropen(tempname,"w");
    fprintf(str,"#! /bin/csh -f\n");
    fprintf(str,"# File created by %s - do not edit\n",getargv0());
    fprintf(str,"stoa %s %s.1\n",iname,tempname);
    fprintf(str,"205-to-gdff %s.1 %s.2\n",tempname,tempname);
    fprintf(str,"gdff-editor %s.2 << END_OF_EDIT >& /dev/null\n",tempname);

    for (i=NEMOPAR; cp=defv[i]; i++) { /* process all args - except first few */
        if (cp==NULL) break;
        if (*cp==NULL) break;
        if (strncmp(cp,"VERSION",7)==0) break;
        strncpy(key,cp,20);
        for(cp=key; *cp; cp++)
            if (*cp=='=') {
                *cp=NULL;
                break;
            }
        cp = getparam(key);
        if (*cp == NULL) continue;
        dprintf(1,"Adding %s=%s\n",key,cp);
        fprintf(str,"add %s %s\n",key,cp);
    }
    fprintf(str,"exit\n");
    fprintf(str,"END_OF_EDIT\n");
    fprintf(str,"%s %s.2\n",pname,tempname);
    fprintf(str,"gdff-to-205 %s.2 %s.3\n",tempname,tempname);
    fprintf(str,"atos %s.3 %s\n",tempname,oname);
    if (getbparam("clean"))
        fprintf(str,"rm -f %s %s.1 %s.2 %s.3\n",
                        tempname,tempname,tempname,tempname);
    strclose(str);

    dprintf(0,"[Created and executing script %s]\n",tempname);
    execl("/bin/csh","csh",tempname,NULL);
    error("execl: %s\n",tempname);
}
