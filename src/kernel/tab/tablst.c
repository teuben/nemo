/*
 * 	List a (binary) table 
 *
 *	27-oct-88	Peter Teuben
 *	 1-dec-88	Added ascii option 	PJT
 *	14-nov-90 	helpvec			PJT
 *       9-jan-95       fixed strcpy() bug, made return address static
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>
#include <history.h>
#include <table.h>

string defv[] = {
	"in=???\n		  input table (filestruct) ",
        "set=\n                   set name table ",
	"col=???\n		  names of columns ",
	"skip=0\n                 skip in front ",
	"stride=1\n               stride through table ",
	"format=%f\n              format for printf ",
	"VERSION=1.3a\n		  1-jan-04 PJT",
	NULL,
};

string usage="list a (binary) table";

#ifndef MAXCOL
#define MAXCOL 64
#endif

extern string *burststring(string,string);

string get_setname(string);

void nemo_main()
{
    stream instr;
    int    i, j, ntab, ncol, skip, stride;
    real   *x[MAXCOL];      /* pointers to the columns */
    char   *set, *col, *fmt;
    string *colname;

    instr = stropen(getparam("in"),"r");
    set = getparam("set");
    if (set==NULL || *set==0)
        set = get_setname(getparam("in"));  /* ??? I/O redir ??? */
    col = getparam("col");
    colname = burststring(col,", ");
    skip = getiparam("skip");
    stride = getiparam("stride");
    fmt = getparam("format");
    
    for (i=0; i<MAXCOL; i++)
        x[i] = NULL;
    
    if (*set) {       /* valid filestruct set */
        get_history(instr);
        get_set (instr,set);
        get_data (instr,"Ntab", IntType, &ntab, 0);
        ncol = 0;                  /* set pointer (to produce columns) */
        while (colname[ncol] != NULL) {
            dprintf(0,"%s ",colname[ncol]);
            x[ncol] = (real *) malloc ( ntab * sizeof(real));
            if (x[ncol] == NULL)
                error("Cannot allocate column");
            get_data_coerced(instr, colname[ncol], RealType, x[ncol], ntab, 0);
            
            ncol++;
        }
        dprintf(0,"\n");
#if 0
        get_tes (instr,set);
        strclose(instr);
#endif
        for (i=0; i<ntab; i++) {
            for (j=0; j<ncol; j++) {
                printf (fmt,*(x[j]+i));
                printf(" ");
            }
            printf("\n");
        }
    } else {
        warning ("ASCII version not yet installed ....");
    }
}


string get_setname (string fname)
{
    permanent char buffer[256];		/* pointer to this space is returned */
    int i;
    stream instr;

    dprintf (0,"Attempt to get setname from file %s\n",fname);
	
	/* SHOULD USE POPEN HERE : */
    sprintf(buffer,"tsf %s | grep set > _get_setname_out",fname);
    system (buffer);

    instr = stropen("_get_setname_out","r");
    if (get_line (instr, buffer) < 0) {
        dprintf (0,"GET_SETNAME: no set found in %s\n",fname);
	buffer[0] = '\0';
        return (buffer);
    }
    if (strncmp(buffer,"set ",4)==0) {
	i=4;
	while (buffer[i]!='\0' && buffer[i]!=' ' && buffer[i]!='\n')
	    i++;
	buffer[i] = '\0';
        return (&buffer[4]);
    } else {
        sprintf(buffer,"Error: no set found in %s",fname);
        return(buffer);
    }
}
