/*
 * 	Select lines from a file
 *
 *      9-mar-99    V1.0    Created
 *     10-mar-99    V1.0a   code cleanup, no pipe check if no selection
 *     14-oct-99    V1.1    added comment= keyword
 *
 */

#include <stdinc.h>
#include <getparam.h>

string defv[] = {
	"in=???\n		input file",
        "out=???\n              output file",
        "select=all\n           lines to select",
	"comment=t\n		count comment lines too?",
        "nmax=10000\n           Default max allocation",
	"VERSION=1.1b\n		8-dec-01 PJT",
	NULL,
};

string usage="Select lines from a file";

#ifndef MAX_LINELEN 
#define MAX_LINELEN  2048
#endif

nemo_main()
{
    stream istr, ostr;
    char line[MAX_LINELEN];
    int nmax,  *select = NULL;
    int nout, next = 0;
    char *selstring=getparam("select");
    bool Qsel = !streq(selstring,"all");
    bool Qcom = getbparam("comment");
    int    i, j;
    string iname = getparam("in");

    if (Qsel) {
        nmax = nemo_file_lines(iname,getiparam("nmax"));
        if (nmax<0) error("Error opening %s",iname);
        if (nmax==0) error("No data in %s?",iname);
        dprintf(1,"NMAX=%d\n",nmax);
        select = (int *) allocate(nmax * sizeof(int));
        nout = nemoinpi(selstring,select,nmax);
        if (nout == -23) error("%s: selected too many lines (max=%d)",
                selstring,nmax);
	if (nout > nmax || select[nout-1] > nmax)
	    warning("Selected too many? input=%d output=%d max=%d",
                    nmax,nout,select[nout-1]);
    } 

    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    i = j = 0;   /* i counts lines, j points into the sorted 'select' array */
    if (Qsel) next = select[j];
    while (fgets(line,MAX_LINELEN,istr) != NULL) {
        i++;
        if (!Qcom && line[0]=='#') i--;
        if (Qsel) {
            dprintf(1,"::: %d %d %d %s\n",i,j,next, i<next ? "skip" : "sel");
            if (i < next) continue;
            fputs(line,ostr);
            if (j<nout-1)
            	next = select[++j];
            else
            	next = nmax+1;
        } else
            fputs(line,ostr);
    }
    if (!Qsel) nout = i;
    strclose(istr);
    strclose(ostr);
    dprintf(1,"Read %d lines, Written %d lines\n",i,nout);
}

