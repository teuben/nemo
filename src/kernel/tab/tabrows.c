/*
 *     TABROWS:	Select rows (lines) from a file - formerly called TABLINES
 *
 *      9-mar-99    V1.0    Created
 *     10-mar-99    V1.0a   code cleanup, no pipe check if no selection
 *     14-oct-99    V1.1    added comment= keyword
 *     10-mar-2022  V1.2    use new table interface
 *      5-may-2022  V2.0    new name (tablines -> tabrows)
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <table.h>

string defv[] = {
	"in=???\n		input file",
        "select=all\n           rows/lines to select",
	"comment=t\n		count comment lines too?",
        "nmax=10000\n           Default max allocation for lines to be picked",
        "out=-\n                output file",
	"VERSION=2.0\n		5-may-2022 PJT",
	NULL,
};

string usage="Select rows/lines from a file";

void nemo_main()
{
    stream istr, ostr;
    table *tptr;
    string s;
    int nmax,  *select = NULL;
    int nout, next = 0;
    char *selstring=getparam("select");
    bool Qsel = !streq(selstring,"all");
    bool Qcom = getbparam("comment");
    int    i, j;
    string iname = getparam("in");

    if (Qsel) {
        // @todo   relic from old table interface, this needs a more dynamic interface
        nmax = nemo_file_lines(iname,getiparam("nmax"));
	dprintf(0,"nmax=%d\n",nmax);
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

    tptr = table_open(istr,1);  // this is a streaming app, no need for pipe support

    i = j = 0;   /* i counts lines, j points into the sorted 'select' array */
    if (Qsel) next = select[j];
    while ( (s=table_line(tptr)) ) {
        i++;
        if (!Qcom && s[0]=='#') i--;
        if (Qsel) {
            dprintf(1,"::: %d %d %d %s\n",i,j,next, i<next ? "skip" : "sel");
            if (i < next) continue;
            fprintf(ostr,"%s\n",s);
            if (j<nout-1)
            	next = select[++j];
            else
            	next = nmax+1;
        } else
            fprintf(ostr,"%s\n",s);	  
    }
    if (!Qsel) nout = i;
    strclose(istr);
    strclose(ostr);
    dprintf(1,"Read %d lines, Written %d lines\n",i,nout);
}

