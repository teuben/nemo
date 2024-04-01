/*
 *     TABROWS:	Select rows (lines) from a file - formerly called TABLINES
 *
 *      9-mar-99    V1.0    Created
 *     10-mar-99    V1.0a   code cleanup, no pipe check if no selection
 *     14-oct-99    V1.1    added comment= keyword
 *     10-mar-2022  V1.2    use new table interface
 *      5-may-2022  V2.0    new name (tablines -> tabrows)
 *     30-mar-2024  V2.2    select table from in-memory if no comments
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <table.h>

string defv[] = {
  "in=???\n		input file",
  "select=\n            rows to select (1=first)",
  "nmax=1000000\n       Default max allocation for lines to be selected [to be deprecated]",
  "comment=f\n          Also count comment lines as part of the rows",
  "out=-\n              output file",
  "VERSION=2.2\n	30-mar-2024 PJT",
  NULL,
};

string usage="Select non-comment rows from a file";

void nemo_main()
{
    stream istr, ostr;
    table *tptr = NULL;
    int i, j, nmax, *select = NULL;
    bool Qsel = hasvalue("select");
    string s, iname = getparam("in");
    bool Qcom = getbparam("comment");
    int next=0, nout=0, tmode = 0;

    if (Qcom) tmode = 1;
    istr = stropen(getparam("in"),"r");
    ostr = stropen(getparam("out"),"w");

    if (Qsel) {
        char *selstring=getparam("select");      
        // @todo   need to grab nmax from tptr->nr if tmode=0
	if (tmode==0) {
	  tptr = table_open(istr,0);
	  nmax = tptr->nr;
	} else 
	  nmax = getiparam("nmax");
	dprintf(1,"nmax=%d\n",nmax);
        if (nmax<0) error("Error opening %s",iname);
        if (nmax==0) error("No data in %s?",iname);
        select = (int *) allocate(nmax * sizeof(int));
        nout = nemoinpi(selstring,select,nmax);
        if (nout == -23) error("%s: selected too many lines (max=%d)",
                selstring,nmax);
	if (nout > nmax || select[nout-1] > nmax)
	    warning("Selected too many? input=%d output=%d max=%d",
                    nmax,nout,select[nout-1]);
    } 

    if (tmode==0) {   // comment lines are ignored here
      if (tptr == NULL)
	tptr = table_open(istr,0);
      if (Qsel) {
	if (nout==1 && select[0] == -1) select[0] = tptr->nr;
      } else {
	nout = tptr->nr;
	select = (int *) allocate(nout * sizeof(int));
	for (i=0; i<nout; i++)
	  select[i] = i+1;
      }
      for (i=0; i<nout; i++) {
	if (select[i] < 1 || select[i] > tptr->nr) error("Row %d illegal. minmax = 1 %d", select[i], tptr->nr);
	fprintf(ostr,"%s\n", table_row(tptr, select[i]-1));
      }
    } else {         // comment lines are not ignored, and part of the (now ordered) line selection
      warning("mode not properly re-implemented yet");

      // check
      j=0;
      for (i=1; i<nout; i++)
	if (select[i] < select[i-1]) j++;
      if (j>0) warning("select= not ordered for comment=t");

      
      // the old code
      tptr = table_open(istr,1);
      i = j = 0;   /* i counts lines, j points into the sorted 'select' array */
      if (Qsel) next = select[j];
      while ( (s=table_line(tptr)) ) {
        i++;
        if (s[0]=='#') i--;   //  don't count comment lines
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
    }

    strclose(istr);
    strclose(ostr);
    dprintf(1,"Read %d lines, Written %d lines\n",tptr->nr,nout);
}

