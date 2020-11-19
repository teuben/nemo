/*
 * TABMATH: general table manipulator, mini spreadsheet calculator
 *
 *      18-may-88  V1.0 created         Peter Teuben
 *       1-Jun-88  V1.1 name changed from 'nemotable' to 'tabmath'
 *                      fixed a small bug for interactive input
 *      13-jun-88  V1.1a        escape character is % now
 *      14-jun-88  V1.2 keyword 'stride' added
 *      23-aug-88  V1.3 keyword 'in2' added as extra option
 *      28-oct-88  V1.4 multiple new columns in 'newcol' now possible, 
 *			%0=linenr. 
 *      10-nov-88  V1.5 allow tab's also as column separators
 *                     a:  changed herinp to nemoinp
 *       5-jul-89      c:  changed nemoinp to nemoinpd
 *	13-nov-90  V1.6 all Xrange -> nemoinpX
 *
 *	18-feb-92  V2.0: using fie() routines to parse, not herinp!
 *                       the fie() routines are now in 'real' (not float)
 *			 this increases the program speed by about 4!
 *                       since no more double <-> float conversions needed
 *                       in= can handle multiple files, makes in2= redundant
 *	19-feb-92	 Deleted the in2= keyword....
 *	18-may-92   (b)	 added seed=0 in case randum numbers are generated
 *	 4-mar-94   (c)  ansi headers
 *	 1-sep-95  V2.1  more prototypes; also allow '!' as starting comment
 *	20-nov-96  V2.1a skip blank lines (old behavior: no more input)
 *	11-jun-98      b fix problem with TAB's in lines
 *      13-jun-98  V3.0  add selfie= to select rows based on column evaluations
 *                       but deleted the skip & stride keywords. 
 *      19-apr-01  V3.1  added comments=
 *       8-sep-01     a  init_xrandom
 *       1-aug-02     b  using nemo_debug()
 *      13-nov-03  V3.3  handle null's (now that herinp sort of knows what to do?)  [unfinished]
 *      31-dec-03  V3.4  added colname=
 *       1-jan-04     a  changed interface to get_line
 *
 *
 */

/**************** INCLUDE FILES ********************************/ 


#include <stdinc.h>
#include <getparam.h>
#include <table.h>
#include <extstring.h>
#include <ctype.h>

/**************** COMMAND LINE PARAMETERS **********************/

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n           input file name(s)",
    "out=???\n          output file name",
    "newcol=\n          formula for new column (fie notation)",
    "delcol=\n          columns to skip while writing",
    "selfie=\n          select rows using fie notation",
    "format=%g\n        format for new output columns",
    "select=all\n	Select lines (not implemented)",
    "seed=0\n           Initial random number",
    "colname=\n         (unchecked) commented column names to add into output",
    "comments=f\n       Pass through comments?",
    "refie=f\n          Re-FIE each output column (not used)",
    "VERSION=3.5b\n     19-nov-2020 PJT",
    NULL
};

string usage = "General table manipulator";



/**************** GLOBAL VARIABLES *******************************/

string *inputs, output;			/* file names */
stream *instr, outstr;			/* file streams */
int   ninput;				/* number of input files */

string fmt;                             /* format of new column */

#define MAXCOL          256             /* MAXIMUM number of columns */
#define MLINELEN       8196		/* linelength of catenated */
#define MNEWDAT          80		/* space needed for one number */

bool   keepc[MAXCOL+1];                 /* columns to keep (t/f) */
int    ndelc;                           /* actual number of skip columns */

string *fies, selfie;                   /* fie pointers */
int    nfies;                           /* number of fie pointers */
bool   Qfie;				/* boolean if multiple fie's loaded */
bool   Qrefie;                          /* recompute each columns via fie ? */
string *colname=NULL;                   /* names of columns */

bool   Qcomment;

local void setparams(void);
local void convert(int, stream *, stream);
local string *burstfie(string);
local void tab2space(char *);

extern  string *burststring(string, string);
extern  int inifie(string);
extern void dofie(real *, int *, real *, real *);
extern void dmpfie(void);
extern int savefie(int);
extern int loadfie(int);

/****************************** START OF PROGRAM **********************/

void nemo_main(void)
{
    int i;

    setparams();

    for (i=0; i<ninput; i++)
        instr[i] = stropen(inputs[i],"r");
    outstr = stropen (output,"w");

    convert (ninput,instr,outstr);

}

local void setparams(void)
{
    string newcol;                          /* formula for new column */
    string delcol;                          /* which columns not to write */
    int    delc[MAXCOL];                    /* columns to skip for output */
    int i;

    inputs = burststring(getparam("in"),", \t");
    ninput = xstrlen(inputs,sizeof(string)) - 1;
    instr = (stream *) allocate(ninput*sizeof(stream));

    output = getparam("out");

    newcol = getparam("newcol");

    delcol = getparam("delcol");
    for (i=0; i<MAXCOL; i++)
        keepc[i]=TRUE;
    if (delcol==NULL || *delcol=='\0')
        ndelc = 0;
    else if (streq(delcol,"all") || streq(delcol,"0"))
        ndelc = -1;
    else {
        ndelc = nemoinpi(delcol,delc,MAXCOL);
        if (ndelc<=0 || ndelc>MAXCOL)
            error("Too many columns given (%d) to delete (MAXCOL %d)",
			ndelc,MAXCOL);
        for (i=0; i<ndelc; i++)
            keepc[delc[i]]=FALSE;
    }
    fmt = getparam("format");
    fies = burstfie(newcol);
    nfies = xstrlen(fies,sizeof(string)) - 1;
    if(nfies)dprintf(1,"%d functions to parse\n",nfies);
    for (i=0; i<nfies; i++) {
	dprintf(1,"Saving: %s\n",fies[i]);
        inifie(fies[i]);
        if (savefie(i+1) < 0) error("Could not save fie[%d]: %s\n",i,fies[i]);
	if(nemo_debug(1)) dmpfie();
    }
    Qfie = nfies > 1;
    selfie = getparam("selfie");
    if (*selfie) {
        inifie(selfie);
	if (savefie(nfies+1) < 0) error("Could not save selfie=%s",selfie);
        Qfie = nfies > 0;
    }
    init_xrandom(getparam("seed"));
    Qcomment = getbparam("comments");
    if (hasvalue("colname"))
      colname = burststring(getparam("colname"),",");
    Qrefie = getbparam("refie");
}



local void convert(int ninput, stream *instr, stream outstr)
{
    char   line[MLINELEN];          /* input linelength */
    real   dval[MAXCOL];            /* number of items (values on line) */
    real   retval;
    char   newdat[MNEWDAT];         /* to store new column in ascii */
    int    nval, i, nlines, one=1;
    string *outv;                   /* pointer to vector of strings to write */
    char   *cp, *seps=", \t";       /* column separators  */
    real   errval=0.0;

    if (colname) {
      nval = xstrlen(colname,sizeof(string))-1;
      fprintf(outstr,"#");
      for (i=0; i<nval; i++)
	fprintf(outstr," %s",colname[i]);

      fprintf(outstr,"\n");
    }
        
    nlines=0;               /* count lines read so far */

    for(;;) {                              /* loop over all lines in file(s) */

        for(i=0, cp=line; i<ninput; i++) { /* append all lines into one line */
            if (get_line(instr[i], cp) < 0)           
                return; 					      /* EOF */
            if(iscomment(cp)) {
	      if (Qcomment)
		fprintf(outstr,"%s",cp);
	      else
		continue;	               	  /* don't use comment lines */
	    }
            if(i<ninput-1) {                         /* if not the last one: */
                cp += strlen(cp);		       /* append a blank for */
                *cp++ = ' ';	      		               /* catenation */
            }
        }
        dprintf(3,"LINE: (%s)\n",line);
        if (iscomment(line)) {
	  if (Qcomment)
	    fprintf(outstr,"%s\n",line);
	  continue;
	}
        nlines++;
        tab2space(line);	          /* work around a Gipsy (?) problem */
        if (nfies>0 || *selfie) {              	/* if a new column requested */
            nval = nemoinpr(line,dval,MAXCOL);         /* split into numbers */
	    if (nval < 0) error("bad parsing in %s",line);
	    /* this could contain some NULL's, so how do we measure this ??? */
            dprintf (3,"nval=%d \n",nval);
            if (nval>MAXCOL)
                error ("Too many numbers: %s",line);
            for(i=0; i<nfies; i++) {
                if (Qfie) loadfie(i+1);
                dofie(dval,&one,&dval[nval+i],&errval);
                dprintf(3," dofie(%d) -> %g\n",i+1,dval[nval+i]); //BUG
                strcat(line," ");
                sprintf(newdat,fmt,dval[nval+i]);
                dprintf (2,"newdat=%s\n",newdat);
                strcat(line,newdat);
            }
        }
	if (*selfie) {                             /* check if row selected */
            if (Qfie) loadfie(nfies+1);
            dofie(dval,&one,&retval,&errval);
            if (retval == 0.0) continue;            
            /* ?? but what about this stride business ?? */
	}
        if (ndelc==0) {                      /* nothing to skip while output */
            strcat (line,"\n");     
            fputs (line,outstr);
        } else {		           /* something to skip while output */
            outv = burststring(line,seps);
            i=0;
            while (outv[i]) {
                if (keepc[i+1] && (ndelc>0 || i>=nval)) {
                    fputs(outv[i],outstr);
                    fputs(" ",outstr);
                }
                i++;
            }
            fputs("\n",outstr);
        }
    } /* for(;;) */
}

/* burstfie(): to be placed with burststring() later on...
 *
 *	18-feb-92	written		PJT
 */

#define MWRD  256	/* max words in list */
#define MSTR  256	/* max chars per word */

local string *burstfie(string lst)
{
    string wrdbuf[MWRD], *wp;
    char strbuf[MSTR], *sp, *lp;
    int level=0;

    wp = wrdbuf;
    sp = strbuf;
    lp = lst;
    do {						/* scan over list */
	if (*lp == 0 || (*lp == ',' && level == 0)){	/*   is this a sep? */
	    if (sp > strbuf) {				/*     got a word? */
		*sp = 0;
		*wp++ = (string) copxstr(strbuf, sizeof(char));
		if (wp == &wrdbuf[MWRD])		/*       no room? */
		    error("burststring: too many words");
		sp = strbuf;				/*       for next 1 */
	    }
	} else {					/*   part of word */
            if (*lp == '(') level++;
            if (*lp == ')') level--;
	    *sp++ = *lp;				/*     so copy it */
	    if (sp == &strbuf[MSTR])			/*     no room? */
		error("burststring: word too long");
	}
    } while (*lp++ != 0);				/* until list ends */
    *wp = 0;
    return ((string *) copxstr(wrdbuf, sizeof(string)));	/*PPAP*/
}

/*
 * small helper function, replaces tabs by spaces before processing.
 * this prevents me from diving into gipsy parsing routines and fix
 * the problem there 
 * PJT - June 1998.
 */

local void tab2space(char *cp)
{
    while (*cp) {
        if (*cp == '\t') *cp = ' ';
        cp++;
    }
}
