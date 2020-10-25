/*
 *   SCANFITS:   scan fits files, optionally repairs (and outputs) them
 *
 *	27-mar-90       V1.0 created under NEMO		PJT
 *      17-jul-90       V1.1 added fix for Prof. Z (bloody Brits)    PJT
 *      18-jul               more: blank=, bloody IRAF programmers  PJT
 *       7-mar-91       V1.3 added Tuscon fixes for Marc Pound  PJT
 *	11-mar-91	     calling burststring() instead	PJT
 *	25-jul-91	V1.4 added IRAF0 option, again for Dr. Z	PJT
 *				and fixed the blocking factor bug	
 *	30-jul-91	     added DECORDER check for those Caltech files
 *				written by XYZ (no names mentioned)
 *	11-apr-92	V1.5 renamed file= to hdu=	PJT
 *	 8-jul-92	V1.6 added new fix=MSDOS 		PJT
 *      13-jan-93       V1.7 added insert= to allow new keywords    PJT
 *	18-jan-93	V1.7b multiple insert= allowed		    PJT
 *	22-feb-94           c ansi header usage
 *	 8-oct-94	    d upgraded for new extra argument in fts_cdata  PJT
 *	14-mar-95	    e added BSWAP option			PJT
 *	23-may-95           f doc
 *       2-dec-98       V1.8  allow out= filename to contains printf()
 *                            structure when split=t
 *                            added fix=promote     XTENSION->SIMPLE    PJT
 *	12-feb-99	    a changes for new fts_cdata			PJT
 *
 */

#include <stdinc.h>
#include <getparam.h>
#include <fits.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits file",
    "out=\n	           Output fits file (optional)",
    "hdu=0\n               Select which HDU's? [0=all]",
    "delete=\n   ...Delete (blank out) headers which match any of these",
    "keep=\n     ...Or keep headers which match any of these",
    "insert=\n             file(s) with header items to insert before END",
    "fix=\n                What kind of fix [ING,TUCSON,LF,CRLF,BLANK]\n\
	ING      Shift left 13 chars special ING card images\n\
	TUCSON   Shift left 18 chars special SINGLDSH card images\n\
	LF       Make ascii editable by inserting a LF in column 80\n\
	BLANK    Delete blank lines (used with keep= and delete=)\n\
	ZERO     Replace zeros by blanks\n\
	DECORDER Delete DECORDER, if present, and then swap the bytes\n\
	BSWAP    swap the bytes anyways\n\
	IRAF0	 Replace HISTORY by COMMENT\n\
	MSDOS	 Replace any non-printable header ASCII by blanks\n\
	FIRST    Retain only first occurance of true keywords\n\
	UNY2K    Un-fix Y2K from yyyy-mm-dd to dd/mm/yy\n\
	PROMOTE  Promote XTENSION to SIMPLE header\n\
	LAST     Retain only last occurance of true keywords",
    "print=\n              Additional items to print from header?",
    "blocking=1,1\n  	   Two blocking factors (blocksize/2880) for i&o",
    "select=header,data\n  Select header, data, ...?",
    "split=f\n             Split fits file into HDU pieces '<out>.#'",
    "VERSION=1.8c\n        25-oct-2020 PJT",
    NULL,
};

string usage="scan fits files, optionally repair (and output) them";


extern bool scanopt(string, string);
extern string *burststring(string, string);

nemo_main()
{
    stream instr, outstr;
    int    i, n, nfile, blocking[2];
    string outfile, hdselect, *insert, *fix, *deletes, *keep, *print;
    char   basename[128];
    struct fits_header fh;
    bool   split, sel_head, sel_data;

    instr = stropen(getparam("in"),"r");    /* open input */
    split = getbparam("split");     /* to split or not to split */

    fix = burststring(getparam("fix"),",");            /* List of  fixes */
    deletes = burststring(getparam("delete"),",");
    keep = burststring(getparam("keep"),",");
    if (*keep != NULL && *deletes != NULL)
        error("Can only handle one of keep= and delete=");
    print =  burststring(getparam("print"),",");
    insert = burststring(getparam("insert"),",");

    
    if (hasvalue("out")) {                       /* optional output */
        outfile = getparam("out");          
        if (!split) 				/* if no split */
            outstr = stropen(outfile,"w");	/* open file now */
	else
	    outstr = NULL;			/* else set NULL */
    } else {
        outstr = NULL;			/* if no output at all */
        split = FALSE;			/* flag these variables */
    }
    nfile = getiparam("hdu");		/* need: nemoinpi() */

    if (hasvalue("select")) {            /* header/data selection */
        hdselect = getparam("select"); 
        sel_head = scanopt(hdselect,"header");
        sel_data = scanopt(hdselect,"data");
    } else {
        sel_data = sel_head = TRUE;
    }
    n = nemoinpi(getparam("blocking"),blocking,2);
    if (n<1 || (outstr!=NULL && n<2))
        error("Not enough blocking factors supplied");

    fts_setiblk(blocking[0]);    /* set input blocking factor */
    fts_setoblk(blocking[1]);    /* set input blocking factor */
	
    for (i=1;;i++) {			             /* loop over all HDU's */
       fts_zero(&fh);			             /* reset header */
       fh.hdu = i;                                   /* keep track of HDU (1=first) */
       
       n = fts_rhead(&fh,instr);	              /* read header */
       if (n<0)				              /* if no data (EOF) .. */
          break;			              /* ... quit */
       if (outstr) dprintf(0,"Working on FITS file %d\n",i);
       fts_dhead(&fh,deletes);                     /* delete= headers */
       fts_khead(&fh,keep);                        /* keep= headers */
       fts_ihead(&fh,insert);                      /* insert= headers */
       fts_fhead(&fh,fix);                         /* fix= headers */
       if (!outstr && !split) fts_phead(&fh,print);  /* ? print header */
       if ((outstr && (nfile==0 || nfile==i)) || split) {
          if (split) {
          	if (strchr(outfile,'%'))
          	    sprintf(basename,outfile,i);
          	else
                    sprintf(basename,"%s.%d",outfile,i);
                dprintf(0,"Working on split FITS file %s\n",basename);
                outstr = stropen(basename,"w");
          }
          if (sel_head) fts_chead(&fh,outstr);        /* copy header */
          if (sel_data) 
              fts_cdata(&fh,instr,outstr,TRUE,TRUE);  /* copy data, with trailing bits */
          else
              fts_sdata(&fh,instr);	              /* skip the data */
          if (split) strclose(outstr);
       } else                                          /* Otherwise just */
          fts_sdata(&fh,instr);		               /* skip the data */
       if (i==nfile)
            break;                                        /* done */
    }
    if (nfile>0 && i>=nfile && n<0)
        stop(1);                    /* see if early eof */

    strclose(instr);                    /* close files */
    if (outstr) strclose(outstr);
}
