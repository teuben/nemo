/*
 *   FITSHEAD:   display a fits header, raw
 *		 or convert a piece of text to fits header format...
 *
 *	9-apr-91        V1.0 created from scanfits		PJT
 *	13-apr-92	V1.1 renamed file= to hdu=		PJT
 *	22-feb-94       V1.1a ansi headers			pjt
 *	 5-apr-94	V1.2  use out= to force header conversion	PJT
 *      23-may-95       V1.2b fixed bug in convert_to_header		pjt
 *	 7-aug-01       V1.3  added counter=t|f
 *                            -- oops, needs a library change --
 */

#include <stdinc.h>
#include <getparam.h>
#include <fits.h>

string defv[] = {	/* Standard NEMO keyword+help */
    "in=???\n              Input fits file",
    "hdu=0\n               Which HDU (0=all, 1=1st etc.)",
    "blocking=1\n          Blocking factor (blocking/2880)",
    "out=\n                Convert input text to output fits header",
    "counter=f\n           Add line counter to output?",
    "VERSION=1.3a\n        10-aug-09 PJT",
    NULL,
};

string usage = "display, or convert to, a fits header";


extern string *burststring(string, string);

nemo_main()
{
    if (hasvalue("out"))
        convert_to_header();
    else
        read_fits_header();
}

read_fits_header()
{
    stream instr, outstr;
    int    i,n,nfile, sel_data, sel_head, blocking, counter;
    size_t dsize;
    string outfile, select, *fix, *delete, *keep, *print;
    char   basename[128];
    struct fits_header fh;
    bool   split, Qcount = getbparam("counter");


    instr = stropen(getparam("in"),"r");    /* open input */

    nfile = getiparam("hdu");		/* need: nemoinpi() */

    fts_setiblk(getiparam("blocking"));
	
    for (i=1;;i++) {			             /* try infinite loop */
       fts_zero(&fh);			             /* clean out header */
       dsize = fts_rhead(&fh,instr);	             /* read header */
       dprintf(1,"%d: dsize=0x%x (%ld)\n",i,dsize,dsize);
       if (dsize < 0)			             /* if no data (EOF) .. */
          break;			             /* ... quit */

       if (nfile==0 || nfile==i)
          fts_thead(&fh);
       fts_sdata(&fh,instr);	                     /* alawys skip the data */
       if (i>=nfile)
            break;                                   /* all done */
    }

    strclose(instr);                    /* close files */
}


#define MAXLEN 256

convert_to_header()
{
    stream instr, outstr;
    char line[MAXLEN+1];
    int i, len, n = 0, nline=0;


    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    while ( fgets(line,MAXLEN,instr) ) {
        n++;
        len = strlen(line);
        if (len > MAXLEN)
            error("line %d way too long (%d); possibly binary file",
                    n,len,MAXLEN);

        if (len > FTSLINSIZ) {
            warning("line %d too long (%d): will be chopped", n,len);
            line[FTSLINSIZ] = 0;
        } else {
            /* correct for newline in last spot */
            for (i=len-1; i<FTSLINSIZ; i++) line[i] = ' ';
            line[FTSLINSIZ] = 0;
        }
        len = fwrite(line,1,FTSLINSIZ,outstr);
        if (len != FTSLINSIZ) error("Error writing fitsheader at line %d",n);
        nline++;
    }

    if (strncmp(line,"END",3) != 0) {
        nline++;
        warning("Appending standard END line to fitsheader at line %d",nline);
	strcpy(line,"END");
	for (i=3; i<FTSLINSIZ; i++) line[i] = ' ';
        len = fwrite(line,1,FTSLINSIZ,outstr);
        if (len != FTSLINSIZ) error("Error writing fitsheader at line %d",nline);
    }

    for (i=0; i<FTSLINSIZ; i++) line[i] = ' ';
    if (nline % FTSLPB) {
        n = FTSLPB - nline % FTSLPB;
        while (n-- > 0) {
            nline++;
            len = fwrite(line,1,FTSLINSIZ,outstr);
            if (len != FTSLINSIZ) error("Error writing fitsheader at line %d",nline);
        }
    }

    strclose(instr);
    strclose(outstr);
}
