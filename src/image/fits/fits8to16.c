/*
 *   FITS8TO16:   convert 8bit fits to 16bit fits
 *
 *	1-jul-04      quick and dirty              PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <fits.h>

string defv[] = {			/* Standard NEMO keyword+help */
    "in=???\n              Input fits file",
    "out=\n	           Output fits file (optional)",
    "VERSION=1.0\n         1-jul-04 PJT",
    NULL,
};

string usage="convert 8bit fits image to a 16bit fits image";


nemo_main()
{
  stream instr, outstr;
  int    i, n, nfile, blocking[2];
  string outfile, hdselect, *insert, *fix, *deletes, *keep, *print;
  char   basename[128];
  struct fits_header fh;
  bool   split, sel_head, sel_data;
  
  instr = stropen(getparam("in"),"r");    /* open input */

  outfile = getparam("out");          
  outstr = stropen(outfile,"w");	/* open file now */

  for (i=1;;i++) {			             /* loop over all HDU's */
    fts_zero(&fh);			             /* reset header */
    n = fts_rhead(&fh,instr);	              /* read header */
    if (n<0)				              /* if no data (EOF) .. */
      break;			              /* ... quit */
    fts_chead816(&fh,outstr);                         /* copy header */
    fts_cdata816(&fh,instr,outstr);         /* copy data, with trailing bits */
                                               /* this also modifies fh->bitpix */
  }
  strclose(outstr);
  strclose(instr);                    /* close files */
}
