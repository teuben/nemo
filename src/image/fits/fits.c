/*
 * FITS.C: various fits-i/o routines
 *
 *      Created  March-90       Peter Teuben - 
 *      
 *      Updates:
 *              27-mar-90       gcount is sometimes floating ???        PJT
 *              24-jun-90       few more keywords, more err checking    PJT
 *              20-jul-90       added fts_cdata(), double, lot more chk PJT
 *		21-feb-91       check_ctype                             PJT
 *		11-feb-91	fts_dhead, fts_khead, fts_fhead bug	PJT
 *		16-mar-91	various fixes, additions etc. for panel PJT
 *		27-jun-91	added more ctype's for checking		PJT
 *		25-jul-91	added IRAF0 fix to _fhead		PJT
 *                              did something about blocking factors    PJT
 *		30-jul-91	DECORDER fix				PJT
 *		 4-aug-91	'fixed' fts_chead; only 1 fwrite()	PJT
 *		   mar-92	   few gcc2 happy's
 *		 8-apr-92	added some TABLE support (reading)	PJT
 *               7-aug-92       fixed reading ctype's and ttype's       PJT
 *				added cunit's as recognized primary array
 *                              wrote check_unit, added tfields output to
 *                              fts_phead
 *              23-sep-92       fixed some (non)(random) GROUPS confusion
 *                              and other minor bugs
 *		13-jan-93	added fts_ihead to insert header cards  PJT
 *              14-jan-93       added fix_multi for redundant options   PJT
 *		18-jan          changed calling params of fts_ihead     PJT
 *				attemp to handle BINTABLE in ptable	PJT
 *		29-mar          ptable: fixed bug with row= in multiple HDUs
 *		 9-aug-93       allow ascii headers at input .....      PJT
 *              17-aug-93       fixed various bintable display bugs     PJT
 *		24-jan-94       patched header output for scanfits      pjt
 *		22-feb-94       ansi (header)
 *		29-jun-94       added fh->nread datamember and support
 *		 1-oct-94 	warn about deprecated *MJD* keywords?	pjt
 *		 7-oct-94	added argument to fts_cdata		pjt
 *		 9-oct-94	error if bad (non-fits) input format	pjt
 *              14-mar-95       added BSWAP to force 'DECORDER'         pjt
 *		23-may-95       fixed bug in fts_thead			pjt
 *               4-mar-96       fixed protype conflict log10()
 *              15-oct-96       fixed some bugs in ptable/bintable      pjt
 *               1-apr-97       added proto
 *              11-feb-99       added xhead, for rawfits.c		pjt
 *              15-oct-99       fixed Y2K problem			pjt
 *              28-nov-00       fixed some NULL -> 0's                  pjt
 *              30-may-01       proper prototypes                       pjt
 *              20-jun-01       a few gcc3 fixes
 *              28-sep-01       'K' format for 64 bit integers (*long*) pjt
 *                              endian swap is done in the show_ routines
 *                              which is bad !!!!
 *    		12-oct-01       new standard added FITS reference
 *              15-jan-03       write floating #'s in header w/ more precision PJT
 *              16-jan-03       warning if one of CROTA1/2 is missing      PJT
 *               3-may-04       add fts_read_img_coord                     PJT
 *               4-may-04       conform more to fits standard              PJT
 *               1-jul-04       quick fix transform BITPIX8 to 16 images   PJT
 *              30-dec-04       don't shorted blank keyval's to appease VALGRIND and common sense ?     PJT
 *                              also made a3..a5 1 char longer to make room for the terminating 0
 *              10-nov-05       newline option in fts_pgroup,fts_ptable 
 *              11-dec-06       store cvsID in saved string
 *              10-aug-09       int -> size_t in a few more places for big files
 *
 * Places where this package will call error(), and hence EXIT program:
 *  - invalid BITPIX
 *  - incomplete I/O (e.g. file too short, full filesystem)
 *  - not enough memory to allocate needed space
 *  - 
 *
 */

#include <stdinc.h>             /* NEMO's stdio.h */
#include <ctype.h>              /* needs: isdigit() */
#include <fits.h>

#if SIZEOF_LONG_LONG==8
typedef long long int8;         /* e.g. i386; sparc <= sol7; ppc ? */
#elif SIZEOF_LONG==8
typedef long int8;              /* e.g. alpha 64's */
#elif SIZEOF_INT==8
typedef int int8;               /* will never happen ? */
#else
#error FITS.C code probably needs to be checked for this architecture
typedef long int8;              /* some stupid fallback, probably wrong */
#endif


local char   *fts_buffer = NULL;        /* buffer for exact header in memory */
local size_t  fts_buflen = 0;

local int ftsblksiz_i=FTSBLKSIZ;	/* blocksize for i and o */
local int ftsblksiz_o=FTSBLKSIZ;
local int ftslpb_i=FTSLPB;		/* lines per block for i and o */
local int ftslpb_o=FTSLPB;

local string cfits1="FITS (Flexible Image Transport System) format is defined in 'Astronomy";
local string cfits2="and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H";
local string cfits3="nemo::fits.c $Id$";


/* local functions in fits.c */
static int fix_ing(fits_header *fh);
static int fix_tucson(fits_header *fh);
static int fix_nl(fits_header *fh, int mode);
static int fix_blank(fits_header *fh);
static int fix_zero(fits_header *fh);
static int fix_msdos(fits_header *fh);
static int fix_uny2k(fits_header *fh);
static int fix_y2k(fits_header *fh);
static int fix_iraf0(fits_header *fh);
static int fix_multi(fits_header *fh, int first);
static int fix_decorder(fits_header *fh);
static int fix_bswap(fits_header *fh);
static int fix_promote(fits_header *fh);
static void set_tbcoln(fits_header *fh);
static int put_ivals(char *key, char *val, char *fitsname, int *nvals, int **ivals);
static int put_dvals(char *key, char *val, char *fitsname, int *nvals, double **dvals);
static int put_cvals(char *key, char *val, char *fitsname, int *nvals, string **cvals);
static int parse_card(int icard, char *card, char *a1, char *a2, char *a3, char *a4, char *a5);
static int atob(char *cp);
static char *atoa(char *s);
static void my_copy(char *src, char *dest, int n);
static void poutline(FILE *fp);
static void check_date(char *key, char *s);
static int check_ctype(char *key, char *val);
static int check_unit(char *key, char *val);
static int my_findstr(string text, string pat, int len);
static void blank_fill(char *card, int len);
static int colitems(char *s);
static char *fmt(char *line, char *format);
static float show_e(float *v, int i);
static double show_d(double *v, int i);
static short show_i(short *v, int i);
static int   show_j(int *v, int i);
static int8  show_k(int8 *v, int i);
static int colmask(int n, string key[],string col[],int colsel[]);


/*
 *  fts_rhead:  reads a fits header, does a lot of checking too
 *              returns:    -1    nothing read
 *                          -k    (k>1) incomplete header
 *                           k    (k>=0) size of data (in bytes) ahead
 *  if you need to fake a header, of fake you've read a header,
 *  use fts_xhead
 */

size_t fts_rhead(fits_header *fh, stream instr)
{
    char buf[FTSLINSIZ+2];                /* buffer to hold one card image */
    char a1[9], a2[2], a3[FTSLINSIZ+1],a4[FTSLINSIZ+1],a5[FTSLINSIZ+1];  /* args */
    int k, n, i, icard;
    struct arglist *p;

    if (fh->naxis >= 0) {         /* we could force a fts_zero here, but ala */
        dprintf(0,"### fts_rhead: fits_header not initialized (use debug>0)\n");
        dprintf(1,"    it seems to still have a naxis=%d\n",fh->naxis);
        dprintf(1,"    unpredictable things might happen\n");
    }
    if (sizeof(int8) != 8)
      warning("Cannot handle BITPIX=64 of tables K format, sizeof(long)=%d",sizeof(long));

    for(icard=1;;icard++) {  /* loop reading cards; count them in 'icard' */

        /* reading a card should always succeed, since the END card  */
        /* breaks this infinite for(icard) loop - if not, return < 0 */
        
        if (fh->ascii) {           /* horrendous misFITS trick */
            if (fgets(buf,FTSLINSIZ,instr) == 0) return -1; /* !! never !! */
            n = strlen(buf);
            if (buf[n-1] == '\n') n--;
            for (i=n; i<FTSLINSIZ; i++) buf[i] = ' ';	/* patch it */
            n = FTSLINSIZ;
        } else {               /* this is of course the legal way to get cards */
          if ((n=fread(buf,1,FTSLINSIZ,instr)) != FTSLINSIZ) {     /* get card */
            if (n==0)
                return -1;         /* error - nothing read -- never -- */
            else 
                return -n-1;       /* return some negative number, but < -1 */
          } else {
            buf[FTSLINSIZ] = 0;      /* be sure it's terminated */
            dprintf(5,"RECORD[%d]%65s\n",icard,buf);
	  }
	}

        fh->hlen += n;   /* keep track of amount of proper header in fts_buf */
        fts_buf((icard+1)*FTSLINSIZ);          /* make sure enough in buffer */
        my_copy( buf, &fts_buffer[(icard-1)*FTSLINSIZ], FTSLINSIZ); /* append */

        if (buf[0] == ' ') {      /* just skip blank lines - old FITS flaw */
            dprintf(2,"### (Card %d) Skipping blank line\n",icard);
            continue;
        }
	buf[FTSLINSIZ] = 0;                         /* make sure it ends */
        n = parse_card(icard,buf, a1,a2,a3,a4,a5);     /* break in parts */
        if (n<=0)               /* for bad looking records */
            continue;           /* go to the next one */
        
        for (p = fitsargs; p->name != NULL; p++) {        /* find keyword */
            if (strncmp(a1,p->name,strlen(p->name))==0) { /* min match */
                if (p->isarray)     /* array keyword, min match is sufficient */
                    break;  
                if (streq(a1,p->name)) /* otherwise need an exact match to match */
                    break;
            }
        }
        if (icard==1) 
            if (p->token != I_SIMPLE && p->token != I_XTENSION)
                error("First FITS keyword not SIMPLE or XTENSION");
        if (icard==2)
            if (p->token != I_BITPIX)
                error("Second FITS keyword not BITPIX");        

        if (p->name == NULL) {
            dprintf(2,"### (Card %d) skipped unknown keyword %s (%s)\n",
                           icard,a1,a3);
            continue;
        }
        if (p->token == I_END)
            break;              /* quit the infinite for-loop */

        dprintf(2,"FTS_HEADER> keyword %s %s %s\n",a1,a2,a3);

        switch (p->token) {	/* parse required and reserved keywords */
        case I_SIMPLE:
            if (icard!=1)			/* 6.2.1.1. */
                dprintf(0,"### SIMPLE keyword in card %d\n",icard);
            if (fh->simple >= 0)
                dprintf(0,"### SIMPLE redefined\n");
            fh->simple = atob(a3);
            if (fh->simple != 1)		/* 6.2.1.1.1. */
                dprintf(0,"Warning: SIMPLE = %d not supported\n",fh->simple);
            break;
        case I_XTENSION:
            if (icard!=1)			/* 5.5 */
                dprintf(0,"### XTENSION keyword in card %d\n",icard);
            fh->xtension = atoa(a3);
            dprintf(1,"Standard extension %s:\n",a3);
                            /* check for standard names */
                            /* check if not in primary header: 6.2.1.3.1 */
            if (fh->simple >=0)
                dprintf(0,"### XTENSION present in primary header\n");
            break;
        case I_BITPIX:
            if (icard!=2)   /* 6.2.1.1 and 6.2.1.3 */
                dprintf(0,"### BITPIX in card %d, should be in 2\n",icard);
            if (fh->bitpix > 0)
                dprintf(0,"### BITPIX redefined\n");
            fh->bitpix = atoi(a3);
            switch (fh->bitpix) {
                case 8:
                case 16:
                case 32:
                case 64:
                case -32:
                case -64:
                    break;      /* any of these values are OK */
                default:        /* otherwise some error:   6.2.1.1.2 */
                    error("Unsupported BITPIX = %d\n",fh->bitpix);
            }                    
            break;
        case I_NAXIS:
            put_ivals(a1,a3,p->name,&fh->naxis, &fh->naxisn);
            if (fh->naxis > MAXNAXIS || fh->naxis < 0)       /* 6.2.1.1.3 */
                dprintf(0,"### Illegal value for NAXIS: %d\n",fh->naxis);
            break;
        case I_MAXIS:
            put_ivals(a1,a3,p->name,&fh->maxis, &fh->maxisn);
            if (fh->maxis > MAXNAXIS || fh->maxis < 0)
                dprintf(0,"### Illegal value for MAXIS: %d\n",fh->maxis);
            break;
        case I_EXTEND:              /* should occur in primary header */
            if (fh->extend >= 0)    /* .. immediately after NAXISn */
                dprintf(0,"Warning: EXTEND redefined\n");
            if (fh->simple < 0)
                dprintf(0,"### EXTEND not in primary header\n");
            fh->extend = atob(a3);
            break;
        case I_BLOCKED:             /* should occur in first 36 cards */
            if (icard>36) dprintf(1,
               "Warning: Keyword BLOCKED appears in card %d > 36\n",icard);
            warning("Keyword BLOCKED deprecated");
            break;
        case I_BLANK:
            if (fh->bitpix < 0)			/* 6.2.2.2.4 */
                dprintf(0,"### Illegal BLANK keyword (%s), possibly ignored\n",
                                           a3);
            else
                fh->blank = atoi(a3);
            break;
        case I_CRPIX:               /* a float array */
            put_dvals(a1,a3,p->name,&fh->naxis, &fh->crpixn);
            break;
        case I_CRVAL:               /* a float array */
            put_dvals(a1,a3,p->name,&fh->naxis, &fh->crvaln);
            break;
        case I_CDELT:               /* a float array */
            put_dvals(a1,a3,p->name,&fh->naxis, &fh->cdeltn);
            break;
        case I_CTYPE:               /* a character array */
            put_cvals(a1,a3,p->name,&fh->naxis, &fh->ctypen);
            check_ctype(a1,a3);
            break;
        case I_CUNIT:               /* a character array */
            put_cvals(a1,a3,p->name,&fh->naxis, &fh->cunitn);
            check_unit(a1,a3);
            break;
        case I_CROTA:               /* a float array */
            put_dvals(a1,a3,p->name,&fh->naxis, &fh->crotan);
            break;
        case I_BSCALE:
            fh->bscale = atof(a3);
            break;
        case I_BZERO:
            fh->bzero = atof(a3);
            break;
        case I_BUNIT:
            fh->bunit = atoa(a3);
            check_unit(a1,a3);
            break;
        case I_RESTFREQ:
            fh->restfreq = atof(a3);
            break;
        case I_DATAMIN:
            fh->datamin = atof(a3);
            break;
        case I_DATAMAX:
            fh->datamax = atof(a3);
            break;
        case I_DATE:
            fh->date = atoa(a3);
            check_date(a1,a3);
            break;
        case I_EPOCH:
            fh->epoch = atoa(a3);
            break;
        case I_DATE_OBS:
            fh->date_obs = atoa(a3);
            check_date(a1,a3);
            break;
        case I_EXTNAME:
            fh->extname = atoa(a3);
            break;
        case I_EXTVER:
            fh->extver = atoi(a3);
            break;
        case I_EXTLEVEL:
            fh->extlevel = atoi(a3);
            break;
        case I_GROUPS:
            if (fh->groups >= 0)
                dprintf(0,"Warning: GROUPS redefined\n");
            fh->groups = atob(a3);
            break;
        case I_PCOUNT:
            fh->pcount = atoi(a3);
            if (fh->pcount < 0)
                dprintf(0,"### Illegal value for PCOUNT: %d\n",fh->pcount);
            break;
        case I_GCOUNT:
            fh->gcount = atoi(a3);
            if (fh->gcount < 0)
                dprintf(0,"### Illegal value for GCOUNT: %d\n",fh->gcount);
            break;
        case I_PZERO:
            put_dvals(a1,a3,p->name,&fh->pcount, &fh->pzeron);
            break;
        case I_PSCAL:
            put_dvals(a1,a3,p->name,&fh->pcount, &fh->pscaln);
            break;
        case I_PTYPE:
            put_cvals(a1,a3,p->name,&fh->pcount, &fh->ptypen);
            break;
        case I_PUNIT:               /* a character array */
            put_cvals(a1,a3,p->name,&fh->pcount, &fh->punitn);
            check_unit(a1,a3);
            break;
        case I_TFIELDS:	/* used in TABLE, BINTABLE */
            fh->tfields = atoi(a3);
            break;
        case I_TFORM:	/* required for TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->tformn);
            break;
        case I_TBCOL:	/* required for TABLE, internally used for BINTABLE */
            put_ivals(a1,a3,p->name,&fh->tfields, &fh->tbcoln);
            break;
        case I_TTYPE:	/* optional reserved for TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->ttypen);
            put_cvals(a1,a5,p->name,&fh->tfields, &fh->ttypen_c);
            break;
        case I_TUNIT:	/* optional reserved for TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->tunitn);
            check_unit(a1,a3);
            break;
        case I_TNULL:	/* optional reserved for TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->tnulln);
            break;
        case I_TZERO:	/* optional reserved for TABLE */
            put_dvals(a1,a3,p->name,&fh->tfields, &fh->tzeron);
            break;
        case I_TDIM:    /* optional mulitdimensional field designation/ TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->tdimn);
            break;
        case I_TDISP:   /* optional suggested display fmt f/ TABLE */
            put_cvals(a1,a3,p->name,&fh->tfields, &fh->tdispn);
            break;
        case I_TSCAL:	/* optional reserved for TABLE */
            put_dvals(a1,a3,p->name,&fh->tfields, &fh->tscaln);
            break;
        case I_INSTRUME:
            fh->instrume = atoa(a3);
            break;
        case I_TELESCOP:
            fh->telescop = atoa(a3);
            break;
        case I_ORIGIN:
            fh->origin = atoa(a3);
            break;
        case I_OBJECT:
            fh->object = atoa(a3);
            break;
        case I_OBSERVER:
            fh->observer = atoa(a3);
            break;
        case I_AUTHOR:
            fh->author = atoa(a3);
            break;
        case I_REFERENC:
            fh->referenc = atoa(a3);
            break;
        default:
            dprintf(2,"(IC %d) Keyword %s unprocessed or unknown\n",icard,a1);
        }
        if (!p->needsvalue && a2[0]=='=')
            dprintf(0,"### Keyword %s does not need an '='\n",a1);
        if (p->needsvalue && a2[0]!='=')
            dprintf(0,"### Keyword %s does not '='\n",a1);
                
    }  /* infinite for (icard) - loop */

    /* assume that END was the last keyword read */    
   
    dprintf(1,"End of header in keyword %d\n",icard);

    if (fh->ascii)
        for (i=0; i<FTSLINSIZ; i++) buf[i] = ' ';
    k = icard;
    while (k % ftslpb_i != 0) {       /* read trailing junk of header */
    	if (fh->ascii == 0)
            if(fread(buf,1,FTSLINSIZ,instr) != FTSLINSIZ)
                error("Error finishing off tail of header");
        k++;
        fts_buf((k+1)*FTSLINSIZ);          /* make sure enough in buffer */
        my_copy( buf, &fts_buffer[(k-1)*FTSLINSIZ], FTSLINSIZ);/* append */
    }
    
    /* process some additional checks */
    fts_lhead(fh);            /* process a few more sanity checks */
    return fts_dsize(fh);     /* return the size of upcoming data */
}

/*
 *  fts_xhead: fake reading a header (see rawfits.c)
 *
 */

size_t fts_xhead(fits_header *fh, stream instr,
	      int bitpix, int naxis, int *naxisn, int skip)
{
    int  i;
    char buf[80];

    if (fh->naxis >= 0) {         /* we could force a fts_zero here, but ala */
        dprintf(0,"### fts_rhead: fits_header not initialized (use debug>0)\n");
        dprintf(1,"    it seems to still have a naxis=%d\n",fh->naxis);
        dprintf(1,"    unpredictable things might happen\n");
    }

    fh->hlen = skip;
    for (i=0; i<skip; i++)
    	if (fread(buf,1,1,instr) != 1) error("Skipping data at byte %d",i+1);

    fh->simple = atob("T");
    fh->bitpix = bitpix;
    fh->naxis = naxis;
    fh->naxisn = naxisn;

    return fts_dsize(fh);     /* return the size of upcoming data */
}



char *fts_shead(fits_header *fh, string keyword)
{
  char key[9];
  int i, j, ncards;
  char *cp;
  static char line[FTSLINSIZ+1];

  sprintf(key,"%-8s",keyword);

  ncards = fh->hlen / FTSLINSIZ;      /* number of cards to process */
  for (i=0; i<ncards; i++) {          /* scan through all cards ... */
    cp = &fts_buffer[i*FTSLINSIZ];    /* point to card */
    if (strncmp(cp,key,8)==0) {       /* if key matches */
      cp = &line[8];
      if (*cp != '=') continue;
      cp++;
      while (*cp == ' ') cp++;
      if (*cp == '\'' || *cp == '"') {
	error("Cannot copy strings yet");
      } else {
	i=0;
	while (*cp != ' ')
	  line[i++] = *cp++;
	line[i] = 0;
	dprintf(0,"fts_shead: %s\n",line);
	return line;
      }
    }
  }
  return NULL;
}

/*
 *   fts_lhead:   (lint) more sanity checks on a fits header
 *                This is not really required, but may prevent
 *                trouble during reading
 */

int fts_lhead(fits_header *fh)
{
    size_t dsize;

    if (fh->simple >= 0) {      /* if its a primary header */
        dsize = fts_dsize(fh);
        if (dsize<0) {
	  warning("Negative size for data? dsize=%ld\n",dsize);
	  return dsize;
        }
    } else if (fh->xtension) {
        if (streq(fh->xtension,"A3DTABLE") ||
            streq(fh->xtension,"BINTABLE") || 
            streq(fh->xtension,"TABLE")) {

            if (fh->bitpix != 8)
                warning("Illegal value of BITPIX (%d) in %s\n",
                                fh->bitpix, fh->xtension);
            if (fh->naxis != 2)
                warning("Illegal value of NAXIS (%d) in %s\n",
                                fh->naxis, fh->xtension);
            if (fh->pcount != 0)
                warning("Illegal value of PCOUNT(%d) in %s\n",
                                fh->pcount, fh->xtension);
            if (fh->gcount != 1)
                warning("Illegal value of GCOUNT(%d) in %s\n",
                                fh->gcount, fh->xtension);
        } else
            dprintf(1,"### fts_lhead: No sanity checks for extension %s\n",
                    fh->xtension);
        
    } else {                                          /* section 5.1 ??? */
        dprintf(0,"### fts_lhead: Not a known/processed FITS element\n");
    }
    return 1;
}

/*
 *  fts_chead:  copy the fits header (the one in memory)
 *
 */

int fts_chead(fits_header *fh, stream outstr)
{
    int n, nread;

    if (ftsblksiz_o > ftsblksiz_i) warning("possibly trailing garbage");

    n = fh->hlen;           /* size of header through END keyword */
    n = ROUNDUP(n,ftsblksiz_o);
    nread = fwrite(fts_buffer,sizeof(char),n,outstr);
    if (nread != n) 
        error("Tried to write %d header, could only write %d",n,nread);
#if 0
    nblank = (ROUNDUP(n,ftsblksiz_o)-n)/FTSLINSIZ;    /* # blank lines */
    while (nblank-- > 0)
        if (fwrite(blank,1,FTSLINSIZ, outstr)!=FTSLINSIZ)
            error("Writing trailing blank header lines, %d to go",nblank);
#endif
    return 1;
}

int fts_chead816(fits_header *fh, stream outstr)
{
    int n, nread;
    char *cp;

    if (ftsblksiz_o > ftsblksiz_i) warning("possibly trailing garbage");

    n = fh->hlen;           /* size of header through END keyword */
    n = ROUNDUP(n,ftsblksiz_o);

    /* patch BITPIX from 8 to 16 */
    cp = strchr(&fts_buffer[80],'8');
    if (cp == NULL) 
      error("BITPIX not 8 ???");
    *(cp-1) = '1';
    *cp = '6';

    nread = fwrite(fts_buffer,sizeof(char),n,outstr);
    if (nread != n) 
        error("Tried to write %d header, could only write %d",n,nread);
    return 1;
}

/*
 *  fts_thead:  type out the fits header (the one in memory)
 *		and do not print trailing blanks. Really meant
 *		for stdout only!!!
 *
 */

int fts_thead(fits_header *fh)
{
    int i, j, ncards;
    char *cp, line[FTSLINSIZ+1];

    ncards = fh->hlen / FTSLINSIZ;      /* number of cards to process */
    dprintf(1,"Printing %d cards\n",ncards);
    for (i=0; i<ncards; i++) {          /* scan through all cards ... */
        cp = &fts_buffer[i*FTSLINSIZ];  /* point to card */
        my_copy(cp,line,FTSLINSIZ);     /* copy a line */
	line[FTSLINSIZ] = 0;
        for (j=FTSLINSIZ-1; j>=0; j--) {   /* find last non-blank */
            if (line[j] != ' ') {
                line[j+1] = 0;
                break;
            } else
            	line[j] = 0;
        }
        if (line[0]) printf("%s\n",line);   /* only print non-blank lines */
    }
    return 1;
}

/* 
 *  fts_dhead:   delete (blank) elements from a fits header (the one in memory)
 *
 */

int fts_dhead(fits_header *fh, string *delete)
{
    char *card, blank[FTSLINSIZ];
    string *sp;
    int ncards, i, ndel;

    if (delete == NULL || *delete == NULL)
        return 0;

    for (i=0; i<FTSLINSIZ; i++)           /* prepare a blank line for patching */
        blank[i] = ' ';
    ncards = fh->hlen / FTSLINSIZ;          /* number of cards to process */
    ndel = 0;                 /* and keep a counter of the number of deletions */

    for(sp=delete; *sp; sp++) { /* process all deletions in order */
        for (i=0; i<ncards; i++) {          /* scan through all cards ... */
            card = &fts_buffer[i*FTSLINSIZ];
            if (my_findstr(card,*sp,FTSLINSIZ) >= 0){  /* got a match */
                dprintf(2,"## Deleting card %d (%-8s)\n",i+1,card);
                ndel++;
                my_copy(blank,card,FTSLINSIZ);     /* overwrite with blanks */
            }
        }
    }
    dprintf(0,"Blanked out a total of %d cards\n",ndel);
    return ndel;
}

/* 
 *  fts_khead:   keep elements from a fits header (the one in memory)
 *		 i.e. blank out the ones not mentioned
 *	Warning: This routine would also blank END, NAXIS, BITPIX etc.
 *		 if not part of the 'keep' array!
 *
 */

int fts_khead(fits_header *fh, string *keep)
{
    char *card, blank[FTSLINSIZ];
    string *sp;
    int ncards, i, ndel;
    bool keepit;

    if (keep == NULL || *keep == NULL)
        return 0;

    for (i=0; i<FTSLINSIZ; i++)           /* prepare a blank line for patching */
        blank[i] = ' ';
    ncards = fh->hlen / FTSLINSIZ;          /* number of cards to process */
    if (ncards<=0) return(0);
    ndel = 0;                 /* keep a counter of the number of deletions */
    for (i=0; i<ncards; i++) {
        card = &fts_buffer[i*FTSLINSIZ];        /* point to card */
        keepit = FALSE;              /* flags the card as to blank */
        for (sp=keep; *sp; sp++) {      /* check if match with keep's */
            if (my_findstr(card,*sp,FTSLINSIZ) >= 0) {
                keepit=TRUE;
                break;                  /* and quit the inner for() loop */
            }         
        }
        if (!keepit) {
            dprintf(2,"## Deleting card %d (%-8s)\n",i+1,card);
            ndel++;
            my_copy(blank,card,FTSLINSIZ);     /* overwrite with blanks */
        }
    }
    return 1;
}

/* 
 *  fts_ihead:   add header cards to memory header
 *
 */

int fts_ihead(fits_header *fh, string *header)
{
    char line[FTSLINSIZ+1], card[FTSLINSIZ];
    int i, n=0;    /* count how many were added */
    stream instr;

    if (header == NULL || *header == NULL)
        return 0;

    fh->hlen -= FTSLINSIZ;      /* point to at the location of END */
    for (i=0; header[i]; i++) {
        instr = stropen(header[i],"r");
        while (fgets(line,FTSLINSIZ,instr)) {
            line[strlen(line)-1] = 0;
            if (line[0] == '#' || line[0] == ';' || line[0] == '\n') continue;
    	    dprintf(2,"FTS_IHEAD: inserting card: %s\n",line);
            strncpy(card,line,FTSLINSIZ);

            blank_fill(card,FTSLINSIZ);
            fts_buf(fh->hlen + FTSLINSIZ);         /* ensure space */
            my_copy(card,&fts_buffer[fh->hlen],FTSLINSIZ);
            fh->hlen += FTSLINSIZ;
  
            n++;
        }
        strclose(instr);
    }

    if (n) {                    /* patch back the END */
        strcpy(card,"END");

        blank_fill(card,FTSLINSIZ);
        fts_buf(fh->hlen + FTSLINSIZ);         /* ensure space */
        my_copy(card,&fts_buffer[fh->hlen],FTSLINSIZ);
        fh->hlen += FTSLINSIZ;
    }

    return n;
}
/* 
 *  fts_fhead:   fix a fits header (the one in memory)
 *
 *    Options (any of the following, separated by comma's):
 *
 *      ING      line which starts with blanks and ING: shift stuff left
 *      TUCSON   line with SINGLDSH shift 18 left
 *      LF       convert to Unix style ascii files
 *      BLANK    delete blank lines, i.e. shift up 
 *	ZERO     replace 0x00 by 0x20
 *      TRIM     put a zero to shorten each line for fts_xhead() [not impl.]
 *	DECORDER Fix Jeff Hester's brilliant DECORDER keyword
 *      BSWAP    same as DECORDER, but forced swap, irrespective
 *	IRAF0	 replace 'HISTORY' by 'COMMENT'
 *      MSDOS    replace all non-printable ascii by blanks (anti CR,LF fix)
 *	PROMOTE  promote XTENSION to SIMPLE=T
 *	UNY2K	 fix down Y2K dates
 *	Y2K	 fix up old dates into Y2K format
 */

int fts_fhead(fits_header *fh, string *fix)
{
    string *sp;

    if (fix == NULL || *fix == NULL)
        return 0;

    for(sp=fix;*sp; sp++) {                 /* process all options in order */
        if (strncmp(*sp,"ING",3)==0)
            fix_ing(fh);                        /* shift keywords */
        else if (strncmp(*sp,"TUCSON",6)==0)
            fix_tucson(fh);                     /* shift keywords */
        else if (strncmp(*sp,"LF",2)==0)
            fix_nl(fh,1);                       /* make header ascii editable */
        else if (strncmp(*sp,"CRLF",4)==0)	/* undocumented/deprecated */
            fix_nl(fh,2);                       /* make header ascii editable */
        else if (strncmp(*sp,"BLANK",5)==0)
            fix_blank(fh);                      /* delete blank lines */
        else if(strncmp(*sp,"ZERO",4)==0)
            fix_zero(fh);                       /* get rid of ZERO's */
        else if(strncmp(*sp,"DECORDER",8)==0)
            fix_decorder(fh);                   /* flip the bytes, if needed */
        else if(strncmp(*sp,"BSWAP",5)==0)
            fix_bswap(fh);                      /* flip the bytes anyways */
        else if(strncmp(*sp,"IRAF0",5)==0)
            fix_iraf0(fh);                      /* HISTORY becomes COMMENT */
        else if(strncmp(*sp,"MSDOS",5)==0)
            fix_msdos(fh);                      /* All non-printable ascii become SPACE */
	else if(strncmp(*sp,"FIRST",5)==0)
	    fix_multi(fh,1);			/* keep only first occurances */
	else if(strncmp(*sp,"LAST",4)==0)
	    fix_multi(fh,0);			/* keep only last occurances */
	else if(strncmp(*sp,"UNY2K",5)==0)
	    fix_uny2k(fh);			/* turn y2k into old style */
	else if(strncmp(*sp,"Y2K",5)==0)
	    fix_y2k(fh);			/* turn old style into y2k */
	else if(strncmp(*sp,"PROMOTE",8)==0)
	    fix_promote(fh);			/* xtension -> simple   */
        else
            warning("tsf_fhead: ... skipping unknown fix %s",*sp);
    }
    return 1;
}

/*  
 * Various examples of fixes are given here.
 */

local int fix_ing(fits_header *fh)
{
    int ncards, i, j, nfix;
    char *cp;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0, nfix=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        if (strncmp(cp,"         ING ",13) != 0)
            continue;
        for (j=13; j<FTSLINSIZ; j++)
            cp[j-13] = cp[j];       /* shift left 13 */
        for (j=0; j<13; j++)
            cp[FTSLINSIZ-j-1] = ' ';    /* patch blank */
        nfix++;
    }
    dprintf(2,"A total of %d card images were ING fixed\n",nfix);
    return 1;
}

local int fix_tucson(fits_header *fh)
{
    int ncards, i, j, nfix;
    char *cp;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0, nfix=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        if (strncmp(cp,"SINGLDSH",8) != 0)
            continue;
        for (j=18; j<FTSLINSIZ; j++)
            cp[j-18] = cp[j];       /* shift left 18 */
        for (j=0; j<18; j++)
            cp[FTSLINSIZ-j-1] = ' ';    /* patch blank */
        nfix++;
    }
    dprintf(2,"A total of %d card images were ING fixed\n",nfix);
    return 1;
}


/*   mode: must be 1 (lf) or 2 (cr/lf) */

local int fix_nl(fits_header *fh, int mode)
{
    int ncards, i;
    char *cp;

    ncards = fh->hlen/FTSLINSIZ;   /* number of cards up and until END */
    ncards = ROUNDUP(ncards,ftslpb_i);

    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];  /* start of card */
        if (mode==2) cp[FTSLINSIZ-2] = 0x0D;    /* for msdos */
        cp[FTSLINSIZ-1] = 0x0A;    /* patch with newline to make it editable */
    }
    return 1;
}

local int fix_blank(fits_header *fh)
{
    int i, j, icard, ncards, b;
    char *cpf, *cpt;

    ncards = fh->hlen / FTSLINSIZ;

    for (i=0, icard=0; i<ncards; i++) {         /* loop over all legal cards */
        cpf = &fts_buffer[i*FTSLINSIZ];         /* line to check */
        for (j=0, b=1; j<FTSLINSIZ; j++) {   /* check if line has blanks */
            if (cpf[j] != ' ') {
                b = 0;  
                break;                       /* found a non-blank, so break */
            }
        }
        if (!b) {                                /* line was not blank */
            if (icard < i) {                    /* see if need to copy */
                cpt = &fts_buffer[icard*FTSLINSIZ];     /* where to go */
                my_copy(cpf,cpt,FTSLINSIZ);         /* copy it */
            }
            icard++;                            /* increase output counter */
        }
    }  /* for(all cards) */
    fh->hlen = icard * FTSLINSIZ;           /* reset header size */

    for (i=icard; i<ncards; i++) {          /* reset all other cards */
        cpt = &fts_buffer[i*FTSLINSIZ];     
        for (j=0; j<FTSLINSIZ; j++)
            cpt[j] = ' ';                   /* to fullly blank */
    }

    dprintf(1,"Out of %d cards %d were copied\n",ncards,icard);
    return fh->hlen;
}

local int fix_zero(fits_header *fh)
{
    int ncards, i, j;
    char *cp;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        for (j=0; j<FTSLINSIZ; j++)
	    if (cp[j]=='\0') cp[j] = ' ';
    }
    return 1;
}

local int fix_msdos(fits_header *fh)
{
    int ncards, i, j;
    char *cp;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        for (j=0; j<FTSLINSIZ; j++)
            if (cp[j] < 0x20 || cp[j] > 0x7E) cp[j] = ' ';
    }
    return 1;
}

local int fix_uny2k(fits_header *fh)
{
    int ncards, i, j;
    char *cp, olddate[32], newdate[32];

    /* DATE    = 'yyyy-mm-dd'  ::: doesn't handle full version yet !!! */
    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        if (strncmp(cp,"DATE",4) == 0) {
            if (cp[10]=='\'' && cp[15]=='-' && cp[18]=='-' && cp[21]=='\'')
                strncpy(olddate,&cp[10],12);
                olddate[12]= 0;
                sprintf(newdate,"'%c%c/%c%c/%c%c'  ",
                    cp[19],cp[20],cp[16],cp[17],cp[13],cp[14]);
                dprintf(0,"WARNING: Y2K fix attempt: %s -> %s\n",
                        olddate,newdate);
                for (j=0; j<strlen(newdate); j++)
                    cp[j+10] = newdate[j];
                
        }
    }
    return 1;
}

local int fix_y2k(fits_header *fh)
{
    int ncards, i, j;
    char *cp, olddate[32], newdate[32];

    /* DATE    = 'yyyy-mm-dd'  ::: doesn't handle full version yet !!! */
    /*           'dd/mm/yy'    */
    
    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
        if (strncmp(cp,"DATE",4) == 0) {
            if (cp[10]=='\'' && cp[13]=='/' && cp[16]=='/' && cp[19]=='\'')
                strncpy(olddate,&cp[10],8);
                olddate[8]= 0;
                sprintf(newdate,"'19%c%c-%c%c-%c%c'",
                    cp[17],cp[18],cp[14],cp[15],cp[11],cp[12]);
                dprintf(0,"WARNING: Y2K fix attempt: %s -> %s\n",
                        olddate,newdate);
                for (j=0; j<strlen(newdate); j++)
                    cp[j+10] = newdate[j];
                
        }
    }
    return 1;
}

local int fix_iraf0(fits_header *fh)
{
    int ncards, i, icard=0;
    char *cp;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
	if (strncmp(cp,"HISTORY ",8)==0) {
	    strcpy(cp,"COMMENT");
            cp[7] = ' ';
            icard++;
        }
    }
    dprintf(1,"IRAF0: Out of %d cards %d were fixed\n",ncards,icard);
    return 1;
}

		/* (i) keep only first (1) or last (0)  */
local int fix_multi(fits_header *fh, int first)
{
    int ncards, i, j, nblank=0;
    char *cardi, *cardj, key[9];

    dprintf(0,"fix_multi %d\n",first);

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cardi = &fts_buffer[i*FTSLINSIZ];
	if (cardi[8] != '=' ||
             strncmp(cardi,"HISTORY ",8)==0 || 
             strncmp(cardi,"COMMENT ",8)==0)
            continue;
	strncpy(key,cardi,8);  key[8]=0;
        for (j=i+1; j<ncards; j++) {
            cardj = &fts_buffer[j*FTSLINSIZ];
            if (cardj[8] != '=') continue;
            if (strncmp(key,cardj,8)==0) {      /* match ! */
                if (first) {
                    dprintf(2,"Deleting non-FIRST card %d: %s\n",j+1,key);
                    cardj[0] = ' ';     /* blank fill */
                    nblank++;
                } else {
                    dprintf(2,"Deleting non-LAST card %d: %s\n",i+1,key);
                    cardi[0] = ' ';
                    nblank++;
                    break;
                }
            }
        }
    }
    dprintf(1,"%s: Out of %d cards %d were fixed\n",
            first ? "FIRST" : "LAST",
            ncards,nblank);    
    return 1;
}

local int fix_decorder(fits_header *fh)
{
    int ncards, i;
    char *cp, *s;

    ncards = fh->hlen / FTSLINSIZ;
    for (i=0; i<ncards; i++) {
        cp = &fts_buffer[i*FTSLINSIZ];
	if (strncmp(cp,"DECORDER",8)==0) {
            s = &cp[8];
            while (*s == ' ' || *s == '=')
                s++;
            if (*s != 'T' || *s != 'F') {
                dprintf(0,"DECORDER: unknown value %c\n",*s);
                break;
            } 
            if (*s == 'F') break;

	    strcpy(cp, "COMMENT");
            cp[7] = ' ';
            fh->flip = 1;
            dprintf(0,"DECORDER: header fixed - data will be byte swapped\n");
            break;
        }
    }
    return 1;
}

local int fix_bswap(fits_header *fh)
{
    fh->flip = 1;
    return 1;
}

local int fix_promote(fits_header *fh)
{
    int ncards, nfix=0;
    char *cp, *cpi;
    static char *promote="SIMPLE  =                    T / promoted";

    ncards = fh->hlen / FTSLINSIZ;
    if (ncards > 0) {
    	cp = &fts_buffer[0];				/* look at 1st card */
        if (strncmp(cp,"XTENSION=",9) == 0) {
             cpi = promote;
             while (*cpi) {
                *cp++ = *cpi++;
             }
             
        }
    }
    dprintf(2,"A total of %d card images were PROMOTEd\n",nfix);
    return 1;

}


/*
 *   fts_phead:   print out header in some human readable format
 */

int fts_phead(fits_header *fh, string *print)
{
    int i,istart,naxis, pcount, gcount, n, ncards;
    size_t dsize;
    char *card, buf[FTSLINSIZ+1];
    char a1[9], a2[2], a3[FTSLINSIZ+1],a4[FTSLINSIZ+1],a5[FTSLINSIZ+1];  /* args */
    string *sp;

    printf("__________________________________________________________\n");
    printf("______________________ FITS HEADER _______________________\n");
    if (fh->xtension) {         /* primary header ? */
      printf("Standard eXTENSION header: %s\n",fh->xtension);
    } else {
      printf("Primary header: ");
      switch(fh->simple) {
        case 1:
            printf("standard SIMPLE\n");
            break;
        case 0:
            printf("Non-standard SIMPLE \n");
            return 0;
        default:
            printf("Illegal value for SIMPLE (%d)",fh->simple);
      }
    }
    if (fh->extend>0)     /* will only be in primary header */
        printf("Possible extension records present\n");
    if (fh->extname)    /* will only be in extension header */
        printf("Name of standard extension = %s\n",fh->extname);
    if (fh->groups == 1)
        printf("Random group format (deprecated)\n");
    naxis = fh->naxis;
    printf("BITPIX = %d\n",fh->bitpix);
    printf("NAXIS  = %d\n",naxis);
    if (naxis > 0) {
        printf("###");
        if (fh->naxisn) printf("\tNAXIS#");                 /* header */
        if (fh->ctypen) printf("\tCTYPE#  ");
        if (fh->cdeltn) printf("\tCDELT#");
        if (fh->crvaln) printf("\tCRVAL#");
        if (fh->crpixn) printf("\tCRPIX#");
        if (fh->crotan) printf("\tCROTA#");
	if (fh->cunitn) printf("\tCUNIT#");
        printf("\n");
        
        istart = ((fh->groups == 1) ? 1 : 0);
	/* BUG BUG ??? perhaps we need to enforce istart=0 ??? GROUPS=F ?? */
        for (i=istart; i<naxis; i++) {                /* array elements */
            printf("%2d:",i+1);
            if (fh->naxisn) printf("\t%d",fh->naxisn[i]);
            if (fh->ctypen) printf("\t%-8s",fh->ctypen[i]);
            if (fh->cdeltn) printf("\t%g",fh->cdeltn[i]);
            if (fh->crvaln) printf("\t%g",fh->crvaln[i]);
            if (fh->crpixn) printf("\t%g",fh->crpixn[i]);
            if (fh->crotan) printf("\t%g",fh->crotan[i]);
            if (fh->cunitn) printf("\t%-8s",fh->cunitn[i]);
            printf("\n");
        }
        if (fh->groups == 1) {             /* groupish format ? */
            pcount = fh->pcount;
            gcount = fh->gcount;
            printf("PCOUNT = %d\n",pcount);
            printf("GCOUNT = %d\n",gcount);
            printf("###");
            if (fh->ptypen) printf("\tPTYPE#");                 /* header */
            if (fh->pzeron) printf("\tPZERO#");
            if (fh->pscaln) printf("\tPSCAL#");
            printf("\n");
            for (i=0; i<pcount; i++) {                  /* array elements */
                printf("%2d:",i+1);
                if (fh->ptypen) printf("\t%-8s",fh->ptypen[i]);
                if (fh->pzeron) printf("\t%20.11g",fh->pzeron[i]);
                if (fh->pscaln) printf("\t%20.11g",fh->pscaln[i]);
                printf("\n");
            }
        }
        if (fh->tfields > 0) {                        /* extension table ? */
            printf("TFIELDS = %d\n",fh->tfields);
            printf("###");
            if (fh->tbcoln) printf("\tTBCOL#");     /* required: */
            if (fh->tformn) printf("\tTFORMN#");
            if (fh->ttypen) printf("\tTTYPEN#");    /* optional: */
            if (fh->tunitn) printf("\tTUNITN#");
            if (fh->tnulln) printf("\tTNULLN#");
            if (fh->tzeron) printf("\tTZERON#");
            if (fh->tscaln) printf("\tTSCALN#");
            printf("\n");
            /*
             *  Todo: the optional ones could have missing members, need
             *        to check for NULL
             */
            for (i=0; i<fh->tfields; i++) {                  /* array elements */
                printf("%2d:",i+1);
                if (fh->tbcoln) printf("\t%d",fh->tbcoln[i]);
                if (fh->tformn) printf("\t%-8s",fh->tformn[i]);
                /* optional: */
                if (fh->ttypen && fh->ttypen[i]) printf("\t%-8s",fh->ttypen[i]);
                if (fh->tunitn && fh->tunitn[i]) printf("\t%-8s",fh->tunitn[i]);
                if (fh->tnulln && fh->tnulln[i]) printf("\t%-8s",fh->tnulln[i]);
                if (fh->tzeron && fh->tzeron[i]) printf("\t%g",fh->tzeron[i]);
                if (fh->tscaln && fh->tscaln[i]) printf("\t%g",fh->tscaln[i]);
                printf("\n");
            }
        }
    }
    /* print additional items ?? */
    if (print!=NULL && *print!=NULL) {      /* if extra to print */
        printf("__________________________________________________________\n");
        printf("##\n");
        ncards = fh->hlen / FTSLINSIZ;
        for (i=0; i<ncards; i++) {          /* loop over all cards */
            card = &fts_buffer[i*FTSLINSIZ];       /* point to the card */
            for (sp=print; *sp; sp++) {         /* loop over all print`s */
                if (strncmp(card,*sp,strlen(*sp))==0) { /* see if key match */
                    my_copy(card,buf,FTSLINSIZ);    /* copy into temp buf */
                    buf[FTSLINSIZ] = 0;  /* make sure ends */
                    n = parse_card(-1,buf,a1,a2,a3,a4,a5);
                    if (n) {
                        if (a2[0]=='=')
                            printf("# %s = %s\n",a1,a3);
                        else if (a2[0]==' ')
                            printf("# %s %s\n",a1,a3);
                        else
                            printf("# %s\n",buf);       /* print whole */
                    }
                }
            }
        }
        printf("##\n");
    }
    dsize = fts_dsize(fh);
    printf("headersize = %d bytes = %d %d-records\n",
            fh->hlen, (fh->hlen - 1)/ftsblksiz_i + 1, ftsblksiz_i);
    printf("datasize = %ld bytes = %ld %d-records\n",
            dsize, (dsize-1)/ftsblksiz_i + 1, ftsblksiz_i);
    printf("__________________________________________________________\n");

    return 1;
}


/* 
 *  fts_read_img_coord:  get the WCS  in the cfitsio style
 *
 */

int fts_read_img_coord(
		       fits_header *fh,            /* (i)  pointer to fits header structure */
		       double *crval1, 
		       double *crval2,
		       double *crpix1, 
		       double *crpix2,
		       double *cdelt1, 
		       double *cdelt2,
		       double *rot,
		       char *proj)
{
  *crval1 = fh->crvaln[0];
  *crval2 = fh->crvaln[1];

  *crpix1 = fh->crpixn[0];
  *crpix2 = fh->crpixn[1];

  *cdelt1 = fh->cdeltn[0];
  *cdelt2 = fh->cdeltn[1];

}


static nnl = 0;
static void printnl(int fnl, int reset) {
  if (fnl==0) return;
  if (reset) {
    nnl=0;
    return;
  }
  nnl++;
  if (nnl%fnl == 0) printf("\n");
}

/* 
 *  fts_ptable:   print a table extensions:   XTENSION="TABLE   "
 *					      XTENSION="BINTABLE"
 *					      XTENSION="A3DTABLE"
 *                if not a table, it returns 0
 *
 */

int fts_ptable(
	       fits_header *fh,            /* (i)  pointer to fits header structure */
	       stream instr, 	           /* (i)  input stream data is associated with */
	       string *col,                /* (i)  optional selection of fields to print */
	       string select,              /* (i)  output mode: header and/or data */
	       int *row,                   /* (i)  list of rows to display; NULL or 1.. */
	       int fnl)                    /* (i)  frequency of a newline  */
{
    char *card, *dp, rowfmt[10];
    int ncards, i, j, k, n, w, len, pos, *colsel;
    int ncolin, ncolout, check, idx; 
    size_t nskip, ntail;
    bool colall, addrow, ascii;
    string *sp;

    if (scanopt(select,"group"))
        return fts_pgroup(fh, instr, col, select, row, fnl);

    if (fh->xtension == NULL) {
        fts_sdata(fh,instr);
        return 0;
    }
    if (!streq(fh->xtension,"TABLE")    &&      /* must be one of these */
        !streq(fh->xtension,"BINTABLE") &&
	!streq(fh->xtension,"A3DTABLE")) {
      fts_sdata(fh,instr);
      return 0;
    }
    ascii = streq(fh->xtension,"TABLE");
    if (!ascii) set_tbcoln(fh);                    /* patch tbcoln/tbitems */

    if (col == NULL || *col == NULL)
        colall = TRUE;
    else
        colall = FALSE;

    len = fh->naxisn[0];		         /* linelength of one row */
    ncards = fh->naxisn[1];	                    /* number of rows */
    ncolin = fh->tfields;
    ncolout = (colall ? ncolin : xstrlen(col,sizeof(string))-1);

    if (scanopt(select,"header") || scanopt(select,"info")) {
      printf("Table: %s rows: %d columns: %d\n",fh->extname,ncards,ncolin);
      if (scanopt(select,"header")) {
        if (fh->author)   printf("Author:    %s\n",fh->author);
	if (fh->referenc) printf("Reference: %s\n",fh->referenc);
        if (ascii) {
          printf("Fld Col Fmt     Name                 Units      Comments\n\n");
          for (i=0; i<fh->tfields; i++) {
            printf("%3d %3d %-7s %-20s %-10s %s\n", 
                    i+1,
                    fh->tbcoln[i],                         /* mandatory */
                    fh->tformn[i],                        /* mandatory */
                    (fh->ttypen ? 
                        (fh->ttypen[i]==NULL ? "-" : fh->ttypen[i]) : "-"),
                    (fh->tunitn ? 
                        (fh->tunitn[i]==NULL ? "-" : fh->tunitn[i]) : "-"),
                    (fh->ttypen_c ? 
                        (fh->ttypen_c[i]==NULL ? "-" : fh->ttypen_c[i]) : "-")
                  );

          }
        } else {
          printf("Fld Fmt     Name                 Units      Comments\n\n");
          for (i=0; i<fh->tfields; i++) {
            printf("%3d %-7s %-20s %-10s %s\n",
                    i+1,
                    fh->tformn[i],                      /* mandatory */
                    fh->ttypen[i],                      /* mandatory */
                    (fh->tunitn ?                       /* optional */
                        (fh->tunitn[i]==NULL ? "-" : fh->tunitn[i]) : "-"),
                    fh->ttypen_c[i]);                    /* private mandatory */
          }
        }
      }
    }

    if (scanopt(select,"data")) {

        addrow = scanopt(select,"row");   /* see if row number be output too? */
        if (addrow) {
            w = log10( (double)ncards ) + 1;
            sprintf(rowfmt,"%%%dd ",w);
            dprintf(1,"w = %d rowfmt = %s\n", w, rowfmt);
        }
        if (!colall) {
            colsel = (int *) allocate(ncolout * sizeof(int));
            check = 0;
            for (i=0; i<ncolout; i++) {  /* loop over all output columns requested */
                colsel[i] = 0;  /* set default to 'no column assigned' */
                for (j=0, sp=fh->ttypen; j<ncolin; j++, sp++) /* check all input cols */
                    if(streq(*sp,col[i])) {
                        colsel[i] = j+1;
                        break;
                    }
                check += colsel[i];
            }
            if(check==0) {
                warning("No columns matched: skipping data");
                fts_sdata(fh,instr);
                return 0;
            }
        }

        if (fh->ttypen) {
        printf("# ");                        /* print header with col names */
        pos = 3;
        for (i=0; i<fh->tfields; i++) {
            if (i>0 && ascii) {
                while(pos<fh->tbcoln[i]) {
                    printf(" ");
                    pos++;
                }
            }
            printf("%s ",fh->ttypen[i]);
            pos += strlen(fh->ttypen[i])+1;
        }
        printf("\n");
        }

        printf("#");                            /* print alignment points */
        pos = 2;
        if (ascii)
          for (i=0; i<fh->tfields; i++) {
            if (i>0) {
                while(pos<fh->tbcoln[i]) {
                    printf(" ");
                    pos++;
                }
            }
            printf("#");
            pos += 1;
          }
        printf("\n");

        card = (char *) allocate(len+1);      /* space for one row of data */
        card[len]='\0';                   /* for ascii: always terminate it */
        for (j=1, idx=0; j<=ncards; j++) {      /* loop over all rows */
            if ( (n=fread(card,1,len,instr)) != len) 
                error("fts_ptable: fread'ing only n=%d",n);
            if (row) {        /* if processing only certain rows */
                if(row[idx]>j) continue;
                if(row[idx]<j) continue;	/* read & break also ok */
                idx++;	                        /* found card to print: */
            }


	    /*
	     *  fix this binary table column selection stuff
	     *  NOTE: this also doesn't use the TDIM stuff.
	     */


            if (colall) {     /* if printing all columns for this row */
                if (addrow) printf(rowfmt,j);
                /* 'BUG': since card could contain NULL's or so: patch it */
                if (ascii) {
		  dprintf(0,"printing all columns in ascii\n");
                  for (i=0; i<len; i++) if (card[i]=='\0') card[i]=' ';
                  printf("%s\n",card);
                } else {
		  dprintf(0,"printing all %d columns from binary\n",fh->tfields);
                  for (i=0; i<fh->tfields; i++) {
                    dp = &card[fh->tbcoln[i]];
                    if (strrchr(fh->tformn[i],'E')) {
                      for (k=0; k<fh->tbitems[i]; k++) {
                        printf(" %f",show_e((float *)dp,k));
			printnl(fnl,0);
		      }
                    } else if (strrchr(fh->tformn[i],'D')) {
                      for (k=0; k<fh->tbitems[i]; k++) {
                        printf(" %lf",show_d((double *)dp,k));
			printnl(fnl,0);
		      }
                    } else if (strrchr(fh->tformn[i],'I')) {
                      for (k=0; k<fh->tbitems[i]; k++) {
                        printf(" %d",show_i((short *) dp,k));
			printnl(fnl,0);
		      }
                    } else if (strrchr(fh->tformn[i],'J')) {
                      for (k=0; k<fh->tbitems[i]; k++) {
                        printf(" %d",show_j((int *) dp,k));
			printnl(fnl,0);
		      }
                    } else if (strrchr(fh->tformn[i],'K')) {
                      for (k=0; k<fh->tbitems[i]; k++) {
                        printf(" %lld",show_k((int8 *) dp,k));     /* portabilty not ok  */
			printnl(fnl,0);
		      }
                    } else if (strrchr(fh->tformn[i],'A')) {
                      printf(" ");
                      for (k=0; k<fh->tbitems[i]; k++)
                        printf("%c",dp[k]);
                    } else if (strrchr(fh->tformn[i],'X')) {
                        ;
                    } else if (strrchr(fh->tformn[i],'L')) {
                        ;
                    } else
                      error ("BINTABLE %s not encoded format yet",fh->tformn[i]);
                  }
                  if (fnl==0) printf("\n");
                }
            } else {
                if (addrow) printf(rowfmt,j);
                for (i=0; i<ncolout; i++) {
                    k = colsel[i]-1;
                    if(k<0) continue;
                    printf("%s ",fmt(&card[fh->tbcoln[k]-1],fh->tformn[k]));
		    printnl(fnl,0);
                }
                if (fnl==0) printf("\n");
            }
        } /* for (j) */
        free(card);

        /* skip over trailing end of data */
        
        nskip = fts_dsize(fh);             
    	if (nskip % ftsblksiz_i) {
        	ntail = ftsblksiz_i - nskip % ftsblksiz_i;
        	nskip += ntail;
    	}
    	nskip -= ncards * len;		/* subtract what already read */
      	dprintf(1,"Skipping %ld bytes for tail\n",ntail);        
        if (nskip < 0) error("fts_ptable: negative tail %ld to skip",nskip);
    	fseek(instr, nskip, 1);
    } else 
        fts_sdata(fh,instr);

    return 1;
}

/*
 *  fts_pgroup: 	print random group parameters of a RG table
 */

int fts_pgroup(
fits_header *fh,            /* (i)  pointer to fits header structure */
stream instr,		    /* (i)  input stream data is associated with */
string *col,                /* (i)  optional selection of fields to print */
string select,              /* (i)  output mode: header and/or data */
int *row,                   /* (i)  list of rows to display; NULL or 1.. */
int fnl)                    /* (*)  frequency of a newline  [not implemented] */
{
    char *card, rowfmt[10];
    int ncards, i, j, k, n, w, len, nitems;
    int bytpix, idx; 
    size_t nskip, ntail;
    bool colall, addrow;
    real pscal, pzero;
    bool Qpar, Qval, Qraw;

    Qpar = scanopt(select,"group");
    Qval = scanopt(select,"data");
    Qraw = scanopt(select,"raw");
    
       
    dprintf(0,"Random Group printing:\n");

    if (fh->naxisn[0] == 0) {       /* NAXIS1=0 is for random groups */
        if (fh->groups != 1) 
            warning("FITS error: NAXIS1 = 0 and GROUPS not T");
    } else {
        fts_sdata(fh,instr);
        return 0;
    }

    if (col == NULL || *col == NULL)
        colall = TRUE;
    else
        colall = FALSE;

    for (i=1, len=1; i<fh->naxis; i++)		 /* compute ... */
        len *= fh->naxisn[i];
    len += fh->pcount;       		         /* linelength of one 'row' */
    bytpix = ABS(fh->bitpix)/8;
    nitems = len;
    len *= bytpix;
    
    ncards = fh->gcount;	                 /* number of rows */

    dprintf(0,"RANDOM GROUP: %d rows with %d (%d-byte) length each; %d pars\n",
		ncards,len/bytpix,bytpix,fh->pcount);

    printf("#");

    
    if (Qpar)
        for (i=0; i<fh->pcount; i++) {
            if (fh->ptypen)
                printf(" %s",fh->ptypen[i]);
            else
                printf(" PTYPE%d",i+1);
        }
    if (Qval)
        printf(" Data(1..)");

    
    printf("\n");
    
    addrow = scanopt(select,"row");   /* see if row number be output too? */
    if (addrow) {
        w = log10( (double)ncards ) + 1;
        sprintf(rowfmt,"%%%dd ",w);
        dprintf(1,"w = %d rowfmt = %s\n", w, rowfmt);
    }

    card = (char *) allocate(len+1);      /* space for one row of data */
    card[len]='\0';                   /* for ascii: always terminate it */

    dprintf(1,"FTELL at data: 0x%x\n",ftell(instr));
    for (j=1, idx=0; j<=ncards; j++) {      /* loop over all rows */
            if ( (n=fread(card,1,len,instr)) != len) 
                error("fts_pgroup: fread'ing only n=%d",n);
            if (row) {        /* if processing only certain rows */
                if(row[idx]>j) continue;
                if(row[idx]<j) continue;	/* read & break also ok */
                idx++;	                        /* found card to print: */
            }

            if (colall) {     /* if printing all columns for this row */
                if (addrow) printf(rowfmt,j);
                /* 'BUG': since card could contain NULL's or so: patch it */
                for (k=0; k< nitems; k++) {
                    if (k <  fh->pcount && !Qpar) continue;
                    if (k >= fh->pcount && !Qval) break;
                    if (k < fh->pcount) {
                        pscal = (fh->pscaln && !Qraw ? fh->pscaln[k] : 1.0);
                        pzero = (fh->pzeron && !Qraw ? fh->pzeron[k] : 0.0);
                    } else {
                        pscal = 1.0;    /* use BSCALE and BZERO instead */
                        pzero = 0.0;
                    }
                    printf(" ");
                    if (fh->bitpix == -32)
                        printf("%g",show_e((float *)card,k)*pscal+pzero);
                    else if (fh->bitpix == -64)
                        printf("%g",show_d((double *)card,k)*pscal+pzero);
                    else if (fh->bitpix == 64)
                        printf("%g **",show_k((int8 *)card,k)*pscal+pzero);
                    else if (fh->bitpix == 32)
                        printf("%g",show_j((int *)card,k)*pscal+pzero);
                    else if (fh->bitpix == 16)
                        printf("%g",show_i((short *)card,k)*pscal+pzero);
                    else
                      error ("RANDOM GROUPS not encoded format yet");
                }
                printf("\n");
            } else {
                error("Cannot print selected columns for GROUPs yet");
            }
    } /* for (j) */

    free(card);

        /* skip over trailing end of data */
        
    nskip = fts_dsize(fh);             
    if (nskip % ftsblksiz_i) {
        	ntail = ftsblksiz_i - nskip % ftsblksiz_i;
        	nskip += ntail;
    }
    nskip -= ncards * len;		/* subtract what already read */
    dprintf(1,"Skipping %ld bytes for tail\n",ntail);        
    if (nskip < 0) error("fts_ptable: negative tail %ld to skip",nskip);
    fseek(instr, nskip, 1);

    return 1;
}

/*
   patch columns 

*/


local void set_tbcoln(fits_header *fh)
{
    int i,k,n,pos;
    
    if (fh->tbcoln)  warning("Overwriting 'tbcoln' structure members");
    fh->tbcoln = (int *) allocate(fh->tfields * sizeof(int));
    if (fh->tbitems) warning("Overwriting 'tbitems' structure members");
    fh->tbitems = (int *) allocate(fh->tfields * sizeof(int));
    pos = 0;
    for (i=0; i<fh->tfields; i++) {
	k = colitems(fh->tformn[i]);
        fh->tbcoln[i] = pos;                /* column where data start */
        fh->tbitems[i] = k;                 /* number of items */
        if (strrchr(fh->tformn[i],'E')) {        /* single prec real */
          n = 4;
        } else if (strrchr(fh->tformn[i],'D')) { /* double prec real */
          n = 8;
        } else if (strrchr(fh->tformn[i],'I')) { /* signed short integer */
          n = 2;
        } else if (strrchr(fh->tformn[i],'J')) { /* signed long integer */
          n = 4;
        } else if (strrchr(fh->tformn[i],'K')) { /* signed long long integer */
          n = 8;
        } else if (strrchr(fh->tformn[i],'C')) { /* single complex real */
          n = 8;
        } else if (strrchr(fh->tformn[i],'M')) {  /* double complex real */
          n = 16;
        } else if (strrchr(fh->tformn[i],'P')) {  /* pointer - var.arr. */
	  n = 8;
        } else if (strrchr(fh->tformn[i],'L')) {  /* logical */
	  n = 1;
        } else if (strrchr(fh->tformn[i],'B')) {  /* byte */
	  n = 1;
        } else if (strrchr(fh->tformn[i],'A')) {  /* character */        
	  n = 1;
        } else if (strrchr(fh->tformn[i],'X')) {  /* bit */
	  n = 0;
        } else
          error ("BINTABLE with Unencoded TFORM%d = %s",i+1,fh->tformn[i]);
	if (n == 0)
	  pos += ROUNDUP(k/8,1);
        else
    	  pos += n*k;
    } /* for (i) all fields */
}


local int put_ivals (
		     char *key, char *val, char *fitsname,
		     int *nvals,
		     int **ivals)
{
    int ivalue, idx, *ip;

    ivalue = atoi(val);
    if (key[strlen(fitsname)]==' ' || key[strlen(fitsname)]=='\0') {
        if (*nvals >= 0)
            error("%s already given a dimension %d\n",key,*nvals);
        *nvals = ivalue;
        ip = *ivals = (int *) allocate((*nvals)*sizeof(int));
        for (idx=0; idx<(*nvals); idx++) ip[idx] = 0;
    } else {
        ip = *ivals;
        if (ip == NULL) {
            if (*nvals < 0)
                    error("%s has no dimension yet\n",key);
            else {
                ip = (int *) allocate((*nvals)*sizeof(int));
                *ivals = ip;
                for (idx=0; idx<(*nvals); idx++) ip[idx] = 0;
            }
	}
        idx = atoi(&key[strlen(fitsname)]) - 1;
        dprintf(3,"DEBUG> (%s) ip[idx = %d] = %d\n",key,idx,ivalue);
        if (idx < *nvals)
            ip[idx] = ivalue;
        else
            dprintf(0,"### Value for %s cannot be saved\n",key);
    }
    return 1;
}
       
local int put_dvals (
		     char *key, char *val, char *fitsname,
		     int *nvals,
		     double **dvals)
{
    int ivalue, idx;
    double dvalue, *fp;

    if (key[strlen(fitsname)]==' ' || key[strlen(fitsname)]=='\0') {
        if (*nvals >= 0)
            error("%s already given a dimension %d\n",key,*nvals);
        ivalue = atoi(val);
        *nvals = ivalue;
        fp = *dvals = (double *) allocate((*nvals)*sizeof(double));
        for (idx=0; idx<(*nvals); idx++) fp[idx] = 0.0;
    } else {
        dvalue = atof(val);
        fp = *dvals;
        if (fp == NULL) {
            if (*nvals < 0)
                    error("%s has no dimension yet\n",key);
            else {
                fp = (double *) allocate((*nvals)*sizeof(double));
                *dvals = fp;
                for (idx=0; idx<(*nvals); idx++) fp[idx] = 0.0;
            }
	}
        idx = atoi(&key[strlen(fitsname)]) - 1;
        dprintf(3,"DEBUG> (%s) fp[idx = %d] = %f\n",key,idx,dvalue);
        if (idx < *nvals)
            fp[idx] = dvalue;
        else
            dprintf(0,"### Value for %s cannot be saved\n",key);
    }
    return 1;
}

local int put_cvals (
		     char *key, char *val, char *fitsname,
		     int *nvals,
		     string **cvals)
{
    int ivalue, idx;
    string svalue, *sp;

    if (key[strlen(fitsname)]==' ' || key[strlen(fitsname)]=='\0') {
        if (*nvals >= 0)
            error("%s already given a dimension %d\n",key,*nvals);
        ivalue = atoi(val);
        *nvals = ivalue;
        sp = *cvals = (string *) allocate((*nvals)*sizeof(string));
        for (idx=0; idx<(*nvals); idx++)
            sp[idx] = NULL;
    } else {
        svalue = val;
        sp = *cvals;
        if (sp == NULL) {
            if (*nvals < 0)
                    error("%s has no dimension yet\n",key);
            else {
                sp = (string *) allocate((*nvals)*sizeof(string));
                *cvals = sp;
                for (idx=0; idx<(*nvals); idx++)
                    sp[idx] = NULL;
            }
	}
        idx = atoi(&key[strlen(fitsname)]) - 1;
        dprintf(3,"DEBUG> (%s) cp[idx = %d] = %s\n",key,idx,svalue);
        if (idx < *nvals)
            sp[idx] = atoa(svalue); 
        else
            dprintf(0,"### Value for %s cannot be saved\n",key);
    }
    return 1;
}

/* 
 * parse_card:    break a card image, pointed to by 'card' into several pieces:
 *  a1 will contain the keyword, if it's there. From position 1-8 in 'buf'.
 *  a2 will contain an '=' if the keyword has a value. From position 9-10.
 *  a3 will contain the value, if present (i.e. if '=') present, or the
 *      the remainder if card[8] was not an '='. If the value field was
 *      a character string, the quotes are stripped, trailing blanks are
 *      stripped, and the string is properly zero terminated.
 *  a4 will contain an optional value. Could happen if the keyword is complex
 *  a5 will contain the comment field, excluding the slash (/)
 *
 *  => It does not do a whole lot of error checking, that could be
 *     be improved if you want to force stricter FITS.
 *
 *  Returns 1 if OK     (always right now)
 *         -1 if error
 */

#define KEYVALID "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_"   /* see 6.1.2.1. */

local int parse_card (int icard, char *card, char *a1, char *a2, char *a3, char *a4, char *a5)
{
    int i;
    char *cp, *buf;

    buf = card;
    *a1 = *a2 = *a3 = *a4 = *a5 = 0;       /* empty all 5 components */
    if (strlen(buf) != FTSLINSIZ)
        dprintf(0,"### Card %d does not contain %d characters, len=%d\n",
                          icard, FTSLINSIZ, strlen(buf));

/* KEYWORD: in bytes 1-8, blankfilled if not all 8 characters used */
/*          we'll read all 8, but zero terminate the a1 string     */

    for (cp=a1, i=0; i<8; i++)  /* copy 8 characters from buf to a1 */
        *cp++ = *buf++;
    *cp = 0;             /* set a1[9] to zero to be sure */
    for (cp=a1, i=0; *cp!=0; cp++) {   /* check keyword in a1 again */
        if (*cp==' ') {     /* aha, keyword less than 8 chars ? */
            *cp = 0;     /* patch/truncate the a1 string */
            cp++;           /* go to next char */
            while (*cp) {   /* and check all non-NULL ones */
                if (*cp != ' ') {    /* they MUST be blank */
                    dprintf(0,"Keyword %s has embedded blanks\n",a1);
                    break;
                }
                cp++;       /* and keep looping until end of string */
            } /* while */
            break;          /* and get out of this for()-loop */
        } else {
            if (strchr(KEYVALID,*cp)==NULL)
                i++;        /* count number of invalid characters */
        }
    } /* a1 */
    if (i>0)
        dprintf(0,"### Keyword %s has %d invalid character(s).\n",a1,i);
            
/* VALUE INDICATOR: IN BYTES 9-10 contain either '= ' or '  ' or comment */
/*      if card[8] does not contain an '=', all text from card[8] on can */
/*      be stuffed into a3, if not, parsing will continue                */

    if (*buf != '=') {
        a2[0] = 0;                   /* zero out the a2 parameter */
        i = 8;                          /* count where we are; starts at i=0 */
        while (*buf == ' ') {           /* skip leading blanks again */
            *buf++;  i++;
        }  
        for (cp=a3; i<FTSLINSIZ; i++)  /* and copy into a3 */
            *cp++ = *buf++;
        return 1;
    } else {
        a2[0] = *buf++;               /* put '=' string in a2 */
        a2[1] = 0;                   /* terminate */
    }
/*  Now we expect a true value, and column 10 SHOULD be a blank */
    if (*buf == ' ')
        buf++;          /* ok, this is what we expected */
    else 
        dprintf(0,"### Card %d (keyword %s) does not have blank in col. 10\n",
                icard, a1);        

/* VALUE / COMMENT FIELD IN BYTES 11-80 */
    /* although FITS says strings to start with a quote in column 11,
       we are liberal here and first skip blanks */

    while (*buf == ' ') {         /* skip initial blanks (liberal) */
        buf++;
    }

    if (*buf == '\'') {         /* character variable: special treatment */
        if (buf-card != 11)
            dprintf(2,"### Card %d (%s) has a string starting in column %d\n",
                        icard,a1,buf-card+2);
        buf++;                              /* skip the quote */
        for (cp=a3; *buf!='\'' && i<FTSLINSIZ ;i++ ) { /* copy buf to a3 */
            *cp++ = *buf++;         /* ... but may not handle quotes well */
        }
        *cp = 0;             /* terminate a3 string */
        if ((int)strlen(a3) < 8)
            dprintf(2,"### Keyword %s = %s has length < 8",a1,a3);
#if 0
	/* VALGRIND unitialized for fully blank keyword values */
        for (i=strlen(a3)-1; i>=0, a3[i]==' '; i--)
            a3[i] = 0;       /* get rid of end spaces in string */
#endif
        if (strlen(a3)==0)
            dprintf(2,"### Keyword %s has effectively zero length\n",a1);
        a4[0] = 0;
        if (*buf != '\'') {     /* would be a strange situation */
            if (buf-card != FTSLINSIZ-1)
                dprintf(0,"### Card %d (%s) has strange string value\n",
                    icard,a1);
        } else
            buf++;              /* skip the ending quote */
    } else {                    /* assume, a regular thingo */
        /* remember *buf is now non-blank, i counted >= 10   */
        cp = a3;                /* copy rest into a3 */
        while (*buf != ' ' && *buf != ',' && *buf != 0) /* until sepr. found */
            *cp++ = *buf++;
        *cp = 0;             /* terminate a3 string properly */
	while (*buf == ' ' || *buf == ',')
	    buf++;
        if (*buf == '/') {   /* check if comment is now coming up */
            a4[0] = 0;
        } else {            /* complex number or something: in a4 */
            cp = a4;
            while (*buf != ' ' && *buf)    /* copy all non-white space */
                *cp++ = *buf++;             /* copy buf to a4 */ 
            *cp = 0;                 /* terminate a4 */
            while (*buf == ' ')     /* skip more whitespace */
                buf++;
        }
    } /* end of parsing value */
    while (*buf == ' ')         /* skip blanks before the comment */
        *buf++;
    if (*buf != '/') {            /* if it is not a comment designator, quit */
        if (buf-card != FTSLINSIZ)
            dprintf(2,"### No comment or ??? in card %d, pos=%d\n",
                icard,buf-card);
        return(1);
    } 
    buf++;
    while (*buf == ' ')     /* skip more blanks - FITS strongly advices 1 */
        buf++;    
    for (cp=a5; buf-card < FTSLINSIZ; )
        *cp++ = *buf++;
    *cp = 0;             /* and terminate the a5 comment string */
    return 1;               /* always 1 for now */
}        


/* 
 * fts_buf:   check if the fits buffer still large enough, else
 *            allocate some extra memory.
 *            This buffer is a direct (byte) image of the fits 
 */

int fts_buf(
	    int len)        /*  (i)  lenght to test for */        
{


    if (fts_buffer==NULL) {
        fts_buflen = ftsblksiz_i;
        fts_buffer = allocate(ftsblksiz_i);
    } else if (len > fts_buflen) {
        fts_buflen += ftsblksiz_i;
        dprintf(4,"******* REALLOC BUFFER to %d bytes ******\n",fts_buflen);
        fts_buffer = reallocate(fts_buffer, fts_buflen);
    }
    return 1;    
}

/*
 *  fts_sdata:  skip the data from current file pointer
 *		29-jun-94: also takes into account data read (fh->nread)
 */

int fts_sdata(
	      fits_header *fh,         /*  (i) pointer to fits header */
	      stream instr)                   /*  (i) file pointer    */
{
    size_t nskip, ntail;

    nskip = fts_dsize(fh);
    dprintf(1,"Skipping %ld bytes for data (but %ld read)\n",
		nskip, fh->nread);

    if (nskip % ftsblksiz_i) {
        ntail = ftsblksiz_i - nskip % ftsblksiz_i;
        dprintf(1,"Skipping %ld bytes for empty tail\n",ntail);        
        nskip += ntail;
    }
    nskip -= fh->nread;
    if (nskip>0) fseek(instr, nskip, 1);
    return 1;
}

/*
 *  fts_cdata:  copy the data from current file pointer
 */
#define CONVBUFLEN   2880*10

int fts_cdata(
	      fits_header *fh,                /*  (i) pointer to fits header */
	      stream instr,           /* (i) file pointer for input */
	      stream outstr,           /* (i) file pointer  for output  */
	      bool trailr,			/* (i) need to read trailing end too ? */
	      bool trailw)			/* (i) need to write trailing end too ? */
{
    int  nread, nwrite, n, ntowrite, ntoread, itemlen, nitems;
    char buffer[CONVBUFLEN];

    itemlen = ABS(fh->bitpix)/8;        /* # bytes in an item */
    n = fts_dsize(fh);                  /* real data size  */
    if (trailr) {
      ntoread  = ROUNDUP(n,ftsblksiz_i);  /* total size to read though */
    } else
      ntoread = n;
    if (trailw) {
      ntowrite = ROUNDUP(n,ftsblksiz_o);  /* .. and to write */
    } else
      ntowrite = n;
    if (n==0) ntoread=0;                /* no data? */
    if (fh->flip) dprintf(1,"fts_cdata: Swapping bytes\n");
    dprintf(2,"fts_data: ntoread=%d   ntowrite=%d\n",ntoread,ntowrite);
    while (ntoread > 0) {
        if (ntoread > CONVBUFLEN)
            n = CONVBUFLEN;
        else
            n = ntoread;
        nread = fread(buffer,sizeof(char),n,instr);
        if (nread != n) /* should be a warning in raw mode */
            error("Tried to read %d, could only read %d",n,nread);
	fh->nread += nread;
	if (fh->flip) {
            nitems = nread/itemlen;
            if (nread%itemlen) error("Bad bufsize %d for byteswap",nread);
            bswap(buffer,itemlen, nitems);
        }
/****FIX****/
        nwrite = fwrite(buffer,sizeof(char),n,outstr);
        if (nwrite != n) 
            error("Tried to write %d, could only write %d\n",n,nwrite);
        ntoread -= nread;
    }
    return 1;
}

int fts_cdata816(
	      fits_header *fh,                /*  (i) pointer to fits header */
	      stream instr,           /* (i) file pointer for input */
	      stream outstr)           /* (i) file pointer  for output  */
{
    int  i,nread, nwrite, n, ntowrite, ntoread, itemlen, nitems;
    char buffer[CONVBUFLEN];
    char zero[1];

    if (fh->bitpix != 8) error("fts_cdata816: illegal bitpix=%d",fh->bitpix);

    n = fts_dsize(fh);                  /* real data size  */
    if (n == 0) return 1;

    zero[0] = 0;
    


    if (n==0) ntoread=0;                /* no data? */
    dprintf(2,"fts_cdata816: ntoread=%d   ntowrite=%d\n",ntoread,ntowrite);
    
    for (i=0; i<n; i++) {
        nread = fread(buffer,sizeof(char),1,instr);
        if (nread != 1) error("Tried to read 1, could only read %d",nread);
	fh->nread += nread;
        nwrite = fwrite(zero,sizeof(char),1,outstr);
        if (nwrite != 1) error("Tried to write 1, could only write %d\n",nwrite);
        nwrite = fwrite(buffer,sizeof(char),1,outstr);
        if (nwrite != 1) error("Tried to write 1, could only write %d\n",nwrite);
    }
    fh->bitpix = 16;
    n *= 2;
    ntowrite = ROUNDUP(n,ftsblksiz_o);  /* .. and to write */
    if (n == ntowrite) return 1;
    n = ntowrite-n;
    for (i=0; i<n; i++) {
      nwrite = fwrite(zero,sizeof(char),1,outstr);
      if (nwrite != 1) error("Tried to write 1, could only write %d\n",nwrite);
    }
    return 1;
}

/*
 *  fts_dsize: return size of data portion (without trailing junk) 
 *             of the fits file. You need to round this up to multiples
 *	       2880 bytes to get the full size.
 *             See also eq. (x.y.z) in the NOST manual
 */

int fts_dsize(fits_header *fh)
{
    int i, size;

    if (fh->naxis > 0) {
        if (fh->naxisn[0] == 0) {       /* NAXIS1=0 is for random groups */
            if (fh->groups != 1)
                warning("FITS error: NAXIS1 = 0 and GROUPS not T");
            size = 1;                   /* must be a group then */
        } else {                        /* else 'regular' good old fits */
            size = fh->naxisn[0];
        }
        for (i=1; i < fh->naxis; i++)   /* accumulate the other axes */
            size *= fh->naxisn[i];
        size += fh->pcount;             /* add the random parameters */
        size *= fh->gcount * ABS(fh->bitpix) / 8;
    } else
        size = 0;
    return size;
}

/*
 *  fts_tsize: return trailing size of useable data portion
 */

int fts_tsize(fits_header *fh)
{
  int n = fts_dsize(fh);
  int tail = ROUNDUP(n,ftsblksiz_o);
  return tail - n;
}

/*
 *  fts_rdata:  read the data into a (double precision) buffer
 *      return 0 if not all the data has been read  (due to buffering)
 *             1 if all data have been read
 *           < 0 some kind of error
 *      This routine also swaps the incoming bytes if the magical
 *      fh->flip is set.... (See DECORDER)
 */

int fts_rdata(
	      fits_header *fh,
	      stream instr,
	      int buflen,
	      char *buf)
{
    double *outp, convbuf[CONVBUFLEN];   /* local buffer used for conversion */
    size_t nread, nskip;      /* number of bytes to read */
    int i, *ip;
    int8 *kp;
    short *sp;
    double *dp;
    float *fp;

    nskip = nread = fts_dsize(fh);      /* total data count (bytes) */
    dprintf(1,"Attempting to read %d data bytes\n",nread);
    if (buflen < nread) {
        dprintf(0,"Buffer not large enough to store data\n");
        return -1;
    }
    if (CONVBUFLEN*sizeof(double) < nread) {
        dprintf(0,"DATA buffer too small - data skipped\n");
        return -1;
    }

    outp = (double *) buf;    /* convert to double's for now */
    fh->nread += nread; 
    switch (fh->bitpix) {       /* read another local buffer */
    case 16:
        nread /= 2;
        sp = (short *) convbuf;
        i = fread(sp,2,nread,instr);
        if (i!=nread) error("Did not read the right amount");
        if (fh->flip) bswap(sp,2,nread);
        for (i=0; i<nread; i++)
            *outp++ = sp[i] * fh->bscale + fh->bzero;
        break;
    case 32:
        nread /= 4;
        ip = (int *) convbuf;     /* local buffer */
        i = fread(ip,4,nread,instr);  /* straight read */
        if (i != nread) error("Did not read right amount");
        if (fh->flip) bswap(sp,4,nread);
        for (i=0; i<nread; i++) 
            *outp++ = ip[i] * fh->bscale + fh->bzero;
        break;
    case 64:
        nread /= 8;
        kp = (int8 *) convbuf;     /* local buffer */
        i = fread(kp,8,nread,instr);  /* straight read */
        if (i != nread) error("Did not read right amount");
        if (fh->flip) bswap(sp,8,nread);
        for (i=0; i<nread; i++) 
            *outp++ = ip[i] * fh->bscale + fh->bzero;
        break;
    case -32:
        nread /= 4;
        fp = (float *) convbuf;
        i = fread(fp,4,nread,instr);  /* straight read */
        if (i != nread) error("Did not read right amount");
        if (fh->flip) bswap(sp,4,nread);
        for (i=0; i<nread; i++) 
            *outp++ = fp[i] * fh->bscale + fh->bzero;
        break;
    case -64:
        nread /= 8;
        dp = (double *) convbuf;
        i = fread(dp,8,nread,instr);  /* straight read */
        if (i != nread) error("Did not read right amount");
        if (fh->flip) bswap(sp,8,nread);
        for (i=0; i<nread; i++) 
            *outp++ = dp[i] * fh->bscale + fh->bzero;
        break;
    default:
        error("Illegal bitpix=%d\n",fh->bitpix);
        return -1;
    }

    fts_sdata(fh,instr);        /* skip remaining (tail) part of data */
    return 0;
}

#if 0

    if (nskip % ftsblksiz_i == 0) /* no need to read tail */
        return 0;

    ntail = ftsblksiz_i - nskip % ftsblksiz_i;
    dprintf(1,"Reading %ld tail bytes\n",ntail);
    for (i=0; i<ntail; i++) {
        if (fread(&c,1,1,instr) != 1) {
            dprintf(0,"*** Error reading tail byte number %d\n",i);
            return -1;
        }
        /* in principle byteswap could occur here ..... not done though */
    }
    fh->nread += ntail;
    return 0;
}
#endif

int fts_rrow(
  fits_header *fh,
  stream instr,
  int nread,
  char *buf)
{
    int n;

    dprintf(1,"fts_rrow: at %d reading %d\n",fh->nread,nread);
    fh->nread += nread;
    n = fread(buf,1,nread,instr);
    if (n!=nread)
        error("fts_rrow: (%d/%d) did not read the right amount",n,nread);

    return 0;
}

int fts_wrow(
  fits_header *fh,
  stream instr,
  int nwrite,
  char *buf)
{
    int n;

    dprintf(1,"fts_wrow: at %d writing %d\n",fh->nwritten,nwrite);
    fh->nwritten += nwrite;
    n = fread(buf,1,nwrite,instr);
    if (n!=nwrite)
        error("fts_wrow: (%d/%d) did not read the right amount",n,nwrite);

    return 0;
}

/*
 *  fts_zero:   reset contents of a fits structure
 *              returns 1 on success, 0 on failure
 */

int fts_zero(
fits_header *fh)     /* (i/o) pointer to fits structure */
{
    if (fh==NULL) return 0;

    fh->hlen     = fh->dlen    = 0;
    fh->nread    = 0;
    fh->nwritten = 0;
    fh->flip     = 0;
    fh->ascii    = 0;
    fh->simple   = -1;        /* -1 or NULL normally signal absence */
    fh->groups   = -1;
    fh->extend   = -1;
    fh->bitpix   = -1;
    fh->naxis    = -1;        fh->maxis  = -1;
    fh->naxisn   = NULL;      fh->maxisn = NULL;
    fh->crvaln   = NULL;
    fh->cdeltn   = NULL;
    fh->ctypen   = NULL;
    fh->cunitn   = NULL;
    fh->crotan   = NULL;
    fh->history  = NULL;
    fh->comment  = NULL;
    fh->xtension = NULL;
    fh->extname  = NULL;
    fh->extver   = 1;
    fh->extlevel = 1;
    fh->bzero    = 0.0;
    fh->bscale   = 1.0;
    fh->blocked  = -1;
    fh->telescop = NULL;
    fh->instrume = NULL;
    fh->crpixn   = NULL;
    fh->pcount   = 0;
    fh->gcount   = 1;
    fh->pzeron   = NULL;          
    fh->pscaln   = NULL;
    fh->ptypen   = NULL;
    fh->punitn   = NULL;
    fh->tfields  = -1;
    fh->ttypen   = NULL;
    fh->ttypen_c = NULL;
    fh->tformn   = NULL;
    fh->tdispn   = NULL;
    fh->tunitn   = NULL;
    fh->tnulln   = NULL;
    fh->tbcoln   = NULL;
    fh->tscaln   = NULL;
    fh->tzeron   = NULL;
    fh->tbitems  = NULL;
    fh->tdimn    = NULL;
    fh->obsra    = 0.0;
    fh->obsdec   = 0.0;

    return 1;
}

/*=======================================================================*/
/*  Local Utilities for the FITS I/O                                     */
/*
 *  FITS keyword must be part of an array, e.g. CDELTx, NAXISn, ...
 *
 */

/* 
 *  atob: convert a FITS logical to an integer (0 or 1)
 *  The rules are over-liberal. No check is done if the T/F is
 *  in the right column (31?)
 *
 *      T (true)  -> 1
 *      F (false) -> 0
 */

local int atob(
char *cp)               /* (i) pointer to character string to convert */
{
    while (*cp == ' ')          /* skip spaces (no tabs allowed) */
        cp++;
    if (*cp == 0)            
        error("End of string where logical (T/F) expected");
    if (*cp == 'T')
        return(1);
    else if (*cp == 'F')
        return(0);
    else
        error("Logical (T/F) expected, but found %s\n",cp);
    return(0);      /* fool compiler, since it will (should) never get here */
}

/*
 *  atoa:   copy a FITS string (normally things between quotes)
 *
 */

local char *atoa(char *s)
{
    char *cp;
    
    cp = allocate(strlen(s)+1);
    strcpy(cp,s);
    return cp;
}

#if 0
local int get_array(buf,arr,narr)
char *buf;
double *arr;
int narr;
{
  char *cp;
  int i;
  
  if ((cp = strpbrk(buf,"0123456789")) == NULL)
    error("FITS keyword not an array: (%s)\n",buf);
  
  i = atoi(cp)-1;     /* point to array location */
}
#endif

local void my_copy (char *src, char *dest, int n)
{
    while (n--)  *dest++ = *src++;
}



local char fline[2*FTSLINSIZ];  /* extra buffer length for screwups */
local int flinetot = 0;         /* keep track how many lines written */

int fts_whead(fits_header *fh, stream ostr)
{
    int i;

    flinetot = 0;                    /* reset line counter */
    if (fh->simple < 0) {
        if (fh->xtension != NULL)
            fts_wvarc(ostr,"XTENSION",fh->xtension,NULL);
        else
            error("Header structure error: SIMPLE and XTENSION missing\n");
    } else
        fts_wvarb(ostr,"SIMPLE",fh->simple,NULL);
    fts_wvari  (ostr,"BITPIX",fh->bitpix,NULL);
    fts_wvari  (ostr,"NAXIS",fh->naxis,NULL);
    fts_wvari_a(ostr,"NAXIS",fh->naxis,fh->naxisn,NULL);
    if (fh->extend >= 0)
      fts_wvarb  (ostr,"EXTEND",fh->extend,NULL);
    if (fh->extname)
      fts_wvarc  (ostr,"EXTNAME",fh->extname,NULL);
    if (fh->groups >= 0) {
        fts_wvarb  (ostr,"GROUPS",fh->groups,NULL);
        fts_wvari  (ostr,"PCOUNT",fh->pcount,NULL);
        fts_wvari  (ostr,"GCOUNT",fh->gcount,NULL);
        if (fh->ptypen)
            fts_wvarc_a(ostr,"PTYPE",fh->pcount,fh->ptypen,NULL);
    }
    if (fh->crpixn)
        fts_wvard_a(ostr,"CRPIX",fh->naxis,fh->crpixn,NULL);
    if (fh->crvaln)
        fts_wvard_a(ostr,"CRVAL",fh->naxis,fh->crvaln,NULL);
    if (fh->cdeltn)
        fts_wvard_a(ostr,"CDELT",fh->naxis,fh->cdeltn,NULL);
    if (fh->crotan)
        fts_wvard_a(ostr,"CROTA",fh->naxis,fh->crotan,NULL);
    if (fh->ctypen)
        fts_wvarc_a(ostr,"CTYPE",fh->naxis,fh->ctypen,NULL);

    fts_wvard(ostr,"DATAMIN",fh->datamin,NULL);
    fts_wvard(ostr,"DATAMAX",fh->datamax,NULL);
    fts_wvard(ostr,"BSCALE",fh->bscale,NULL);
    fts_wvard(ostr,"BZERO",fh->bzero,NULL);
    fts_wvar(ostr,"COMMENT",cfits1);
    fts_wvar(ostr,"COMMENT",cfits2);
    fts_wvar(ostr,"COMMENT",cfits3);
    if (fh->comment) 
        for (i=0; fh->comment[i] != NULL; i++)
            fts_wvar(ostr,"COMMENT",fh->comment[i]);
    if (fh->history)
        for (i=0; fh->history[i] != NULL; i++)
            fts_wvar(ostr,"HISTORY",fh->history[i]);
    fts_wvar   (ostr,"END",NULL);

    while (flinetot % ftslpb_o)        /* while not at end of a '2880' block */
        fts_wvar(ostr," ",NULL);        /* ... fill with blank lines */
    return 1;
}


/* 
 *   poutline:  padd the fline[] string with blanks 
 *              and send it to the output file stream
 */

local void poutline(FILE *fp)
{
    int i;

    for (i=strlen(fline); i<FTSLINSIZ; i++)     /* padd blanks */
        fline[i] = ' ';
    fwrite(fline,1,FTSLINSIZ,fp);
    flinetot++;         /* poutline() BUG BUG::: should be in struct FITS */
}

int fts_wvar(                  /*  write simple variable, no '=' */
	     stream fp,
	     string key,
	     string comment)
{
    sprintf(fline,"%-8s ",key);       /* keyword takes up 8, and a space */
    if (comment)                    /* comment starts on 10, if present */
        strcat(fline,comment);    
    poutline(fp);
    return 1;
}

int fts_wvarc(    /* write string variable */
	      stream fp,
	      string key,
	      string value,
	      string comment)
{
    if ((int)strlen(value) > 8)
        dprintf(0,"### Character value for key %s (%s) has length > 8\n",
                        key,value);
    sprintf(fline,"%-8s= \'%s\' ",key,value);       /* value at %-8s or %s ??? */
    if (comment) {
        strcat(fline,"/ ");
        strcat(fline,comment);
    }
    poutline(fp);
    return 1;
}

int fts_wvarc_a(    /* write string variable array */
		stream fp,
		string key,
		int n,
		string *values,
		string comment)
{
    int i;
    char keyi[10];

    if ((int)strlen(key)>5)
        error("Length of array keyword %s must be less then 5\n",key);

    for (i=0; i<n; i++) {
        if (values[i]) {
            sprintf(keyi,"%s%d",key,i+1);
            fts_wvarc(fp,keyi,values[i],comment);
        }
    }
    return 1;
}


int fts_wvarb(   /* write single logical variable */
	      stream fp,
	      string key,
	      int value,
	      string comment)
{
    sprintf(fline,"%-8s=                    %c ",key,(value==1)? 'T' : 'F');
    if (comment) {
        strcat(fline,"/ ");
        strcat(fline,comment);
    }
    poutline(fp);
    return 1;
}

int fts_wvarb_a(    /* write logical variable array */
		stream fp,
		string key,
		int n,
		int *values,
		string comment)
{
    int i;
    char keyi[10];

    if ((int)strlen(key)>5)
        error("Length of array keyword %s must be less then 5\n",key);

    for (i=0; i<n; i++) {
        sprintf(keyi,"%s%d",key,i+1);
        fts_wvarb(fp,keyi,values[i],comment);
    }
    return 1;
}


int fts_wvari(    /* write single int variable */
	      stream fp,
	      string key,
	      int value,
	      string comment)
{
    char stmp[21];

    sprintf(stmp,"%d",value);
    sprintf(fline,"%-8s= %20s ",key,stmp);
    if (comment) {
        strcat(fline,"/ ");
        strcat(fline,comment);
    }
    poutline(fp);
    return 1;
}

int fts_wvari_a(   /* write int variable array */
		stream fp,
		string key,
		int n,
		int *values,
		string comment)
{
    int i;
    char keyi[10];

    if ((int)strlen(key)>5)
        error("Length of array keyword %s must be less then 5\n",key);

    for (i=0; i<n; i++) {
        sprintf(keyi,"%s%d",key,i+1);
        fts_wvari(fp,keyi,values[i],comment);
    }
    return 1;
}
   /* write single float variable */
int fts_wvarf(stream fp, 
	      string key, 
	      float value, 
	      string comment)
{
    char stmp[21];

    sprintf(stmp,"%g",value);
    sprintf(fline,"%-8s= %20s ",key,stmp);
    if (comment) {
        strcat(fline,"/ ");
        strcat(fline,comment);
    }
    poutline(fp);
    return 1;
}

int fts_wvarf_a(    /* write float variable array */
		stream fp,
		string key,
		int n,
		float *values,
		string comment)
{
    int i;
    char keyi[10];

    if ((int)strlen(key)>5)
        error("Length of array keyword %s must be less then 5\n",key);

    for (i=0; i<n; i++) {
        sprintf(keyi,"%s%d",key,i+1);
        fts_wvarf(fp,keyi,values[i],comment);
    }
    return 1;
}


int fts_wvard(   /* write single double variable */
	      stream fp,
	      string key,
	      double value,
	      string comment)
{
    char stmp[21];
    sprintf(stmp,"%18.11g",value);
    sprintf(fline,"%-8s= %20s ",key,stmp);
    if (comment) {
        strcat(fline,"/ ");
        strcat(fline,comment);
    }
    poutline(fp);
    return 1;
}

int fts_wvard_a(   /* write double variable array */
		stream fp,
		string key,
		int n,
		double *values,
		string comment)
{
    int i, ia0, ia1;
    char keyi[10];

    if ((int)strlen(key)>5)
        error("Length of array keyword %s must be less then 5\n",key);

    if (streq(key,"CROTA")) {
      /* notice: should really check which of the two are sky angles */
      if (n>=2) {
	ia0 = 0;
	ia1 = 1;
	if (values[ia0]!=0 && values[ia1]==0) {
	  warning("setting CROTA2=%g since it was missing",values[ia0]);
	  values[ia1] = values[ia0];
	} else if (values[ia1]!=0 && values[ia0]==0) {
	  warning("setting CROTA1=%g since it was missing",values[ia1]);
	  values[ia0] = values[ia1];
	}
      }
    }

    for (i=0; i<n; i++) {
        sprintf(keyi,"%s%d",key,i+1);
        fts_wvard(fp,keyi,values[i],comment);
    }
    return 1;
}


int fts_wdata(
	      fits_header *fh,
	      stream ostr,
	      int n,
	      char *cp)
{
    if (fwrite(cp,sizeof(char),n,ostr) != n)
        error("Error writing %d data bytes",n);
    fh->nwritten += n;
    dprintf(4,"Accumulated %d written bytes so far\n",fh->nwritten);
    return 1;
}

/*
 * check_date:  check if a string of the format  'dd/mm/yy'
 *		formally dd-mm-yy is not allowed
 * as per 1-jan-1998 the new ISO standard may go into effect
 * this format is
 *        ccyy-mm-dd
 *        ccyy-mm-ddThh:mm:ss.ss
 * with an optional time indication
 */

local void check_date(char *key, char *s)
{
    int err=0, old=0, new=0;
    char *val;

    val = s;
    if ((int)strlen(s) < 8)
        err++;
    else {
       if (!isdigit(*s++)) err++;       /*  s[0]  */
       if (!isdigit(*s++)) err++;       /*  s[1]  */
       if (*s == '/') {                 /*  s[2]  */
          old++;
          s++;
          if (!isdigit(*s++)) err++;       /*  s[3]  */
          if (!isdigit(*s++)) err++;       /*  s[5]  */
          if (*s++ != '/') err++;
          if (!isdigit(*s++)) err++;       /*  s[7]  */
          if (!isdigit(*s++)) err++;       /*  s[8]  */
       } else if (isdigit(*s)) {
          new++;
          s++;
          if (!isdigit(*s++)) err++;
          if (*s++ != '-') err++;
          if (!isdigit(*s++)) err++;
          if (!isdigit(*s++)) err++;
          if (*s++ != '-') err++;
          if (!isdigit(*s++)) err++;
          if (!isdigit(*s++)) err++;
          if (*s) {
             if (*s++ != 'T') err++;
             if (!isdigit(*s++)) err++;
             if (!isdigit(*s++)) err++;
             if (*s++ != ':') err++;             
             if (!isdigit(*s++)) err++;
             if (!isdigit(*s++)) err++;
             if (*s++ != ':') err++;             
             if (!isdigit(*s++)) err++;
             if (!isdigit(*s++)) err++;
          }
       } else
          err++;
    }
    if (err || (old>0 && new>0))
        dprintf(0,"### Keyword %s has non-standard date structure (%s)\n",
                key,val);
}

/* 
 *  check_ctype: check if string is a known ctype...
 */

local char *ctype0[] = { "RA", "DEC", "GLON", "GLAT", "ELON", "ELAT", 
                         "DX", "DY",
                          "VELO", "FELO", "FREQ", "STOKES", "LL", "MM",
                          ".   ", NULL};
local char *ctype5[] = { "TAN", "SIN", "ARC", "NCP", "STG", "AIT", "GLS", "CAR",
			  "MER", "LSR", "HEL", "OBS", "ATF",/* ATF to be conmfirmed by */
                          ".   ", NULL};

local int check_ctype(char *key, char *val)
{
    int i;

    if ((int)strlen(val) < 4) {     
        dprintf(2,"### check_ctype: Cannot check %s = %s\n",key,val);
        return 0;
    }
    if ((int)strlen(val) >= 4) {         /* check first part */
        for (i=0; ctype0[i] != NULL; i++) {
            if (strncmp(&val[0],ctype0[i],strlen(ctype0[i]))==0)
                break;              /* found match - quit loop */
        }
        if (ctype0[i] == NULL)
            dprintf(2,"### Warning: %s = %s unknown CTYPE(1:4)\n",key,val);
    } else
        dprintf(0,"### Warning: %s = %s cannot check CTYPE(1:4)\n",key,val);
    if ((int)strlen(val) >= 8) {         /* check second part */
        for (i=0; ctype5[i] != NULL; i++) {
            if (strncmp(&val[5],ctype5[i],strlen(ctype5[i]))==0)
                break;              /* found match - quit loop */
        }
        if (ctype5[i] == NULL)
            dprintf(2,"### Warning: %s = %s unknown CTYPE(5:8)\n",key,val);
    } else
        dprintf(2,"### Warning: %s = %s cannot check CTYPE(5:8)\n",key,val);
    return 1;
}

local char *units[] = { 
    "DEGREES", "RADIANS", "arcsec", "arcmin", "deg", "rad",      /* angles */
    "W.U.", "JY", "KELVIN", "mag",                              /* brightness */
    "ANGSTROMS", "micron",                                     /* wavelenght */
    "MHZ", "KHZ", "HZ",                                        /* frequency */
    "SQSEC",                                                    /* areas */
    "Mpc", "kpc", "pc", "au", "km", "m", "cm", "mm",            /* length */
    "year", "yr", "month", "day", "min", "sec", "s",            /* time */
    "Msol", "kg", "g",                                          /* mass */
    "erg", "W",                                                 /* energy */
    "counts",                                                   /* numbers */
    NULL};


local int check_unit(char *key, char *val)
{
    int i;

    /* first check the unit as a whole */
    for (i=0; units[i] != NULL; i++) {
        if (strncmp(val,units[i],strlen(units[i]))==0)
                return 1;              /* found match - quit loop */

    } 

    /* since no match found so far, break it up, using '/' and '*' as
     * field separators, and try again 
     *  ==> still to be implemented
     */


    /* code arrives here if all else fails, i.e. units are not known
     * to the program 
     */
    dprintf(1,"### Warning: %s = %s unknown units\n",key,val);
    return 0;
}

int fts_setiblk(int n)
{
  if (n<=0) error("fts_setiblk: Illegal blocking factor %d\n",n);
  ftsblksiz_i=n * FTSBLKSIZ;	/* new blocksize for input */
  ftslpb_i=n * FTSLPB;		/* new lines per block for input */
  dprintf(1,"Input Blocking factor %d set",n);
  return n;
}

int fts_setoblk(int n)
{
  if (n<=0) error("fts_setoblk: Illegal blocking factor %d\n",n);
  ftsblksiz_o=n * FTSBLKSIZ;	/* new blocksize for output */
  ftslpb_o=n * FTSLPB;		/* new lines per block for output */
  dprintf(1,"Output Blocking factor %d set",n);
  if (ftsblksiz_i != ftsblksiz_o) 
    warning("Cannot handle different values for blocksize i/o");
  if (ftsblksiz_i > ftsblksiz_o) 
    warning("but ok since blocksize(input) > blocksize(output)");
  return n;
}


local int my_findstr(string text, string pat, int len)
{
  register string s;
  int nch;
        
  nch = strlen(pat);

  for (s = text; *s && s-text < len; s++)
    if (strncmp(s, pat, nch) == 0)
       return (s - text);
  return -1;

}

local void blank_fill(char *card, int len)
{
    int i;

    for (i=0; i<len; i++)    /* do not use strlen, since there may be no 0 */
        if (card[i] == 0) break;
    while (i<len)           /* right fill with blanks if need be */
        card[i++] = ' ';
}
                                        
/* local routine, made for fts_ptable: probably useful in library */

local int colmask(
	    int n,              /* (i) number of keys */
	    string key[],       /* (i) array of 'n' keys to check against */
	    string col[],       /* NULL terminated array of strings to check for */
	    int colsel[])        /* (o) array of 'n' int's : 0=not selected 1=selected */
{
    string *s;
    int i;
    
    if (col==NULL || *col==NULL) {
        for(i=0; i<n; i++) colsel[i] = 1;
        return 1;
    }

    for(i=0; i<n; i++) colsel[i] = 0;

    for (s=col; *s; s++) {
        for (i=0; i<n; i++)
            if (colsel[i]==0) colsel[i] = streq(*s,key[i]);
    }
    
    for (i=0; i<n; i++)
        if (colsel[i]) return 1;
    return 0;
}


/*
 * string to integer conversion, primarely intended for decoding
 * the TFORMnnn's of BINTABLE's
 */

local int colitems(char *s)
{
    int len;
    char *num, *cp;

    if (s==NULL || *s==0) return 0;

    len = strlen(s);
    num = (char *) allocate(len+1);
    strcpy(num,s);
    cp = num;
    while (*cp && isdigit(*cp))
        cp++;
    *cp = 0;
    if (*num)
        return atoi(num);   /* format like: '2E' with preceding count */
    else
        return 1;           /* format like:  'E' no preceding count */
    
}
/*
 *  fmt: crude fortran FMT to string conversion: 
 *	It handles things like A20, F20.10, D12.7 etc.
 *	It will not handle repeat keys, such as 2X or 3F20.10
 *
 *       NOTE: returns space to statically allocate space!!!
 */

local char *fmt(char *line, char *format)
{
    char *cp = format;
    static char sample[80];
    int n;

#if 1
    /* skip numerical count */   
    while (*cp && strchr("0123456789",*cp))
        cp++;
#endif

    switch (*cp) {
      case 'A':
      case 'I':
      case 'F':
      case 'E':
      case 'D':
        cp++;
        break;
      default:
        warning("Unrecognized format: %s",format);
        cp++;
    }
    n = atoi(cp);       /* will also correctly parse '4.1' as '4' */
    strncpy(sample,line,n);
    sample[n] = '\0';
    return sample;
}

/*
 * return atomic entities in native format
 * for machines which require things to be aligned, more work needs to be done
 */

local float show_e(float *v, int i)
{
#if 1
    float f;
    int j;
    char *dest = (char *) &f;
    char *src = (char *) &v[i];
  
#ifndef WORDS_BIGENDIAN 
  bswap(&v[i], 4, 1);
#endif  

    for (j=0; j<sizeof(float); j++)
        dest[j] = src[j];
    return f;
#else
    return v[i];
#endif
} 

local double show_d(double *v,int i)
{
#if 1
    double d;
    int j;
    char *dest = (char *) &d;
    char *src = (char *) &v[i];

#ifndef WORDS_BIGENDIAN 
  bswap(&v[i], 8, 1);
#endif  

    for (j=0; j<sizeof(double); j++)
        dest[j] = src[j];
    return d;
    
#else
    return v[i];
#endif    
}

local short show_i(short *v,int i)
{
#ifndef WORDS_BIGENDIAN 
  bswap(&v[i], 2, 1);
#endif  
    return v[i];
} 

local int show_j(int *v,int i)
{
#ifndef WORDS_BIGENDIAN 
  bswap(&v[i], 4, 1);
#endif  
    return v[i];
} 

local int8 show_k(int8 *v,int i)
{
#ifndef WORDS_BIGENDIAN 
  bswap(&v[i], 8, 1);
#endif  
    return v[i];
} 

