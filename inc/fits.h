/*
 *  definitions and declarations for nemo's FITS.C
 *  22-feb-94  ansi- and C++ safe
 *  10-aug-09  size_t instead of int for 2GB+ files
 */

#ifndef _fits_h_
#define _fits_h_

#define FTSLINSIZ       80      /* line size in fits header */
#define FTSBLKSIZ       2880    /* block size (excluding blocking factor) */
#define FTSLPB  (FTSBLKSIZ/FTSLINSIZ)   /* lines per block */
#define MAXNAXIS        999     /* largest NAXIS possible */


enum fits_token { I_SIMPLE, I_GROUPS, I_EXTEND, I_XTENSION, I_BITPIX,
                  I_BLOCKED,
                  I_NAXIS, I_BZERO, I_BSCALE, I_BUNIT,
                  I_MAXIS,
                  I_EPOCH, I_EQUINOX, I_DATE, I_DATE_OBS,
                  I_RESTFREQ,
                  I_EXTNAME, I_EXTVER, I_EXTLEVEL,
                  I_CRVAL, I_CRPIX, I_CDELT, I_CTYPE, I_CUNIT, I_CROTA,
                  I_DATAMIN, I_DATAMAX, 
                  I_BLANK, I_OBSRA, I_OBSDEC,
                  I_PCOUNT, I_GCOUNT, 
                  I_PZERO, I_PSCAL, I_PTYPE, I_PUNIT,
                  I_TFIELDS, I_TTYPE, I_TFORM, I_TUNIT, I_TNULL, I_TDIM, 
                  I_TBCOL, I_TSCAL, I_TZERO, I_TDISP, 
                  I_COMMENT, I_HISTORY, I_COMMAND,
                  I_INSTRUME, I_TELESCOP, I_ORIGIN, I_OBJECT, I_OBSERVER,
                  I_AUTHOR, I_REFERENC,
                  I_END, I_DONE } ;

static struct arglist {
        string name;        /* name of keyword (max 8 chars) */
        bool  needsvalue;   /* keyword, i.e. needs an equals with a value? */
        bool  isarray;      /* is an array? */
	enum fits_token token;   /* token of this keyword for rapid selection */
} fitsargs[] =  {
        {"SIMPLE",  TRUE,  FALSE,  I_SIMPLE},      /* t/f */
        {"GROUPS",  TRUE,  FALSE,  I_GROUPS},      /* t/f */
        {"EXTEND",  TRUE,  FALSE,  I_EXTEND},      /* t/f */
        {"XTENSION",TRUE,  FALSE,  I_XTENSION},    /* name of extension type */
        {"BITPIX",  TRUE,  FALSE,  I_BITPIX},      /* 8, 16, +/- 32, +/- 64 */
        {"BLOCKED", TRUE,  FALSE,  I_BLOCKED},     /* blocked mode? */
        {"NAXIS",   TRUE,  TRUE,   I_NAXIS},       /* number of axes */
        {"BZERO",   TRUE,  FALSE,  I_BZERO},       /* scaling of data matrix */
        {"BSCALE",  TRUE,  FALSE,  I_BSCALE},      /* scaling of data matrix */
        {"BUNIT",   TRUE,  FALSE,  I_BUNIT},       /* units in matrix */
        {"DATAMIN", TRUE,  FALSE,  I_DATAMIN},     /* minimum in map */
        {"DATAMAX", TRUE,  FALSE,  I_DATAMAX},     /* maximum in map */
        {"BLANK",   TRUE,  FALSE,  I_BLANK},

	{"MAXIS",   TRUE,  TRUE,   I_MAXIS},       /* number of multidim axes */

        {"EPOCH",   TRUE,  FALSE,  I_EPOCH},       /* epoch (deprecated) */
        {"EQUINOX", TRUE,  FALSE,  I_EQUINOX},     /* equinox */
        {"DATE",    TRUE,  FALSE,  I_DATE},        /* creation of FITS map */
        {"DATE-OBS",TRUE,  FALSE,  I_DATE_OBS},    /* observation date */

        {"RESTFREQ",TRUE,  FALSE,  I_RESTFREQ},    /* some frequency */
        
        {"EXTNAME", TRUE,  FALSE,  I_EXTNAME},
        {"EXTVER",  TRUE,  FALSE,  I_EXTVER},
        {"EXTLEVEL",TRUE,  FALSE,  I_EXTLEVEL},

        {"CRVAL",   TRUE,  TRUE,   I_CRVAL},       /* ref pixel value */
        {"CRPIX",   TRUE,  TRUE,   I_CRPIX},       /* ref pixel location */
        {"CDELT",   TRUE,  TRUE,   I_CDELT},       /* pixel increment */
        {"CTYPE",   TRUE,  TRUE,   I_CTYPE},       /* and ID name of axes */
        {"CUNIT",   TRUE,  TRUE,   I_CUNIT},       /* units of the axis */
        {"CROTA",   TRUE,  TRUE,   I_CROTA},       /* rotation angle */


        {"OBSRA",   TRUE,  FALSE,  I_OBSRA},
        {"OBSDEC",  TRUE,  FALSE,  I_OBSDEC},

        {"PCOUNT",  TRUE,  FALSE,  I_PCOUNT},      /* number of parameters */
        {"GCOUNT",  TRUE,  FALSE,  I_GCOUNT},      /* number of random groups */
        {"PZERO",   TRUE,  TRUE,   I_PZERO},       /* scaling of parameters */
        {"PSCAL",   TRUE,  TRUE,   I_PSCAL},       /* scaling of parameters */
        {"PTYPE",   TRUE,  TRUE,   I_PTYPE},       /* name of parameter */
        {"PUNIT",   TRUE,  TRUE,   I_PUNIT},       /* unit of parameter */

        {"TFIELDS", TRUE,  FALSE,  I_TFIELDS},     /* number of columns in table */
        {"TTYPE",   TRUE,  TRUE,   I_TTYPE},       /* name of column */
        {"TFORM",   TRUE,  TRUE,   I_TFORM},       /* format of column */
        {"TUNIT",   TRUE,  TRUE,   I_TUNIT},       /* units of column */
        {"TNULL",   TRUE,  TRUE,   I_TNULL},       /* what is a blank */
        {"TBCOL",   TRUE,  TRUE,   I_TBCOL},       /* starting column */
        {"TSCAL",   TRUE,  TRUE,   I_TSCAL},       /* scaling */
        {"TZERO",   TRUE,  TRUE,   I_TZERO},       /* offset */
	{"TDIM",    TRUE,  TRUE,   I_TDIM},        /* multidimensional array */
	{"TDISP",   TRUE,  TRUE,   I_TDISP},       /* suggested disp format */


        {"COMMENT", FALSE, FALSE,  I_COMMENT},
        {"HISTORY", FALSE, FALSE,  I_HISTORY},
        {"COMMAND", FALSE, FALSE,  I_COMMAND},

        {"INSTRUME",TRUE,  FALSE,  I_INSTRUME}, 
        {"TELESCOP",TRUE,  FALSE,  I_TELESCOP}, 
        {"ORIGIN",  TRUE,  FALSE,  I_ORIGIN},
        {"OBJECT",  TRUE,  FALSE,  I_OBJECT},
        {"OBSERVER",TRUE,  FALSE,  I_OBSERVER},
        {"AUTHOR",  TRUE,  FALSE,  I_AUTHOR},
        {"REFERENC",TRUE,  FALSE,  I_REFERENC},
                                    /* REQUIRED ENDING OF PARAMETERS */
        {"END",     FALSE, FALSE,  I_END},         /* end of header */

        {NULL,      FALSE, FALSE,  I_DONE}	   /* End Of Table signature */
};


/*   datatypes  */

#define FTS_NULL   0    /* 0: not defined */
#define FTS_BYTE   1    /* 8 */
#define FTS_SHORT  2    /* 16 */
#define FTS_LONG   3    /* 32 */
#define FTS_REAL   4    /* -32 */
#define FTS_DOUBLE 5    /* -64 */


/*      The fits structure to hold primary header */

typedef struct fits_header {
    /* Optionally we could store the buffer in the first structure too... */
    char *buffer;       /* pointer to exact header copy */
    size_t buflen;      /* current length of buffer */
    size_t hlen;        /* length of header (bytes, through the END keyword) */
    size_t dlen;        /* length of data (bytes) */
    size_t nwritten;    /* accumulated in fts_wdata() only */
    size_t nread;	/* accumulated in various routines ... */
    int flip;           /* flips bytes when read in fts_rdata; see DECORDER */
    int ascii;		/* allow unfilled ascii headers to be read ? */
    int simple;         /*  0=FALSE 1=TRUE  */
    int groups;         /*  0=FALSE 1=TRUE  */
    int extend;         /*  0=FALSE 1=TRUE  */
    int blocked;        /*  0=FALSE 1=TRUE  */
    int bitpix;         /*  8, 16, +/-32, +/-64 */
    int naxis;          /* number of axes */
    int *naxisn;        /* length of axes (array or NULL) */
    int maxis;
    int *maxisn;
    int blank;           /* only used for positive bitpix */
    double *crvaln;      /* coordinate at reference pixel (array or NULL) */
    double *crpixn;      /* reference pixel (array or NULL) */
    double *cdeltn;      /* pixel increment (array or NULL) */ 
    double *crotan;      /* rotation angle of axis */
    double datamin;      /* maximum in map */
    double datamax;      /* minimum in map */
    double obsra;
    double obsdec;
    double restfreq;
    char *instrume;
    char *telescop;
    char *origin;
    char *epoch;
    char *equinox;
    char *date;
    char *date_obs;
    char *object;
    char *bunit;        /* units of (primary?) matrix value */
    char *observer;
    char *author;
    char *referenc;
    char **ctypen;      /* coordinate type (array or NULL) */
    char **cunitn;      /* units (array or NULL) */
    char **history;     /* array of history items */
    char **comment;     /* array of comment items */
    char **command;     /* array of command items */
    int pcount;         /* parameter count in group extensions */
    int gcount;         /* group count in group extensions */
    double bzero;        /* scale offset factor for matrix values */
    double bscale;       /* scale factor: true value=BZERO+BSCALE*value */
    double *pzeron;      /* scale offset factor for parameters */
    double *pscaln;      /* scale factor */
    char **ptypen;      /* name of parameter */
    char **punitn;
    int tfields;
    int *tbcoln;
    int *tbitems;
    double *tscaln;
    double *tzeron;
    char *xtension;
    char *extname;
    int extver;
    int extlevel;
    char **tnulln;
    char **ttypen;
    char **tformn;
    char **tunitn;
    char **ttypen_c;  /* additional comments for ttypen'n */
    char **tdimn;
    char **tdispn;
} fits_header;


size_t fts_rhead    (fits_header *, stream);
char *fts_shead    (fits_header *, string);
int fts_lhead    (fits_header *);
int fts_chead    (fits_header *, stream);
int fts_thead    (fits_header *);
int fts_dhead    (fits_header *, string *);
int fts_khead    (fits_header *, string *);
int fts_ihead    (fits_header *, string *);
int fts_fhead    (fits_header *, string *);
int fts_phead    (fits_header *, string *);
int fts_whead    (fits_header *, stream);
size_t fts_xhead    (fits_header *, stream, int, int, int *, int);

int fts_ptable   (fits_header *, stream, string *, string, int *, int);
int fts_pgroup   (fits_header *, stream, string *, string, int *, int);

int fts_buf      (int);

int fts_sdata    (fits_header *, stream);
int fts_cdata    (fits_header *, stream, stream, bool, bool);
int fts_dsize    (fits_header *);
int fts_tsize    (fits_header *);
int fts_rdata    (fits_header *, stream, int, char *);
int fts_wdata    (fits_header *, stream, int, char *);
int fts_rrow     (fits_header *, stream, int, char *);
int fts_wrow     (fits_header *, stream, int, char *);
int fts_zero     (fits_header *);

int fts_chead816 (fits_header *, stream);
int fts_cdata816 (fits_header *, stream, stream);

int fts_wvar     (stream, string, string);
int fts_wvarc    (stream, string, string, string);
int fts_wvarc_a  (stream, string, int, string *, string);
int fts_wvarb    (stream, string, int, string);
int fts_wvarb_a  (stream, string, int, int *, string);
int fts_wvari    (stream, string, int, string);
int fts_wvari_a  (stream, string, int, int *, string);
int fts_wvarf    (stream, string, float, string);
int fts_wvarf_a  (stream, string, int, float *, string);
int fts_wvard    (stream, string, double, string);
int fts_wvard_a  (stream, string, int, double *, string);

int fts_setiblk  (int);
int fts_setoblk  (int);




#endif /* __FITS_H___ */


#if 0
static struct arglist {
        char  *name;        /* name of keyword (max 8 chars) */
        bool  needsvalue;   /* keyword, i.e. needs an equals with a value? */
        bool  isarray;      /* is an array? */
        enum {    I_SIMPLE, I_GROUPS, I_EXTEND, I_XTENSION, I_BITPIX,
                  I_BLOCKED,
                  I_NAXIS, I_BZERO, I_BSCALE, I_BUNIT,
                  I_MAXIS,
                  I_EPOCH, I_EQUINOX, I_DATE, I_DATE_OBS,
                  I_RESTFREQ,
                  I_EXTNAME, I_EXTVER, I_EXTLEVEL,
                  I_CRVAL, I_CRPIX, I_CDELT, I_CTYPE, I_CUNIT, I_CROTA,
                  I_DATAMIN, I_DATAMAX, 
                  I_BLANK, I_OBSRA, I_OBSDEC,
                  I_PCOUNT, I_GCOUNT, 
                  I_PZERO, I_PSCAL, I_PTYPE, I_PUNIT,
                  I_TFIELDS, I_TTYPE, I_TFORM, I_TUNIT, I_TNULL, I_TDIM, 
                  I_TBCOL, I_TSCAL, I_TZERO, I_TDISP, 
                  I_COMMENT, I_HISTORY, I_COMMAND,
                  I_INSTRUME, I_TELESCOP, I_ORIGIN, I_OBJECT, I_OBSERVER,
                  I_AUTHOR, I_REFERENC,
                  I_END, I_DONE } 
              token;
}              
#endif

