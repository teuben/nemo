#define MAXIMAGES 5
#define MAXAXES 4

/* values for use in axisname.type */
#define UNDEF 0
#define SPATIAL 1
#define VELOCITY 2
#define FREQUENCY 3
#define WAVELENGTH 4
#define TIME 5
#define LATTYPE 0x40
#define OFFSETTYPE 0x80

/* Axis units */
# define DEGREES 0x0			/* degrees */
# define RADIANS 0x1			/* radians */
# define HOURS 0x2			/* hours */
# define ARCMINUTES 0x3			/* arc minutes */
# define NONE 4
# define SECONDS 5
# define MEGAHERTZ 6
# define MICROMETERS 7
# define KMPERSEC 8

/* logical values in a FITS header */
#define SIMPLE 0
#define simple l[SIMPLE]
/* END doesn't really have a value, so set aside no storage for it */
#define END 1
#define NUMFITSLOG 1

/* integer values in a FITS header */
#define BITPIX 0
#define bitpix i[BITPIX]
#define NAXIS 1
#define naxis i[NAXIS]
#define NAXIS1 2
#define naxis1 i[NAXIS1]
#define NAXIS2 3
#define naxis2 i[NAXIS2]
#define NAXIS3 4
#define naxis3 i[NAXIS3]
#define NAXIS4 5
#define naxis4 i[NAXIS3]
#define BLANK 6
#define blank i[BLANK]
#define NUMFITSINT 7

/* string values in a FITS header */
#define OBJECT 0
#define object s[OBJECT]
#define TELESCOP 1
#define telescop s[TELESCOP]
#define DATE_MAP 2
#define date_map s[DATE_MAP]
#define BUNIT 3
#define bunit s[BUNIT]
#define CTYPE1 4
#define ctype1 s[CTYPE1]
#define CTYPE2 5
#define ctype2 s[CTYPE2]
#define CTYPE3 6
#define ctype3 s[CTYPE3]
#define CTYPE4 7
#define ctype4 s[CTYPE3]
#define NUMFITSSTR 8

/* Long strings from HISTORY COMB lines in a FITS header */
#define TITLE 0
#define title h[TITLE]
#define PARAM 1
#define param h[PARAM]
#define NUMFITSHIST 2

/* double values in a FITS header */
#define CRVAL1 0
#define crval1 d[CRVAL1]
#define CDELT1 1
#define cdelt1 d[CDELT1]
#define CRPIX1 2
#define crpix1 d[CRPIX1]
#define CROTA1 3
#define crota1 d[CROTA1]
#define CRVAL2 4
#define crval2 d[CRVAL2]
#define CDELT2 5
#define cdelt2 d[CDELT2]
#define CRPIX2 6
#define crpix2 d[CRPIX2]
#define CROTA2 7
#define crota2 d[CROTA2]
#define CRVAL3 8
#define crval3 d[CRVAL3]
#define CDELT3 9
#define cdelt3 d[CDELT3]
#define CRPIX3 10
#define crpix3 d[CRPIX3]
#define CROTA3 11
#define crota3 d[CROTA3]
#define CRVAL4 12
#define crval4 d[CRVAL3]
#define CDELT4 13
#define cdelt4 d[CDELT3]
#define CRPIX4 14
#define crpix4 d[CRPIX3]
#define CROTA4 15
#define crota4 d[CROTA3]
#define FEPOCH 16
#define fepoch d[FEPOCH]
#define BSCALE 17
#define bscale d[BSCALE]
#define BZERO 18
#define fbzero d[BZERO]
#define DATAMAX 19
#define datamax d[DATAMAX]
#define DATAMIN 20
#define datamin d[DATAMIN]
#define NUMFITSDOUBLE 21

/* file types */
#define TMP 0
#define EXTERNAL 1
#define FITS 4

struct of_imagehdr {
	char fname[64];			/* name of image file */
	int fd;				/* file descriptor */
	int ftype;			/* file type (tmp, external, ...) */
	int dataStart;			/* byte offset of first pixel in file */
	int buflen;			/* data buffer len in lines */
	int bufline;			/* offset (in lines) of buffer from
					 * bottom line in image */
	float *buf;			/* points to data buffer */
	char s[NUMFITSSTR][12];		/* strings from the FITS header */
	char h[NUMFITSHIST][60];	/* long strings from HISTORY COMB */
	char l[NUMFITSLOG];		/* logical values from the header */
	int i[NUMFITSINT];		/* integers from FITS header */
	double d[NUMFITSDOUBLE];	/* floats and doubles from header */
};
extern struct of_imagehdr *imageHdr[MAXIMAGES];

/*
 * Structure to keep the names of keywords we use and their position in the
 * imagehdr structure.
 */
#define LOGICAL 0
#define INT 1 
#define STR 2
#define DOUBLE 3
#define HISTORY 4

struct of_fitskey {
	char keyword[9];
	char type;
	short int index;
	short int sequence;
	short int special;
};

struct AXISNAMES {
	float factor;		/* multiply header values by factor */
	char unit;		/* e.g. ARCMINUTES preferred for display */
	char type;		/* whether spatial, longitude type, etc. */
	char fname[6];		/* Name to recognize in CTYPEn */
	char len;		/* number of chars to check in fname */
	char lname[13];		/* name to use in axis labels */
};


/* fits.c*/

extern void OpenFits(/*imnumber,fname*/);
extern void NewImage(/*imnumber*/);
extern void ReleaseImageBuff(/*imnumber*/);
extern void CloseImage(/*imnumber*/);
extern float * GetImageLine(/*imnumber,y,nlines,ydir*/);
extern void ImageXGrid(/*imnumber,lowx,highx,deltax,xunit*/);
extern void ImageYGrid(/*imnumber,lowy,highy,deltay,yunit*/);
extern void ImNumErr(/*imnumber*/);
extern struct AXISNAMES FitsAxisUnits(/*ctype*/);
extern char *ChkAlloc(/*size,name*/);
extern error(/*__builtin_va_alist*/);
extern round(/*d*/);
extern cifdif(/*a,b,fuzz*/);
