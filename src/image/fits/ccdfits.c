/* 
 *      CCDFITS: convert CCD image to a fits file
 *
 *	Also contains a controversial 'long strings' implementation
 *
 *      29-apr-88  V1.0 created (again) PJT
 *       2-jun-88  V1.1 renamed 'wfits' to 'ccdfits'  
 *                      new filestruct: code is almost same
 *      23-dec-88  V1.2 added velocity information to header
 *      30-jan-89  V1.3 velocity information is now in Z
 *       6-mar-89  V1.3a 2D cube written when NZ=1
 *      28-jul-90  V2.0 Bob's new fitsio.c routines, finally
 *       1-oct-90  V2.1 default scale is now 1, not 1/60   PJT
 *                      crpix is now real !!
 *	11-oct-90  V2.2 blocking factor introduced 	   PJT
 *      14-feb-91  V2.3 write axes names for AIPS          PJT
 *	 7-mar-92  V2.4 gcc2.0 happy - new style	   PJT
 *	23-feb-93  V2.5  axisnames (CTYPE's) now ok	   PJT
 *      18-may-93  V2.6 history trial following sci.astro.fits proposal  PJT
 *      21-feb-94  V2.6a ansi
 *      22-nov-94      b ???
 *	20-mar-95      c long histories now written multiline		PJT
 *	 8-oct-95      d force CTYPEs for some WCS freaks		pjt
 *      28-mar-98      e write out Headline[] as COMMENT                PJT
 *      11-mar-99      e fixed bug in the above change                  pjt
 *	21-mar-99  V2.7  fixed scaling problems in bitpix=8		pjt
 *      18-may-99  V2.8  cdmatrix optional output                       pjt
 *      15-oct-99  V2.9  force a more silly RA---SIN/DEC--SIN axis      pjt
 *      23-mar-01  V3.0  allow fits reference image to inherit WCS from PJT
 *       8-apr-01      a fixed SINGLEPREC operation
 *      23-may-01  V3.1  added dummy= to be able to not write dummy axes  PJT
 *                        NOT WORKING YET
 *
 *  TODO:
 *      reference mapping has not been well tested, especially for 2D
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>
#include <history.h>
#include <fitsio.h>

string defv[] = {
        "in=???\n        Input image filename",
        "out=???\n       Output fits filename",
	"bitpix=-32\n	 FITS bitpix value {16,32,-32}",
        "scale=1,1,1\n   Extra scalefactor for cdelt&crval",
        "iscale=1,0\n    Scale and Offset Intensity conversion factor",
        "object=\n       Object name",
        "comment=\n      Extra comments",
	"cdmatrix=f\n    Use standard CD matrix instead of CDELT?",
	"blocking=1\n	 Blocking factor for output (blocksize/2880)",
	"refmap=\n       reference map to inherit WCS from",
	"refpix=\n       reference pixel, if different from default",
	"radecvel=f\n    Enforce reasonable RA/DEC/VEL axis descriptor",
	"dummy=t\n       Write dummy axes also ? ",
        "VERSION=3.1\n   23-may-01 PJT",
        NULL,
};

string usage = "convert image to a fits file";

stream  instr, outstr;                         /* file streams */

imageptr iptr=NULL;                     /* image, allocated dynamically */

double scale[3];        /* scale conversion for FITS (CDELT) */
double iscale[2];	/* intensity rescale */
char *object;           /* name of object in FITS header */
char *comment;          /* extra comments */
char *headline;         /* optional NEMO headline, added as COMMENT */
bool Qcdmatrix;         /* writing out new-style cdmatrix ? */
bool Qradecvel;         /* fake astronomy WCS header */
bool Qrefmap;
bool Qdummy;            /* write dummy axes ? */

int   nref = 0;
FLOAT ref_crval[3], ref_crpix[3], ref_cdelt[3];
char  ref_ctype[3][80];

void setparams(void);
void write_fits(string,imageptr);
void stuffit(FITS *, char *);
void set_refmap(string);

void nemo_main()
{
    setparams();                               /* set cmdln par's */
    instr = stropen (getparam("in"), "r");     /* open image file */
    read_image (instr,&iptr);                  /* read image */
    headline = ask_headline();                 /* possible headline */
    strclose(instr);                           /* close image file */
    write_fits(getparam("out"),iptr);          /* write fits file */
    free_image(iptr);
}

void setparams(void)
{
  int i,n;
  real tmpr[3];

  switch (nemoinpd(getparam("scale"),scale,3)) {
  case 1:
    scale[1] = scale[0];
    scale[2] = 1.0;
    break;
  case 2:
    scale[2] = 1.0;
    break;
  case 3:
    break;
  case 0:
    scale[0] = scale[1] = scale[2] = 1.0;
    break;
  default:
    error("parsing error scale=%s",getparam("scale"));
  }
  if (nemoinpd(getparam("iscale"),iscale,2) != 2)
    error("parsing error scale=%s - must be 2 numbers",getparam("scale"));
  object = getparam("object");
  comment = getparam("comment");
  Qcdmatrix = getbparam("cdmatrix");
  Qradecvel = getbparam("radecvel");
  Qrefmap = hasvalue("refmap");
  if (Qrefmap) {
    set_refmap(getparam("refmap"));
    if (hasvalue("refpix")) {
      nref =  nemoinpr(getparam("refpix"),tmpr,3);
      for (i=0; i<nref; i++)
	ref_crpix[i] = tmpr[i];
    }
  }
  Qdummy = getbparam("dummy");
}

static string ctypes[3] = { "CTYPE1", "CTYPE2", "CTYPE3" };
static string cdelts[3] = { "CDELT1", "CDELT2", "CDELT3" };
static string crvals[3] = { "CRVAL1", "CRVAL2", "CRVAL3" };
static string crpixs[3] = { "CRPIX1", "CRPIX2", "CRPIX3" };

void write_fits(string name,imageptr iptr)
{
    int  nx, ny, nz, i, j, k, bitpix;
    FLOAT tmpr,xmin,ymin,zmin,dx,dy,dz,mapmin,mapmax;   /* fitsio FLOAT !!! */
    char *cp;
    FITS *fitsfile;
    string *hitem;
    float *buffer, *bp;
    int ndim, keepaxis[3], naxis[3];   /* at most 3D cubes for now */
    double bscale, bzero;

    
    nx = naxis[0] = Nx(iptr);
    ny = naxis[1] = Ny(iptr);
    nz = naxis[2] = Nz(iptr);   if (nz <= 0) nz = 1;
    xmin = Xmin(iptr)*scale[0];
    ymin = Ymin(iptr)*scale[1];
    zmin = Zmin(iptr)*scale[2];
    dx = Dx(iptr)*scale[0];
    dy = Dy(iptr)*scale[1];
    dz = Dz(iptr)*scale[2];
    mapmin = MapMin(iptr);
    mapmax = MapMax(iptr);
    for (i=0; i<3; i++) {
      keepaxis[i] = ((naxis[i] > 1) ? 1 : 0);
      if (Qdummy) keepaxis[i] = 1;
    }

    dprintf(1,"NEMO Image file written to FITS disk file\n");
    dprintf(1,"%d %d %d   %f %f %f   %f %f %f   %f %f \n",
        nx,ny,nz,xmin,ymin,zmin,dx,dy,dz,mapmin,mapmax);

    ndim = (nz > 1 ? 3 : 2);    /* Make it 2D or real 3D map/cube */

    fit_setblocksize(2880*getiparam("blocking"));
    bitpix = getiparam("bitpix");
    fit_setbitpix(bitpix);
    if (bitpix == 16) {      /* scale from -2^(bitpix-1) .. 2^(bitpix-1)-1 */
        bscale = (mapmax - mapmin) / (2.0*32768.0 - 1.0);
        bzero = mapmax - bscale*32767.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 32) {
        bscale = (mapmax - mapmin) / (2.0*2147483648.0 - 1.0);
        bzero = mapmax - bscale*2147483647.0;
        fit_setscale(bscale,bzero);
    } else if (bitpix == 8) {
    	bscale = (mapmax - mapmin) / (2.0*128.0 - 1.0);
    	bzero = mapmin;
    	fit_setscale(bscale,bzero);
    }
    dprintf(1,"bscale,bzero=%g %g\n",bscale,bzero);

    fitsfile = fitopen(name,"new",ndim,naxis);
    if (fitsfile==NULL) error("Could not open fitsfile %s for writing\n",name);

    if (Qrefmap) {
      fitwrhdr(fitsfile,"CRPIX1",ref_crpix[0]);       
      fitwrhdr(fitsfile,"CRPIX2",ref_crpix[1]);       
      if (nz>1) fitwrhdr(fitsfile,"CRPIX3",ref_crpix[2]);
    } else {
      fitwrhdr(fitsfile,"CRPIX1",1.0);        /* CRPIX = 1 by Nemo definition */
      fitwrhdr(fitsfile,"CRPIX2",1.0);
      if (nz>1) fitwrhdr(fitsfile,"CRPIX3",1.0);
    }
    if (Qrefmap) {
      fitwrhdr(fitsfile,"CRVAL1",ref_crval[0]);
      fitwrhdr(fitsfile,"CRVAL2",ref_crval[1]);
      if (nz>1) fitwrhdr(fitsfile,"CRVAL3",ref_crval[2]);
    } else {
      fitwrhdr(fitsfile,"CRVAL1",xmin);
      fitwrhdr(fitsfile,"CRVAL2",ymin);
      if (nz>1) fitwrhdr(fitsfile,"CRVAL3",zmin);
    }

    if (Qcdmatrix) {
      fitwrhdr(fitsfile,"CD1_1",dx);    
      fitwrhdr(fitsfile,"CD2_2",dy);    
      if (nz>1) fitwrhdr(fitsfile,"CD3_3",dz);    
    } else {
      if (Qrefmap) {
	fitwrhdr(fitsfile,"CDELT1",ref_cdelt[0]*scale[0]);
	fitwrhdr(fitsfile,"CDELT2",ref_cdelt[1]*scale[1]);
	if (nz>1) fitwrhdr(fitsfile,"CDELT3",ref_cdelt[2]*scale[2]);
      } else {
	fitwrhdr(fitsfile,"CDELT1",dx);    
	fitwrhdr(fitsfile,"CDELT2",dy);    
	if (nz>1) fitwrhdr(fitsfile,"CDELT3",dz);
      }
    }
    if (Qradecvel) {
      dprintf(0,"[Axes names written as RA-SIN, DEC-SIN, VELO-LSR]\n");
      fitwrhda(fitsfile,"CTYPE1","RA---SIN");
      fitwrhda(fitsfile,"CTYPE2","DEC--SIN");
      if (nz>1) fitwrhda(fitsfile,"CTYPE3","VELO-LSR");
    } else {
      if (Qrefmap) {
        fitwrhda(fitsfile,"CTYPE1",ref_ctype[0]);
        fitwrhda(fitsfile,"CTYPE2",ref_ctype[1]);
        if (nz>1) fitwrhda(fitsfile,"CTYPE3",ref_ctype[2]);
      } else {
        if (Namex(iptr))
            fitwrhda(fitsfile,"CTYPE1",Namex(iptr));
        else
    	    fitwrhda(fitsfile,"CTYPE1","X");
        if (Namey(iptr))
            fitwrhda(fitsfile,"CTYPE2",Namey(iptr));
        else
    	    fitwrhda(fitsfile,"CTYPE2","Y");
        if (nz>1) {
            if (Namez(iptr))
                fitwrhda(fitsfile,"CTYPE3",Namez(iptr));
            else
      	        fitwrhda(fitsfile,"CTYPE3","Z");
      	}
      }
    }

    fitwrhdr(fitsfile,"DATAMIN",mapmin);
    fitwrhdr(fitsfile,"DATAMAX",mapmax);
    fitwrhda(fitsfile,"ORIGIN","NEMO");

    cp = getenv("USER");                                /* AUTHOR */
    if (cp && *cp)
        fitwrhda(fitsfile,"AUTHOR",cp);
    else
        fitwrhda(fitsfile,"AUTHOR","NEMO");

    if (*object)                                        /* OBJECT */
        fitwrhda(fitsfile,"OBJECT",object);

    if (*comment)                                       /* COMMENT */
        fitwra(fitsfile,"COMMENT",comment);
#if 1
    if (headline && *headline)
        fitwra(fitsfile,"COMMENT",headline);
#endif

    hitem = ask_history();                              /* HISTORY */
    fitwra(fitsfile,"HISTORY","NEMO: History in reversed order");
    for (i=0, cp=hitem[0]; cp != NULL; i++) {
    	stuffit(fitsfile,cp);
        cp = hitem[i+1];
    }

    buffer = (float *) allocate(nx*sizeof(float));  /* MEMLEAK */

    for (k=0; k<nz; k++) {          /* loop over all planes */
        fitsetpl(fitsfile,1,&k);
        for (j=0; j<ny; j++) {      /* loop over all rows */
            for (i=0, bp=buffer; i<nx; i++, bp++)
                *bp =  iscale[0] * CubeValue(iptr,i,j,k) + iscale[1];
            fitwrite(fitsfile,j,buffer);
        }
    }
    free(buffer);
    fitclose(fitsfile);
}


/*
 * stuff a character string accross the 80-line boundary
 * NOTE: the exact implementation of this routines is still controversial

 From: pence@tetra.gsfc.nasa.gov Sat May 15, 1993 10:37 "Continuation Keywords"

 Option 1.   Blank Continuation Keyword Convention

 A continued string value is recognized by a backslash (\)
 as the last character in the first keyword string, followed
 by a keyword with a blank name followed by the continuation
 of the keyword value enclosed in quotes.  Example:

 KEYNAME = 'This is a very long keyword string value which \'
         'continues over several lines\'
                 ' of the FITS header.'

 *
 */
                 
void stuffit(FITS *fitsfile, char *cp)
{
    char line[81], *hp;
    int i;
    int maxlen = 69;	/* could be 71 if you want to fill the whole card */
			/* but for the sake of CR-LDF patchers we do 69 */

    hp = cp;
    strncpy(line,hp,maxlen);
    line[maxlen] = 0;
    fitwra(fitsfile,"HISTORY",line);
    while ((int)strlen(hp) > maxlen) {
        hp += maxlen;
        strncpy(line,hp,maxlen);
        line[maxlen] = 0;
        fitwra(fitsfile," ",line);
    }
}

void set_refmap(string name)
{
  FITS *fitsfile;
  FLOAT tmpr, one = 1.0;
  int ndim = 3;
  int naxis[3], tmpi;


  fitsfile = fitopen(name,"old",ndim,naxis);
  dprintf(0,"[Reading reference map %s [%d,%d,%d]\n",
	  name,naxis[0],naxis[1],naxis[2]);

  fitrdhdr(fitsfile,"CDELT1",&ref_cdelt[0], one);
  fitrdhdr(fitsfile,"CDELT2",&ref_cdelt[1], one);
  fitrdhdr(fitsfile,"CDELT3",&ref_cdelt[2], one);

  fitrdhdr(fitsfile,"CRPIX1",&ref_crpix[0], one);
  fitrdhdr(fitsfile,"CRPIX2",&ref_crpix[1], one);
  fitrdhdr(fitsfile,"CRPIX3",&ref_crpix[2], one);

  fitrdhdr(fitsfile,"CRVAL1",&ref_crval[0], one);
  fitrdhdr(fitsfile,"CRVAL2",&ref_crval[1], one);
  fitrdhdr(fitsfile,"CRVAL3",&ref_crval[2], one);

  fitrdhda(fitsfile,"CTYPE1",ref_ctype[0],"");
  fitrdhda(fitsfile,"CTYPE2",ref_ctype[1],"");
  fitrdhda(fitsfile,"CTYPE3",ref_ctype[2],"");

  fitclose(fitsfile);


}
