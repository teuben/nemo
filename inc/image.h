/*
 * IMAGE.H: structured binary file definitions for standard 2- and 3D-image cubes
 * 
 *   Note:   data[ix][iy][iz] is the notation, but can be stored in C or Fortran fashion
 *         
 *
 *  30-Jun-87    V2.0 structures format for 2D maps, also allows series of maps  PJT
 *  23-dec-88	 V2.3 velocity added to header
 *  18-jan-89	 V3.0 added 3D
 *  30-jan-89    V3.1 Vel deleted from image struct
 *   1-feb-89    V4.0 Fortran definition to store matrices (it's simple to
 *		      to switch : use #define FORDEF or CDEF (default is FORDEF)
 *  16-may-92    -- made include file re-entrant
 *  12-aug-92   toyed with alternate storage method : pointers to pointers (ala NumRecC)
 *              works only for 2D images now
 *  30-jul-93    V4.3 Added Unit member - should be upward compatible
 *  21-feb-94    V4.4 ANSI
 *  13-apr-96    V5.0 support for non-linear axes
 *  18-may-99	 V6.0 placement for the CD_i_j matrix used in FITS
 *  21-feb-00    V6.1 added mapX_image() routines to return pointer arrays 
 *   9-sep-02    V6.2 added copy_image()
 *  13-nov-02    V6.3 added image masks as a boolean image
 *   7-may-04    V7.0 added the notion of a reference point and value
 *                    by (for now optionally) using
 *                    why was i so lazy and didn't do this in 1987.....
 *   6-jan-05         added prototypes for wcsio.c
 *  19-oct-11    V8.0 optional intelligent array from Karma, see USE_IARRAY
 *                    i->frame[i->x[ix] + i->y[iy]+ i->z[iz]]
 *                    which seems 15% slower..... (ccdstat EVLA benchmark)
 *  22-may-21         added Object
 *  13-dec-22         added various frequently used FITS header items for fitsccd-ccdfits conversions
 *  14-sep-22         Also allow more common names in FITS (CDELTi,CRVALi,CRPIXi)
 */
#ifndef _h_image
#define _h_image

#include <matdef.h>

typedef struct {          // image_axis
  int    nr;              // FITS 'naxisN' - length of axis
  string name;            // FITS 'ctype'
  string unit;            // FITS 'cunit'

  int  beamtype;
  real beamsize;

  int   type;		  /* axis type: 0 not defined, 1 linear, 2 array */
    
  real  rmin;             /* 'crval' at 0 ; linear axis */
  real  dr;               /* 'cdelt' */
  real  refpix;           /* normally 0, but like 'crpix' in FITS */
  real  refval;           /* normally rmin, but like 'crval' in FITS */

  real *val;              /* array of length 'nr' for coordinates if non-linear */
} image_axis, *image_axisptr;

typedef struct {        // image_mask
  bool    *frame;	/* pointer to a contiguous block of data */
  bool   **matrix;      /* 2D special case: pointers to pointers */
  bool  ***cube;        /* 3D special case: ptr to ptr to ptr's  */
  // @todo: needs USE_IARRAY
  int   nx;		/* dimensions in X, Y and Z */
  int   ny;
  int   nz;
} image_mask, *image_maskptr;


typedef struct {        // image
    real    *frame;     /* pointer to a contiguous block of data */
    real   **matrix;    /* 2D special case: pointers to pointers */
    real  ***cube;      /* 3D special case: ptr to ptr to ptr's  */
#ifdef USE_IARRAY
    int   *x;           /* iarray offset X array into *frame */
    int   *y;           /* iarray offset Y array into *frame */
    int   *z;           /* iarray offset Z array into *frame */
#endif

    int   axis;         /* new style axis (image_axis x,y,z) ??  */
                        // FITS can also use WCSAXES=3

    int   nx;		/* dimensions in X, Y and Z */
    int   ny;
    int   nz;
  
    real  xmin;         /* coordinates at xref etc. */
    real  ymin;	        /*   --- which is at center of a cell/voxel --- !!  */
    real  zmin;         /* this is close to the CRVAL keyword */ 
  
    real  dx;           /* grid spacing (same for all pixels) */
    real  dy;
    real  dz;
                        /* the 'ref' coordinates only used for axis>0 */
    real  xref;         /* fake corner (used to hardcode 0,0,0 in lower left) */
    real  yref;         /* for new style axis */
    real  zref;         /* this is close to the CRPIX keyword */
  
    real  restfreq;     // FITS 'RESTFRQ'
    real  vlsr;         // FITS
                        // FITS 'VELREF' e.g. 257   / 1 LSR, 2 HEL, 3 OBS, +256 Radio

    image_axis  ax;     /* new axis descriptors -- not used yet */
    image_axis  ay;
    image_axis  az;

    char proj[16];      /* standard FITS WCS projection types */
    real rotang;        /* FITS: CROTA  - not used here yet */

#if 0
    matrix cd;          /* WCS: note, this can only handle NDIM by NDIM */
#endif
    real   map_min;	/* data min and max in data */
    real   map_max;
    int    beamtype;	/* beams - not very well used yet */
    real   beamx;       /* smoothing beams */
    real   beamy;
    real   beamz;
    string namex;       /* name of axes (could be a NULL) */
    string namey;
    string namez;
    string unitx;       /* unit of the axis (FITS: CUNITx) */
    string unity;   
    string unitz;   
    string unit;        /* FITS 'BUNIT'   units (could be a NULL) */
                        /* FITS 'BTYPE' (e.g. 'Intensity') not used yet */
    string object;      /* object name */
    string telescope;   // FITS 'TELESCOP'
    string instrument;  // FITS 'INSTRUME'
    string observer;    // FITS
    real   equinox;     // FITS, e.g. 2000.0
    string radecsys;    // FITS (e.g. FK5)
    string specsys;     // FITS (e.g. LSRK)
    string timesys;     // FITS (e.g. UTC)
    string ssysobs;     // GBT uses 'TOPOCENT'
    real   time;	/* time tag */
    string storage;	/* array stored in Fortran or C definition */
  
    image_mask *mask;   /* optional image mask */
} image, *imageptr;

typedef struct {        // new_image
    void    *frame;  	/* pointer to a contiguous block of data */
    void   **matrix;    /* 2D special case: pointers to pointers */
    void  ***cube;      /* 3D special case: ptr to ptr to ptr's  */
    // @todo: needs USE_IARRAY
    int	type;		/* data-type of the data */

    image_axis  ax;      /* axis descriptors */
    image_axis  ay;
    image_axis  az;

    real  map_min;	/* data min and max in cube */
    real  map_max;

    string unit;        /* map units (could be a NULL) */
    real   time;	/* time tag */
    string storage;	/* array stored in Fortran or C definition */
} new_image, *new_imageptr;


typedef struct {        // region
  int mode;             /* default mode is rectangular region blc-trc */
  int blc[3];           /* bottom lower (boundingbox) corner  in ix,iy,iz */
  int trc[3];           /* top right (boundingbox) corner in ix,iy,iz */
} region,  *regionptr;

// future?
// #define Nx_new(iptr)    ((iptr)->ax.nr)


#define Frame(iptr)	((iptr)->frame)
#define Axis(iptr)      ((iptr)->axis)
#define Nx(iptr)	((iptr)->nx)
#define Ny(iptr)	((iptr)->ny)
#define Nz(iptr)	((iptr)->nz)
#define IDx(iptr)       ((iptr)->x)
#define IDy(iptr)       ((iptr)->y)
#define IDz(iptr)       ((iptr)->z)
#define Xmin(iptr) 	((iptr)->xmin)
#define Ymin(iptr) 	((iptr)->ymin)
#define Zmin(iptr) 	((iptr)->zmin)
#define CRVAL1(iptr) 	((iptr)->xmin)
#define CRVAL2(iptr) 	((iptr)->ymin)
#define CRVAL3(iptr) 	((iptr)->zmin)
#define Dx(iptr)	((iptr)->dx)
#define Dy(iptr)	((iptr)->dy)
#define Dz(iptr)	((iptr)->dz)
#define CDELT1(iptr) 	((iptr)->dx)
#define CDELT2(iptr) 	((iptr)->dy)
#define CDELT3(iptr) 	((iptr)->dz)
#define Xref(iptr)      ((iptr)->xref)
#define Yref(iptr)      ((iptr)->yref)
#define Zref(iptr)      ((iptr)->zref)
#define CRPIX1(iptr) 	((iptr)->xref)
#define CRPIX2(iptr) 	((iptr)->yref)
#define CRPIX3(iptr) 	((iptr)->zref)
#define MapMin(iptr)	((iptr)->map_min)
#define MapMax(iptr)	((iptr)->map_max)
#define BeamType(iptr)	((iptr)->beamtype)
#define Beamx(iptr)	((iptr)->beamx)
#define Beamy(iptr)	((iptr)->beamy)
#define Beamz(iptr)	((iptr)->beamz)
#define Namex(iptr)     ((iptr)->namex)
#define Namey(iptr)     ((iptr)->namey)
#define Namez(iptr)     ((iptr)->namez)
#define CUNIT1(iptr)    ((iptr)->namex)
#define CUNIT2(iptr)    ((iptr)->namey)
#define CUNIT3(iptr)    ((iptr)->namez)
#define Unitx(iptr)     ((iptr)->unitx)
#define Unity(iptr)     ((iptr)->unity)
#define Unitz(iptr)     ((iptr)->unitz)
#define Unit(iptr)      ((iptr)->unit)
#define Object(iptr)    ((iptr)->object)
#define Telescope(iptr) ((iptr)->telescope)
#define Restfreq(iptr)  ((iptr)->restfreq)
#define Vlsr(iptr)      ((iptr)->vlsr)
#define Time(iptr)	((iptr)->time)
#define Storage(iptr)   ((iptr)->storage)
#define Mask(iptr)      ((iptr)->mask)

#define BLC(rptr)       ((rptr)->blc)
#define TRC(rptr)       ((rptr)->trc)



#ifdef USE_IARRAY
#define MapValue(i,ix,iy)	(*( (i)->frame + (i)->x[ix] + (i)->y[iy]))
#define CubeValue(i,ix,iy,iz)   (*( (i)->frame + (i)->x[ix] + (i)->y[iy] + (i)->z[iz]))
#else


/* row major:   data[iy][ix]   data[r][c]   offset = row*NUMCOLS + col    */
#if defined(CDEF)
#define MapValue(i,ix,iy)	(*( (i)->frame + iy + Ny(i)*(ix)))
#define CubeValue(i,ix,iy,iz)   (*( (i)->frame + (size_t)iz + (size_t)Nz(i)*((size_t)iy + (size_t)Ny(i)*((size_t)ix))))
#endif

/* column major:    data(ix,iy)  data(r,c)   offset = col*NUMROWS + row    */
#if defined(FORDEF)
#define MapValue(iptr,ix,iy)	 (*( (iptr)->frame + ix + Nx(iptr)*(iy)) )
#define CubeValue(iptr,ix,iy,iz) (*( (iptr)->frame + ix + Nx(iptr)*(iy+Ny(iptr)*(iz))))
#endif

#endif

/* 
 *  BeamTypes -- not really used though
 */

#define NONE	       -1
#define ANYBEAM		0
#define	HANNING		1
#define GAUSS		2

/* 
 *  DataTypes -- for new 5.0 images
 */

#define BYTE_IMAGE	1
#define SHORT_IMAGE	2
#define INT_IMAGE	3
#define LONG_IMAGE	4
#define FLOAT_IMAGE	5
#define DOUBLE_IMAGE	6
#define REAL_IMAGE	7

/*
 * Item tags for Image components.
 */

#define ImageTag		"Image"

#define   ParametersTag		"Parameters"
#define     NxTag		"Nx"
#define     NyTag		"Ny"
#define     NzTag		"Nz"
#define     DxTag		"Dx"
#define	    DyTag		"Dy"
#define	    DzTag		"Dz"
#define	    XminTag		"Xmin"
#define	    YminTag		"Ymin"
#define	    ZminTag		"Zmin"
#define	    XrefTag		"Xrefpix"
#define	    YrefTag		"Yrefpix"
#define	    ZrefTag		"Zrefpix"
#define	    MapMinTag		"MapMin"
#define	    MapMaxTag		"MapMax"
#define	    BeamTypeTag		"BeamType"
#define	    BeamxTag		"Beamx"
#define     BeamyTag	  	"Beamy"
#define     BeamzTag	  	"Beamz"
#define     NamexTag		"Namex"
#define     NameyTag		"Namey"
#define     NamezTag		"Namez"
#define     UnitxTag		"Unitx"
#define     UnityTag		"Unity"
#define     UnitzTag		"Unitz"
#define     UnitTag             "Unit"
#define     ObjectTag           "Object"
#define     TelescopeTag        "Telescope"
#define     RestfreqTag         "Restfreq"
#define     VlsrTag             "VLSR"
#define	    TimeTag		"Time"		/* note: from snapshot.h  */
#define     StorageTag	        "Storage"
#define     AxisTag             "Axis"

#define     MapTag		"Map"
#define     MapValuesTag	"MapValues"

int minmax_image       (imageptr);
int write_image        (stream, imageptr);
int read_image         (stream, imageptr *);
int free_image         (imageptr);
int create_image       (imageptr *, int, int);
int create_image_mask  (imageptr, image_maskptr *);
int create_cube        (imageptr *, int, int, int);
int create_header      (imageptr);
int copy_image         (imageptr, imageptr *);
int copy_header        (imageptr, imageptr, int);

real  **map2_image(imageptr);
real ***map3_image(imageptr);

/* worldpos.c */
int worldpos(double xpix, double ypix, double xref, double yref, double xrefpix, double yrefpix, double xinc, double yinc, double rot, char *type, double *xpos, double *ypos);
int xypix(double xpos, double ypos, double xref, double yref, double xrefpix, double yrefpix, double xinc, double yinc, double rot, char *type, double *xpix, double *ypix);

/* wcsio.c */
void wcs_f2i(int ndim, double *crpix, double *crval, double *cdelt,            image *iptr);
void wcs_i2f(image *iptr, 	     int ndim, double *crpix, double *crval, double *cdelt);


/* get_nan.c */
void get_nanf(float *x);
void get_nand(double *x);

#endif
