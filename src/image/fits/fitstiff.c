/*
 *   FITSTIFF:   convert fits files  to tiff files, needs libtiff.a
 *
 *	 6-nov-91	V0    Crude version                     PJT
 *      15-nov-92       V0.1  renamed bmin/bmax -> datamin/max  PJT
 *                            and wrote smart getminmax(). Also
 *                            wrote getcolormap() to load lut's
 *	24-may-95	V0.1a fixed fitsio convention Y=0..NY-1	pjt
 *			 0.2  and added compress=t
 *                            and enforce FITS to be read contiguously
 *			      if not compressed
 *	 5-oct-95       0.3     ??
 *	18-jul-01	added usage line
 *      17-dec-03       0.4   optional 
 *
 * Note: the compres=f may now have a bug due to random access of TIFF ???
 */

#include <stdinc.h>
#include <getparam.h>
#include <filefn.h>     /* for pathopen, defext */
#include <image.h>
#include <fitsio.h>
#include <tiffio.h>	/* try -I$NEMOINC/tiff */

string defv[] = {
    "in=???\n		Input fits file",
    "out=???\n		Output tiff file",
    "planes=\n          Planes to select from fits cube [all]", 
    "datamin=\n         Min, if to override FITS 'DATAMIN'",
    "datamax=\n         Max, if to override FITS 'DATAMAX'",
    "lut=\n             Associated colormap; greyscale if none",
    "invert=f\n		Invert the colortable?",
    "compress=t\n	LZW Compression turned on?",
    "VERSION=0.3a\n	19-jul-01 PJT",
    NULL,
};

string usage="convert fits files to tiff files";

#define MAXPLANES 512
#define MAXCOLOR  256

local bool getminmax(FITS *, real* , int);
local bool getcolormap(string, short *, short *, short *, bool);

nemo_main()
{
    stream outstr;
    FITS *fitsfile, *rawopen();
    TIFF *tif;
    int ndim=3, naxis[3], nx, ny, nz, i, j, j1, j2, k, npl, p, planes[MAXPLANES];
    int nbval=0;
    real tmpa[3], crpix1, crpix2, crpix3, bval_out, bmin, bmax;
    FLOAT cdelt[3], *buffer, *bp, bval_in, tmpr;  /* fitsio- is in FLOAT !!! */
    imageptr iptr;
    string mode, tmps;
    bool   Qblank, Qcompress;
    unsigned char *outbuf, *cp;
    short r[MAXCOLOR], g[MAXCOLOR], b[MAXCOLOR];

    npl = nemoinpi(getparam("planes"),planes,MAXPLANES);
    Qcompress = getbparam("compress");


    /* Process FITS file header */

    fitsfile = fitopen(getparam("in"),"old",ndim,naxis);
    if (fitsfile==NULL) error("Could not open fitsfile %s",getparam("in"));
    if (hasvalue("datamin"))
        bmin = getdparam("datamin");
    else if (!getminmax(fitsfile,&bmin,1))
        error("FITS file does not contain DATAMIN: you need to supply one");

    if (hasvalue("datamax"))
        bmax = getdparam("datamax");
    else if (!getminmax(fitsfile,&bmax,0))
        error("FITS file does not contain DATAMAX: you need to supply one");

    dprintf(0,"Data will be scaled in the range %g - %g\n",bmin,bmax);

    nx = naxis[0];
    ny = naxis[1];
    nz = (npl>0) ? npl : naxis[2];
    if (nx*ny*nz==0) error("Bad fits image: nx*ny*nz=0");
    if (nx==1) warning("Fits image: nx=1 ???");
    if (ny==1) warning("Fits image: ny=1 ???");

    /* 
     * Open TIFF file for writing; note it will overwrite a
     * previously existing dataset,
     * and set some initial TIFF tags.
     */
 
    tif = TIFFOpen(getparam("out"),"w");
    if (tif==NULL) error("Could not open TIFF file %s",getparam("out"));
    TIFFSetField(tif,TIFFTAG_IMAGEWIDTH, (short) nx);
    TIFFSetField(tif,TIFFTAG_IMAGELENGTH, (short) ny);
    TIFFSetField(tif,TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif,TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif,TIFFTAG_PLANARCONFIG, 2);
    if (Qcompress)
        TIFFSetField(tif,TIFFTAG_COMPRESSION, 5);	/* lzw */

    if (hasvalue("lut")) {
        TIFFSetField(tif,TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_PALETTE);
        if (getcolormap(getparam("lut"),r,g,b,getbparam("invert"))) {
            TIFFSetField(tif,TIFFTAG_COLORMAP,r,g,b);
        }
    } else {
        TIFFSetField(tif,TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    }


    buffer = (FLOAT *) malloc(nx*sizeof(FLOAT));
    if (buffer==NULL) error("No memory for 1D buffer of size %d",nx);
    outbuf = (unsigned char *) malloc(nx);

    for (k=0; k<nz; k++) {                    /* loop over all/selected planes */
	if (k>0) {
            warning("Can write only one plane");
	    break;   /* TODO: dunno what to do with multiple planes yet */
        }
        p = (npl>0) ? planes[k]-1 : k;        /* select plane number */
        dprintf(1,"Reading plane %d\n",p+1);
        fitsetpl(fitsfile,1,&p);
        for (j=0; j<ny; j++) {                 /* loop over all rows */
            j1 = Qcompress ? ny-1-j : j;
            j2 = Qcompress ? j : ny-1-j;
            fitread(fitsfile,j1,buffer);     /* read line from fits file */
            for (i=0, bp=buffer, cp=outbuf; i<nx; i++, cp++, bp++) { 
                *cp = (char) (255*(*bp-bmin)/(bmax-bmin));
            }
            TIFFWriteScanline(tif,outbuf,j2,k);	/* write scaled to TIFF */
        }
    }

    fitclose(fitsfile);
    TIFFClose(tif);
}


#if 0
/* Capture TIFFError and TIFFWarning */

void TIFFError(mod, fmt, arg)
#endif


/*
 *  Try and get Min or Max from a FITS header only
 */

local bool getminmax(FITS *ff, real *val, int needmin)
{
    FLOAT tmp, bscale, bzero;
    real datamin, datamax;
    int bitpix;

    /*
     *  First figure out if DATAMIN/MAX were present in
     *  the header. If so, return these.
     */

    if (needmin) {
        if (fitexhd(ff,"DATAMIN")) {
            fitrdhdr(ff,"DATAMIN",&tmp,0.0);
            *val = tmp;
            return TRUE;
        }
    } else {
        if (fitexhd(ff,"DATAMAX")) {
            fitrdhdr(ff,"DATAMAX",&tmp,0.0);
            *val = tmp;
            return TRUE;
        }
    }

    /*
     * Else reconstruct them from assuming, for
     * positive bitpix (integer data) only, that the 
     * data scales over the full range of the unsigned 
     * integer
     */

    fitrdhdi(ff,"BITPIX",&bitpix,0);
    if (bitpix < 0) return FALSE;
    fitrdhdr(ff,"BSCALE",&bscale,1.0);
    fitrdhdr(ff,"BZERO",&bzero,0.0);
    if (bitpix == 8) {
        datamin = bzero - bscale * 128.0;
        datamax = bzero + bscale * 127.0;
    } else if (bitpix == 16) {
        datamin = bzero - bscale * 32768.0;
        datamax = bzero + bscale * 32767.0;
    } else if (bitpix == 32) {
        datamin = bzero - bscale * 2147483648.0;
        datamax = bzero + bscale * 2147483647.0;
    } else
        error("Illegal value of BITPIX = %d in fitsfile",bitpix);

    *val = (needmin ? datamin : datamax);
 
    return TRUE;
}


#define NOCOLOR(x)  ((x)<0||(x)>1)

local bool getcolormap(string name,short *r,short *g,short *b, bool invert)
{
    stream cstr;
    int i, j, ncolors=0;
    float red, green, blue;
    char line[256], *n, lutpath[256];
    short c;

    if (name==NULL || *name == NULL) return NULL;

    n = getenv("NEMODAT");
    if (n)
        sprintf(lutpath,".:%s/lut",n);
    else
        sprintf(lutpath,".");

    cstr = pathopen(lutpath,name, "r");
    if (cstr==NULL)
        cstr = pathopen(lutpath,defext(name,".lut"),"r");
    if (cstr==NULL)
        error("Cannot find LUT table %s in %s\n",name,lutpath);
   
    for (i=0; i<MAXCOLOR; i++)  
        r[i] = g[i] = b[i] = 0;
    while (fgets(line,256,cstr)) {
        if (line[0] == '#' || line[0] == ';' || line[0] == '\n') continue;
        if(ncolors>=MAXCOLOR) error("(%d/%d): Too many colors in %s",
              ncolors+1, MAXCOLOR, name);
        sscanf(line,"%f %f %f",&red, &green, &blue);
        if (NOCOLOR(red) || NOCOLOR(green) || NOCOLOR(blue)) {	
          warning("Skipping RGB=%f %f %f",red,green,blue);	/* ask KGB */
          continue;
        }
        dprintf(1,"LUT(%d): %f %f %f\n",ncolors,red,green,blue);
        r[ncolors] = (unsigned short) (65535*red);
        g[ncolors] = (unsigned short) (65535*green);
        b[ncolors] = (unsigned short) (65535*blue);
        ncolors++;
    }
    strclose(cstr);
    dprintf(0,"Colortable %s: %d entries\n", name, ncolors);
    if (invert) {
    	for (i=0, j=ncolors-1; i<ncolors/2; i++, j--) {
            c = r[i]; r[i] = r[j]; r[j] = c;
            c = g[i]; g[i] = g[j]; g[j] = c;
            c = b[i]; b[i] = b[j]; b[j] = c;
        }
    }

    return TRUE;
}
