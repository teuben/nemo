/*
 * subfile used by NEMO's fitsgids.c
 *   kept in separate file in order to avoid collisions between
 *   NEMO & GIPSY's conventions
 *
 *	14-may-93	created		pjt
 *	20-jul-94	caching data - one line at a time was too slow 	pjt
 *			this improved loading 512*512 from 1'43" to 3" !!
 *	12-jan-95       using local myerror, which calls error
 *       1-feb-98       activating on myerror()
 *
 * mods from lookc.c:
 *      - local defs of gipsy datatypes because of cc/gcc difficulties
 *	- using NEMO's error() routine
 *	- local main() variables now global static's
 *	- split main() into look_init(), look_loop() and look_finis()
 *
 */

#include <stdio.h>

#define MAXDATA	16384		/* cache size */

/* local GIPSY definitions ...  */
typedef long fint;
typedef struct { char *a; fint  l; } fchar;

/* now local variables */
static   fchar fstr;
static   char  string[1];
static   fint idata[MAXDATA], ndata;
static   fint gid=0, giderr=0;
static   fint minc, maxc, nc, bc;
static   fint blo[2], bhi[2];
static   fint one;
static   fint i, j, k;
static   fint isize, jsize, n;
static   fint mrec, nrec, irec=0;
static   float bscale, bzero, dmin, dmax;
static	 int run=1;
static   int nrecords, records[1024];
static   char emsg[80];

int look_init( filename, nx, ny, datamin, datamax, first_id)
   char *filename;
   int nx, ny;
   float datamin, datamax;
   int first_id;
{
   run = first_id;
   if (run < 0) run=0;  	/* bypass all gdi calls - for debug */
   ndata = 0;			/* pointer into idata[] array */
   dmin = datamin;
   dmax = datamax;
   one = 1;
   blo[0] = 1;
   blo[1] = 1;
   bhi[0] = nx;
   bhi[1] = ny;
   fstr.a = string;
   fstr.l = 1;
   string[0] = ' ';

if(run)	
   gid = gdi_open_c( fstr );					/* open GIDS */
   if (gid < 0) {
        myerror(emsg, gid);
	error("GDI_OPEN error %d", gid );
   }
if(run)
   giderr = gdi_cinfo_c( &gid, &minc, &maxc, &nc, &bc ); /* color table info */
   if (giderr < 0) error("GDI_CINFO error %d", giderr );
   dprintf(0,"Color range: %d - %d Blank value: %d\n",
                minc, maxc, bc);
   if (nc<1) nc=1;                
   bscale = (datamax-datamin)/nc;
   bzero = datamin - bscale*minc;
   dprintf(0,"bscale=%g bzero=%g\n",bscale,bzero);
if(run)
   giderr = gdi_defimg_c( &gid, blo, bhi, &bscale, &bzero ); /* set size and */
   if (giderr < 0) error("GDI_DEFIMG error %d", giderr );	/* and scale */
if(run)
   giderr = gdi_rinfo_c( &gid, &nrec, &mrec );          /* get recording cap */
   if (giderr < 0) error("GDI_RINFO error %d", giderr );
   dprintf(0,"Max number of images that can be recorded: %d\n",mrec);
   dprintf(0,"%d images have been recorded so far\n",nrec);

   fstr.a = "FITSGIDS";
   fstr.l = strlen( fstr.a );
if(run)
   giderr = gdi_immid_c( &gid, fstr );		         /* identify program */
   if (giderr < 0 ) error("GDI_IMMID error %d", giderr );

   fstr.a = filename;
   fstr.l = strlen( fstr.a );
if(run)
   giderr = gdi_imsid_c( &gid, fstr );			   /* identify image */
   if (giderr < 0 ) error("GDI_IMSID error %d", giderr );
}

look_loop(j, buffer, nx)
   int j;
   float *buffer;      /* float is really FLOAT from nemo's fitsio.h !!! */
   int nx;
{
   int i, idat, nb=0;

   for ( i = 0; i < nx; i++) {
      if (buffer[i] < dmin) {		  /* clipping */
         idat = bc;
         nb++;
      } else if (buffer[i] > dmax) {      /* blanking */
         idat = maxc;				/* ???? or use bc ???? */
         nb++;
      } else
         idat = minc + (int) ((buffer[i]-dmin)*(maxc-minc)/(dmax-dmin));
      dprintf(1,"dat: %d %d %g %d\n",i,j,buffer[i],idat);
      idata[ndata++] = idat;
      if (ndata==MAXDATA) look_flush();             /* flush row */
   }

   return nb;
}

look_flush()
{
        if(run)
          giderr = gdi_imwrite_c( &gid, idata, &ndata, &one );
        if (giderr < 0) error("GDI_IMWRITE error %d", giderr );

        ndata=0;
}

look_record(i)
   int i;
{
   irec = i;
if(run)
   giderr = gdi_record_c( &gid, &irec);
   if (giderr>=0) records[nrecords++] = i;
   return giderr>=0;
}

look_finis()
{
if(run)
   giderr = gdi_sequence_c(&gid, records, &nrecords);
if(run)
   giderr = gdi_close_c( &gid );				/* close connection */
   if (giderr < 0) error("GDI_CLOSE error %d", giderr );
   return 0;
}
   
#if 1

myerror(char *fmt, int error_code)
{
   char msg[128];
   
   gdi_error_c(&error_code, msg);
   warning(msg);
/*   error(fmt,error_code); */
}

#endif
