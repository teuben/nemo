/* 
 * CCDDIFFRACT: diffract a 2D image, if 3D each slice done indepedantly
 *
 *	 9-may-10  V0.1 adapted from ccdsmooth     PJT
 */

#include <stdinc.h>
#include <string.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <image.h>

string defv[] = {
	"in=???\n               Input filename",
	"out=???\n              Output filename",
	"gauss=0.1\n            FWHM gaussian beamwidth for the core image",		
	"fraction=0.1\n         Spike fraction (rest in gaussian beam)",
	"cutoff=0\n             Cutoff pixel value below which no diffraction spikes",
	"spike=0.5\n            Spike sinc width",
	"dir=xy\n               Smoothing direction(s)",
	"noise=0\n              Add poisson noise",
	"bad=\n			Optional ignoring this bad value",
	"VERSION=0.1\n          9-may-10 PJT",
	NULL,
};

string usage = "add diffraction spikes to an image";

#ifndef HUGE
# define HUGE 1.0e20
#endif

string	infile, outfile;			/* file names */
stream  instr, outstr;				/* file streams */

#define MSIZE  8196		      /* maximum # pixels along one dimension */
#define MSMOOTH 501 		    /* maximum full beam-size (has to be odd) */
	              /* because of symmetry, you could try and be smart here */

imageptr iptr=NULL;			/* will be allocated dynamically */
imageptr iptr1=NULL;
int    nx,ny,nz,nsize;			/* actual size of map */
real   xmin,ymin,zmin,dx,dy,dz;
real   size;				/* size of frame (square) */
real   cell;				/* cell or pixel size (square) */

real   smooth[MSMOOTH];			/* full 1D beam */
int    lsmooth;				/* actual smoothing length */
int    nsmooth;				/* number of smoothings */
real   gauss_fwhm;			/* beam size in case gauss */
real   gauss_cut;                       /* cutoff of gaussian */
string dir;                             /* direction to smooth 'xyz' */
real   fraction;                        /* fraction in spike */
real   spike;                           /* spike width */
real   cutoff;                          /* no spikes below this */
real   noise;                           /* add noise */

bool   Qbad;                            /* ignore smoothing for */
real   bad;                             /* this value */

void setparams(), report_minmax(), smooth_it(), make_gauss_beam();
int convolve_cube (), convolve_x(), convolve_y(), convolve_z();
int spike_x(), spike_y();
real sinc2();


void nemo_main ()
{
    setparams();			/* set par's in globals */

    instr = stropen (infile, "r");
    read_image (instr,&iptr);
				/* set some global paramters */
    nx = Nx(iptr);	
    ny = Ny(iptr);
    nz = Nz(iptr);
    xmin = Xmin(iptr);
    ymin = Ymin(iptr);
    zmin = Zmin(iptr);
    dx = Dx(iptr);
    dy = Dy(iptr);
    dz = Dz(iptr);
    create_cube(&iptr1, nx, ny, nz);

    report_minmax("old");

    if(hasvalue("gauss"))
        make_gauss_beam(getparam("dir"));

    outstr = stropen (outfile,"w");
    smooth_it();
    write_image (outstr,iptr);

    strclose(instr);
    strclose(outstr);
}

void setparams()
{
    infile = getparam ("in");
    outfile = getparam ("out");

    if(hasvalue("gauss")) {         /* gaussian beam */
        if(nemoinpr(getparam("gauss"),&gauss_fwhm,1) != 1)
            error ("gauss=%s : beam must have 1 dimension",
                    getparam("gauss"));
        gauss_cut = 0.0;
    } else
      error("need gauss=");
    spike = getdparam("spike");
    fraction = getdparam("fraction");
    cutoff = getdparam("cutoff");
    noise = getdparam("noise");

    nsmooth = 1;
    dir = getparam("dir");
    Qbad = hasvalue("bad");
    if (Qbad) bad = getdparam("bad");
}


void report_minmax(string t)
{
    real m_min, m_max, brightness, total;
    int    i, ix, iy, iz, kounter, idir;
    char   *cp;

    m_max = -HUGE;                      /* determine new min/max */
    m_min =  HUGE;
    for (ix=0; ix<Nx(iptr); ix++)   	
    for (iy=0; iy<Ny(iptr); iy++)
    for (iz=0; iz<Nz(iptr); iz++) {
          brightness = CubeValue(iptr,ix,iy,iz);
	  total += brightness;
          m_max = MAX(m_max, brightness);
          m_min = MIN(m_min, brightness);
    }
    dprintf (0,"%s min and max in map are: %f %f\n",t,m_min,m_max);
    dprintf (0,"%s total values: %f\n",t,total);
    dprintf (0,"%s total surden: %f\n",t,total*Dx(iptr)*Dy(iptr));
}

int spike_x(d,ix,iy,iz)
     real d;
     int ix,iy,iz;

{
  int i, idx;
  real x;
  static real norm = -1.0;

  if (norm < 0.0) {
    norm = 0.0;
    for (i=-Nx(iptr); i<Nx(iptr); i++)
      norm += sinc2(i*Dx(iptr)/spike);
  }
  
  d /= norm;
  for (i=0; i<Nx(iptr); i++) {
    idx = i-ix;
    x = idx*Dx(iptr)/spike;
    x = ABS(x);
    CubeValue(iptr1,i,iy,iz) += d*sinc2(x);
  }
}

int spike_y(d,ix,iy,iz)
     real d;
     int ix,iy,iz;

{
  int j, jdy;
  real y;
  static real norm = -1.0;

  if (norm < 0.0) {
    norm = 0.0;
    for (j=-Ny(iptr); j<Ny(iptr); j++)
      norm += sinc2(j*Dy(iptr)/spike);
  }
  
  d /= norm;
  for (j=0; j<Ny(iptr); j++) {
    jdy = j-iy;
    y = jdy*Dy(iptr)/spike;
    y = ABS(y);
    CubeValue(iptr1,ix,j,iz) += d*sinc2(y);
  }
}

void smooth_it()
{
  real m_min, m_max, brightness, total, d;
  int    i, ix, iy, iz, kounter, idir, nd;
  char   *cp;

  if (noise > 0.0) {
    dprintf(0,"Adding sqr(gausian_noise(%g))\n",noise);
    for (ix=0; ix<Nx(iptr); ix++) {
      for (iy=0; iy<Ny(iptr); iy++) {
	for (iz=0; iz<Nz(iptr); iz++) {
	  CubeValue(iptr1,ix,iy,iz) = sqr(grandom(0.0,noise));
	}
      }
    }
  } else {
    for (ix=0; ix<Nx(iptr); ix++) {
      for (iy=0; iy<Ny(iptr); iy++) {
	for (iz=0; iz<Nz(iptr); iz++) {
	  CubeValue(iptr1,ix,iy,iz) = 0.0;
	}
      }
    }
  }
  
  if (fraction > 0.0) {
    dprintf(0,"Taken fraction %g of intensity in diffraction\n",fraction);
    nd = 0;
    for (ix=0; ix<Nx(iptr); ix++) {
      for (iy=0; iy<Ny(iptr); iy++) {
	for (iz=0; iz<Nz(iptr); iz++) {
	  brightness = CubeValue(iptr,ix,iy,iz);
	  d = fraction*brightness;
	  if (brightness <= cutoff) {
	    CubeValue(iptr1,ix,iy,iz) += d;
	    continue;
	  }
	  CubeValue(iptr,ix,iy,iz) = (1-fraction)*brightness;
	  nd++;
	  d *= 0.5;
	  spike_x(d,ix,iy,iz);
	  spike_y(d,ix,iy,iz);
	}
      }
    }
    dprintf(0,"Diffracted %d pixels above %g\n",nd,cutoff);
  }

  m_min = HUGE;
  m_max = -HUGE;
  total = 0.0;
  
  dprintf (0,"Convolving %s with %d-length beam: \n",dir,lsmooth);
  for (i=0; i<lsmooth; i++)
    dprintf (1," %f ",smooth[i]);
  convolve_cube (Frame(iptr),nx,ny,nz,smooth,lsmooth,1);  /* smooth in X */
  convolve_cube (Frame(iptr),nx,ny,nz,smooth,lsmooth,2);  /* smooth in Y */

  m_max = -HUGE;                      /* determine new min/max */
  m_min =  HUGE;
  for (ix=0; ix<Nx(iptr); ix++)   	
    for (iy=0; iy<Ny(iptr); iy++)
      for (iz=0; iz<Nz(iptr); iz++) {
	CubeValue(iptr,ix,iy,iz) += CubeValue(iptr1,ix,iy,iz);
	brightness = CubeValue(iptr,ix,iy,iz);
	total += brightness;
	m_max = MAX(m_max, brightness);
	m_min = MIN(m_min, brightness);
      }
  MapMin(iptr) = m_min; 		/* update map headers */
  MapMax(iptr) = m_max;
    	
  dprintf (0,"New min and max in map are: %f %f\n",m_min,m_max);
  dprintf (0,"New total brightness/mass is %f\n",total*Dx(iptr)*Dy(iptr));
}


void make_gauss_beam(sdir)
char *sdir;		/* smoothing direction: 'x', 'y' or 'z' */
{
    real sigma, sum, x, fx;
    int    i;

    sigma = gauss_fwhm/2.355;   /* FWHM = 2 sqrt(2 ln(2)) * sigma */
    lsmooth=1;              /* count how many we need for the 'smooth' ARRAY */
    x=0;                    /* minimum is 1, at x=0 */
    sum=1.0;                /* unnormalized value at x=0 */
    if (*sdir == 'x')
        cell = ABS(dx);
    else if (*sdir == 'y')
        cell = ABS(dy);
    else if (*sdir == 'z')
        cell = ABS(dz);
    if (cell==0.0)    
        error("make_gauss_beam: cellsize zero for dir=%c",*sdir);
    dprintf(0,"dir=%c sigma=%g cell=%g\n",*sdir,sigma,cell);

    if (gauss_fwhm == 0.0) {
        warning("No smoothing done in %c",*sdir);
        lsmooth = 1;
        smooth[0] = 1.0;
        return;
    }

    do {
        x += cell;
        lsmooth += 2;
        fx = exp(-0.5*sqr(x/sigma));            
        sum += 2*fx;
    } while ((fx>gauss_cut) && (lsmooth<MSMOOTH));  /* cutoff beam */
    sum = 1.0/sum;
    dprintf (1,"Gaussean beam will contain %d points at cutoff %g \n",
                lsmooth, gauss_cut);
    if ((lsmooth==MSMOOTH) && (fx>gauss_cut)) 
        dprintf(0,"Warning: beam cutoff %g at %d points\n",
                gauss_cut, MSMOOTH);
    for (i=0; i<lsmooth; i++) {
        x = (i-(lsmooth-1)/2)*cell/sigma;
        smooth[i] = exp(-0.5*x*x) * sum;   /* normalize */
    }
}
                

int convolve_cube (a, nx, ny, nz, b, nb, idir)
real *a, b[];
int  nx,ny,nz,nb,idir;
{
    int ix,iy,iz, ier;
    
    if (idir==1)
        for (iy=0; iy<ny; iy++)
        for (iz=0; iz<nz; iz++)
            ier += convolve_x (a,iy,iz,nx,ny,nz,b,nb);
    else if (idir==2)
        for (ix=0; ix<nx; ix++)
        for (iz=0; iz<nz; iz++)
            ier += convolve_y (a,ix,iz,nx,ny,nz,b,nb);
    else if (idir==3)
        for (ix=0; ix<nx; ix++)
        for (iy=0; iy<ny; iy++)
            ier += convolve_z (a,ix,iy,nx,ny,nz,b,nb);
    else
        return 0;
    return ( ier==0 ? 1 : 0 );
}


int convolve_x (a, iy, iz, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, iy, iz;
{
	real c[MSIZE];
	int    ix, jx, kx, offset;
	
	if (nx>MSIZE) {
		warning("convolve_x: MSIZE=%d array too small for input %d\n",
			MSIZE,nx);
                return 0;
	}
	
	for (ix=0; ix<nx; ix++) {
		offset = iz + iy*nz + ix*nz*ny;	/* CDEF */
		c[ix] = *(a+offset);		/* copy array */
		*(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (ix=0; ix<nx; ix++) 
	    for (jx=0; jx<nb; jx++) {
	        kx = ix + jx - (nb-1)/2;
		if (kx>=0)
		    if (kx<nx) {
                        if (Qbad && c[ix]==bad) continue;
		        *(a+iz+iy*nz+kx*ny*nz) += b[jx]*c[ix];
		    } else
			continue;
	    }
	return 1;
}

int convolve_y (a, ix, iz, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, ix, iz;
{
	real c[MSIZE];
	int    iy, jy, ky, offset;
	
	if (ny>MSIZE) {
		warning("convolve_y: MSIZE=%d array too small for input %d\n",
			MSIZE,ny);
                return 0;
	}
	
	for (iy=0; iy<ny; iy++) {
	    offset = iz+iy*nz+ix*nz*ny;     /* CDEF */
	    c[iy] = *(a+offset);		/* copy array */
	    *(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (iy=0; iy<ny; iy++) 
    	    for (jy=0; jy<nb; jy++) {
	        ky = iy + jy - (nb-1)/2;
		if (ky>=0)
		    if (ky<ny){
                        if (Qbad && c[iy]==bad) continue;
		        *(a+iz+ky*nz+ix*ny*nz) += b[jy]*c[iy];
		    } else
			continue;
	    }
	return 1;
}

int convolve_z (a, ix, iy, nx, ny, nz, b, nb)
real *a, b[];
int    nx, ny, nz, nb, ix, iy;
{
	real c[MSIZE];
	int    iz, jz, kz, offset;
	
	if (nz>MSIZE) {
		warning("convolve_z: MSIZE=%d array too small for input %d\n",
			MSIZE,nz);
		return 0;
	}
	
	for (iz=0; iz<nz; iz++) {
            offset = iz+iy*nz+ix*nz*ny;         /* CDEF */
            c[iz] = *(a+offset);		/* copy array */
            *(a+offset) = 0.0;		/* reset for accumulation */
	}
		
	for (iz=0; iz<nz; iz++) 
            for (jz=0; jz<nb; jz++) {
		kz = iz + jz - (nb-1)/2;
		if (kz>=0)
	            if (kz<nz) {
                        if (Qbad && c[iz]==bad) continue;
			*(a+kz+iy*nz+ix*nz*ny) += b[jz]*c[iz];
		    } else
			continue;
	    }
	return 1;
}


real sinc2(x)
     real x;
{
  real y;
  if (x==0.0) return 1.0;
  if (x>PI)   return 0.0;
  y = sin(x)/x;
  dprintf(1,"sinc2: %g %g\n",x,y);
  return y*y;
}
