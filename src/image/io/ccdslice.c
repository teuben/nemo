/* 
 *	CCDSLICE: take slices from a cube
 *
 *	 6-may-95   V1.0 created 	Peter Teuben
 *      29-dec-01   V1.0a   also compute MapMin/Max
 *      27-feb-2021 V1.2    zslabs= and zscale implemented
 *       1-may-2022 V1.3    fix WCS for sampling
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>

#include <image.h>

string defv[] = {	/* keywords + help string for user interface */
    "in=???\n       Input filename (image)",
    "out=???\n      Output filename (image)",
    "zvar=z\n       Slice variable (x,y,z)",
    "zrange=\n      Slices to select (1..n)",
    "zslabs=\n      Zmin,Zmax pairs in WCS to select",
    "zscale=1\n     Scaling applied to zslabs",
    "select=t\n     Select the planes for output (t) or de-select those (f)",
    "VERSION=1.3\n  1-may-2022 PJT",
    NULL,
};

string usage="takes slices from a cube";


#define X_SLICE 0
#define Y_SLICE 1
#define Z_SLICE 2

void ax_copy(imageptr i0, imageptr i1);
void slice(imageptr i, imageptr o, int mode, int *planes);
int  get_planes(int nslabs, real *slabs, int *planes, int na, real aref, real amin, real da, real ascale);

void nemo_main(void)
{
    imageptr iptr=NULL, optr=NULL;
    stream instr, outstr;
    int *planes, mode, i;
    real *slabs, zscale;
    string zvar;
    int nx, ny, nz, maxplane;
    bool Qsel = getbparam("select");

    instr = stropen (getparam("in"),"r");
    read_image (instr,&iptr);            
    strclose(instr);                     

    if (Axis(iptr)) warning("axis=1 not fully supported yet");
    if (!Qsel) error("select=f not implemented yet");

    zvar = getparam("zvar");
    if (streq(zvar,"x")) {
        mode = X_SLICE;
        nx = Ny(iptr);
        ny = Nz(iptr);
        nz = Nx(iptr);
	maxplane = nx;
    } else if (streq(zvar,"y")) {
        mode = Y_SLICE;
        nx = Nx(iptr);
        ny = Nz(iptr);
        nz = Ny(iptr);
	maxplane = ny;
    } else if (streq(zvar,"z")) {
        mode = Z_SLICE;
        nx = Nx(iptr);
        ny = Ny(iptr);
        nz = Nz(iptr);
	maxplane = nz;	
    } else 
        error("Illegal slice axis %s, must be x,y,z",zvar);
    warning("slice mode %d", mode);

    planes = (int *)  allocate(maxplane * sizeof(int));
    slabs  = (real *) allocate(maxplane * sizeof(real));
    
    zscale = getrparam("zscale");

    if (hasvalue("zrange")) {
        nz = nemoinpi(getparam("zrange"),planes,maxplane);
        if (nz<1) 
            error("Bad syntax %d for zrange=%s",nz,getparam("zrange"));
	if (nz>1 && ABS(planes[1]-planes[0]) > 1)
	  warning("Stepping assumed to be regular as %d for WCS correction", planes[1]-planes[0]);
    } else if (hasvalue("zslabs")) {
        nz = nemoinpr(getparam("zslabs"),slabs,maxplane);
        if (nz<1) 
            error("Bad syntax %d for zrange=%s",nz,getparam("zslabs"));
	if (mode == X_SLICE)
	  nz = get_planes(nz, slabs, planes, Nx(iptr), Xref(iptr), Xmin(iptr), Dx(iptr), zscale);
	else if (mode == Y_SLICE)
	  nz = get_planes(nz, slabs, planes, Ny(iptr), Yref(iptr), Ymin(iptr), Dy(iptr), zscale);
	else if (mode == Z_SLICE) 
	  nz = get_planes(nz, slabs, planes, Nz(iptr), Zref(iptr), Zmin(iptr), Dz(iptr), zscale);
	else
	  error("mode not implemented for slabs");
    } else {
        for (i=0; i<nz; i++)
            planes[i] = i+1;
    }
    create_cube(&optr,nx,ny,nz);
    ax_copy(iptr,optr);
    slice(iptr,optr,mode,planes);

    outstr = stropen(getparam("out"),"w");
    write_image(outstr,optr);
    strclose(outstr);
}


/* from miriad::imbin WCS correction:
 *  crpixo(i) = 1 + (crpixi(i)-blc(i))/bin(1,i)
 *  cdelto(i) = bin(1,i) * cdelti(i)
 *  crpixo(i) = crpixo(i) - 0.5*(bin(2,i)-1)/bin(2,i)
 *
 *  where bin(1) = step  and bin(2) = width
 *
 */


void slice(imageptr i, imageptr o, int mode, int *planes)
{
    int x, y, z, iz;
    real ival;


    if (mode==X_SLICE) {
        warning("Code for X not converted to fix reference pixel value");
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Nx(i)) 
                error("%d: illegal plane in x, max is %d",z,Nx(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++) {
		ival = CubeValue(i,z,x,y);
                CubeValue(o,x,y,iz) = ival;
	      }
        }
        Namex(o) = Namey(i);
        Namey(o) = Namez(i);
        Namez(o) = Namex(i);
        Xmin(o) = Ymin(i);  // fix
        Ymin(o) = Zmin(i);
        Zmin(o) = Xmin(i);
        Dx(o) = Dy(i);      // fix
        Dy(o) = Dz(i);
        Dz(o) = Dx(i);
    } else if (mode==Y_SLICE) {
        warning("Code for Y not converted to fix reference pixel value");      
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Ny(i)) 
                error("%d: illegal plane in y, max is %d",z,Ny(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++) {
		ival = CubeValue(i,x,z,y);
                CubeValue(o,x,y,iz) = ival;
	      }
        }
        Namex(o) = Namex(i);
        Namey(o) = Namez(i);
        Namez(o) = Namey(i);
        Xmin(o) = Xmin(i);
        Ymin(o) = Zmin(i);
        Zmin(o) = Ymin(i);
        Dx(o) = Dx(i);
        Dy(o) = Dz(i);
        Dz(o) = Dy(i);
    } else if (mode==Z_SLICE) {
        for(iz=0; iz<Nz(o); iz++) {
            z = planes[iz]-1;
            if (z<0 || z>Nz(i)) 
                error("%d: illegal plane in z, max is %d",z,Nz(i));
            for (y=0; y<Ny(o); y++)
	      for (x=0; x<Nx(o); x++)
                CubeValue(o,x,y,iz) = CubeValue(i,x,y,z);
        }
	real width_step = planes[1]-planes[0];
	dprintf(0,"Z_SLICE: %d %g\n",planes[0],width_step);
	Zref(o) = (Zref(i)-planes[0])/width_step;
	Dz(o) = Dz(i) * width_step;
	Zref(o) = Zref(o) - 0.5*(width_step - 1.0)/width_step;
    }
    minmax_image(o);
}

/*
 *  convert slabs to 1..N based index array
 *
 */

int get_planes(int nslabs, real *slabs, int *planes, int na, real aref, real amin, real da, real ascale)
{
  int i0, i1, i, j, nz = 0;
  real axis_min, axis_max;

  dprintf(0,"Axis: %d %g  %g %g\n",na,aref,amin,da);
  axis_min =       -aref*da + amin;
  axis_max = (na-1-aref)*da + amin;
  dprintf(0,"Axis min/max: %g %g\n",axis_min,axis_max);
  
  if (nslabs%2) error("Not an even number of slab values");

  for (i=0; i<nslabs; i+=2) {
    i0 = (int) ((ascale*slabs[i]  -amin)/da + aref);
    i1 = (int) ((ascale*slabs[i+1]-amin)/da + aref);
    if (i0<0) {
      warning("%d < 0: ",i0);
      i0=0;
    }
    if (i1>=na) {
      warning("%d >= %d: ",i1,na);
      i1=na-1;
    }
    dprintf(0,"slab %g %g -> %d %d\n",slabs[i],slabs[i+1],i0,i1);
    if (i0<i1)
      for (j=i0; j<=i1; j++)
	planes[nz++] = j+1;    
    else
      for (j=i1; j<=i0; j++)
	planes[nz++] = j+1;
  }
  dprintf(0,"Total planes: %d\n",nz);
  return nz;
}
  
void ax_copy(imageptr i0, imageptr i1)
{
  Dx(i1) = Dx(i0);
  Dy(i1) = Dy(i0);
  Dz(i1) = Dz(i0);
  Xmin(i1) = Xmin(i0);
  Ymin(i1) = Ymin(i0);
  Zmin(i1) = Zmin(i0);
  Xref(i1) = Xref(i0);
  Yref(i1) = Yref(i0);
  Zref(i1) = Zref(i0);
  Beamx(i1) = Beamx(i0);
  Beamy(i1) = Beamy(i0);
  Beamz(i1) = Beamz(i0);
  if (Namex(i0))  Namex(i1) = strdup(Namex(i0));
  if (Namey(i0))  Namey(i1) = strdup(Namey(i0));
  if (Namez(i0))  Namez(i1) = strdup(Namez(i0));
  Axis(i1) = Axis(i0);
}
