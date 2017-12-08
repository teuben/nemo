/*
 * VELFOUR: 
 *
 *	VELFOUR computes the fourier harmonics to an image,
 *      allowing for a center and inclined disk.
 *   
 *   5-dec-2017   Cloned off velfit (and fitting code via snapfour)
 *
 */

#include <nemo.h>
#include <image.h>
#include <spline.h>

string defv[] = {
  "in=???\n       input velocity field",
  "radii=\n       radii of the ring boundaries (Nring+1)",
  "pa=0\n         position angle of disk",
  "inc=45\n       inclination angle of disk",
  "vsys=0\n       systemic velocity at center to be subtracted",
  "center=\n      rotation center (mapcenter if left blank, 0,0=lower left)",

  "frang=0\n      free angle around minor axis (2*frang is the total)",
  "blank=0.0\n    Value of the blank pixel to be ignored",
  "coswt=1\n      power of cos(theta) weighting",
  "mode=vtan\n    Output mode {vtan,vmod,vres,vtan/r,ome,vrad,dv/dr}",

  "cos=0:4:1\n    List of noj-zero cos(m.phi) terms",
  "sin=1:4:1\n    List of non-zero sin(m.phi) terms",
  "amode=t\n      Display sin/cos amps or amp/phase if possible?",

  "out=\n         Optional output map of converted rotation speeds []",
  "tab=\n         Optional output table of radii, velocities etc.",
  "VERSION=0.1\n  5-dec-2017 PJT",
  NULL,
};

string usage="compute fourier harmonics of an image";

string cvsid="$Id$";


#define MAXRING    4096


bool Qden = FALSE;
bool Qout = FALSE;
bool Qtab = FALSE;
imageptr denptr = NULL, velptr = NULL, outptr = NULL;

real rad[MAXRING];
int nrad;

real pixe[MAXRING], flux[MAXRING], vsum[MAXRING], vsqu[MAXRING], wsum[MAXRING];
real vrot[MAXRING], radius[MAXRING], coeff[3*MAXRING];

real pa, inc, vsys, xpos, ypos;
real undf;

    /* outmodes: the order of these is important, see mode= */
    /*    mode:    0    1    2    3      4   5    6     */
string outmodes = "vtan,vmod,vres,vtan/r,ome,vrad,dv/dr";

int string_index(string options, string s)
{
  int i=0;
  string *sa = burststring(options,",");

  while (sa[i]) {
    if (streq(sa[i],s))  {
      freestrings(sa);
      return i;
    }
    i++;
  }
  freestrings(sa);	   
  return -1;
}

/*
 *  -1:     internal
 *   0:     first ring   (rad[0]..rad[1])
 *   nring: last ring
 *   -2:    outside
 */

int ring_index(int n, real *r, real rad)
{
  int i;
  if (rad < r[0]) return -1;
  if (rad > r[n-1]) return -2;
  for (i=0;i<n;i++)
    if (rad >= r[i] && rad < r[i+1]) return i;
  error("ring_index: should never gotten here %g in [%g : %g]",
	rad,r[0],r[n-1]);
}

void nemo_main()
{
  stream denstr, velstr, outstr, tabstr;
  real center[2], cospa, sinpa, cosi, sini, sint, cost, costmin, 
    x, y, xt, yt, r, den, vrot_s, dvdr_s;
  real vr, wt, frang, dx, dy, xmin, ymin, rmin, rmax, fsum, ave, tmp, rms;
  real sincosi, cos2i, tga, dmin, dmax, dval, vmod;
  int i, j, k, nx, ny, ir, nring, nundf, nout, nang, nsum, coswt;
  string outmode;
  int mode = -1;

  warning("PJT CODE NOT FINISHED - this is still bareboned velfit");

  velstr = stropen(getparam("in"),"r");

  read_image(velstr,&velptr);
  nx = Nx(velptr);
  ny = Ny(velptr);

  if (hasvalue("out")) {
    outmode = getparam("mode");
    mode = string_index(outmodes, outmode);
    if (mode < 0) error("Illegal mode=%s [%d], valid:",outmode,mode,outmodes);
    warning("New out= mode mode=%s [%d]",outmode,mode);
    Qout = TRUE;
    outstr = stropen(getparam("out"),"w");
    copy_image(velptr,&outptr);
  } 

  if (hasvalue("tab")) {
    Qtab = TRUE;
    tabstr = stropen(getparam("tab"),"w");
  }

  nrad = nemoinpd(getparam("radii"),rad,MAXRING);
  if (nrad < 2) error("got no rings, use radii=");
  nring = nrad-1;
  if (hasvalue("center")) {
    if (nemoinpd(getparam("center"),center,2) != 2)
      error("not enuf for center=");
    xpos = center[0];
    ypos = center[1];
  } else {
    xpos = (Nx(velptr)-1.0)/2.0;
    ypos = (Ny(velptr)-1.0)/2.0;
  }
  pa = getdparam("pa");
  inc = getdparam("inc");
  vsys = getdparam("vsys");
  undf = getdparam("blank");
  frang = getdparam("frang");
  coswt = getiparam("coswt");

  cospa   = cos(pa*PI/180.0);
  sinpa   = sin(pa*PI/180.0);
  sini    = sin(inc*PI/180.0);
  cosi    = cos(inc*PI/180.0);
  costmin = sin(frang*PI/180.0);
  sincosi = sini*cosi;
  cos2i   = cosi*cosi;
    
  for (i=0; i<nring; i++)
    pixe[i] = flux[i] = vsum[i] = vsqu[i] = wsum[i] = 0.0;
  nundf = nout = nang = 0;

  ymin = Ymin(velptr);
  xmin = Xmin(velptr);
  dx = Dx(velptr);       dx = ABS(dx);    dx = -dx;
  dy = Dy(velptr);       dy = ABS(dy);
  rmin = -nx*dx*10.0;
  rmax = 0.0;

  dmin = dmax = vsys;

  /* loop over the map, accumulating data for fitting process */

  for (j=0; j<ny; j++) {
    y = (j-ypos)*dy;
    for (i=0; i<nx; i++) {
      if (MapValue(velptr,i,j) == undf) {
	nundf++;
	if (Qout) MapValue(outptr,i,j) = undf;
	continue;
      }
      x = (i-xpos)*dx;
      yt = x*sinpa + y*cospa;
      xt = x*cospa - y*sinpa;
      if (Qout) {      /* record the circular speed in diskplane */
	dval = MapValue(velptr,i,j)-vsys;
	if (dval < 0) dval = -dval;
	if (yt==0.0) {
	  dval = undf;
	} else {
	  tga = xt/yt;        /* X and Y are reversed  */
	  dval *= sqrt(cos2i+tga*tga)/sincosi;
	}
	if (mode==0)                      /* mode=vtan */
	  MapValue(outptr,i,j) = dval;
	else if (mode==3 || mode==4) {    /* mode=vtan/r or mode=ome */
	  MapValue(outptr,i,j) = dval/sqrt(sqr(xt/cosi)+sqr(yt));
	} 
      }
      xt /= cosi;
      r  = sqrt(xt*xt+yt*yt);
      rmin = MIN(r,rmin);
      rmax = MAX(r,rmax);
      ir = ring_index(nrad,rad,r);
      dprintf(2,"r=%g ir=%d  (x,y)=%g,%g  (xt,yt)=%g,%g\n",
	      r,ir,x,y,xt,yt);
      if (ir < 0) {
	nout++;
	continue;
      }
      cost = yt/r;
      vr = (MapValue(velptr,i,j)-vsys)/cost/sini;  /* rotation speed */
      den = 1.0;
      wt = den;
      for (k=0; k<coswt; k++)  wt *= cost;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	pixe[ir] += 1.0;
	flux[ir] += den;
	vsum[ir] += wt*vr;
	vsqu[ir] += wt*vr*vr;
	wsum[ir] += wt;
      } else
	nang++;
    }
  }

  /* find the rotation curve */

  for (i=0; i<nring; i++) {
    if (pixe[i] < 1) continue;
    if(wsum[i] != 0)
      vrot[i] = vsum[i]/wsum[i];
    else
      vrot[i] = 0.0;
    /* reset arrays for next round of accumulations */
    wsum[i] = vsum[i] = vsqu[i] = 0.0;
  }

  /* set up rotation curve for differentation and get a spline */
  for (i=0; i<nring; i++) {
    radius[i] = 0.5*(rad[i]+rad[i+1]);
  }
  spline(coeff, radius, vrot, nring);


 /* loop over the map again, computing the residuals */

  for (j=0; j<ny; j++) {
    y = (j-ypos)*dy;
    for (i=0; i<nx; i++) {
      if (MapValue(velptr,i,j) == undf) continue;
      x = (i-xpos)*dx;
      yt =  x*sinpa + y*cospa;
      xt = (x*cospa - y*sinpa)/cosi;
      r  = sqrt(xt*xt+yt*yt);
      rmin = MIN(r,rmin);
      rmax = MAX(r,rmax);
      ir = ring_index(nrad,rad,r);
      if (ir < 0) continue;
      cost = yt/r;
      sint = xt/r;
      if (ABS(cost) < costmin) {
	if (Qout) MapValue(outptr,i,j) = undf;
	continue;
      }
      den = 1.0;
      wt = den;
      for (k=0; k<coswt; k++)  wt *= cost;
      wt = ABS(wt);
      if (ABS(cost) > costmin) {
	vrot_s = seval(r, radius, vrot, coeff, nring);  
	dvdr_s = spldif(r, radius, vrot, coeff, nring);
	vmod = vsys + vrot[ir]*cost*sini;     /* model */
	vr = MapValue(velptr,i,j)-vmod;       /* residual:   vobs-vmod */
	if (Qout) {
	  if (mode==1)                        /* mode=vmod */
	    MapValue(outptr,i,j) = vmod;
	  else if (mode==2)                   /* mode=vres */
	    MapValue(outptr,i,j) = vr;
	  else if (mode==5)                   /* mode=vrad */
	    MapValue(outptr,i,j) = vr/(sint*sini);
	  else if (mode==6)                   /* mode=dv/dr */
	    MapValue(outptr,i,j) = dvdr_s;
	}
	wsum[ir] += wt;
	vsum[ir] += wt*vr;
	vsqu[ir] += wt*vr*vr;
	dprintf(1,"%d %g %g %g %g\n",ir,xt,yt,vr,wt);
      } 
    }
  }

  /* report */

  if (Qtab) fprintf(tabstr,"# r v rms N tmp sqrt(F/N),  ave, vsum, wsum\n");
  fsum = 0.0;
  nsum = 0;
  for (i=0; i<nring; i++) {
    if (wsum[i] == 0.0) continue;
    nsum++;
    r = 0.5*(rad[i] + rad[i+1]);
    tmp = flux[i]/pixe[i];
    ave = vsum[i]/wsum[i];
    rms = vsqu[i]/wsum[i]-ave*ave;
    if (pixe[i]==0 || rms <  0) rms=0.0;
    rms = sqrt(rms);
    fsum += rms*rms;
    if (Qtab) fprintf(tabstr,"%g %g %g %g %g %g  %g %g %g\n",
	   r,vrot[i],rms,pixe[i],tmp,sqrt(fsum/nsum),    ave,vsum[i],wsum[i]);
  }

  dprintf(0,"Nundf=%d/%d Nout=%d Nang=%d (sum=%d)\n",
	  nundf,nx*ny,nout,nang,nout+nundf+nang);
  dprintf(0,"Rmin/max = %g %g\n",rmin,rmax);


  /* write output map(s), if needed */
  if (Qout) {
    dmin = dmax = MapValue(outptr,0,0); 
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
	dval = MapValue(outptr,i,j);
	if (dval > dmax) dmax = dval;
	if (dval < dmin) dmin = dval;
      }
    }
    dprintf(0,"Data min/max = %g %g\n",dmin,dmax);    
    MapMin(outptr) = dmin;
    MapMax(outptr) = dmax;
    write_image(outstr,outptr);
  }


}

#if 0

/*
 * SNAPFOUR.C:   fourier coefficients of an N-body distribution
 *
 *      30-nov-90       V1.0    Created - after Kevin Long's talk       PJT
 *	20-dec-90	V1.0b	added printed header			PJT
 *      17-feb-92       V1.1    added weight=                           PJT
 *	22-feb-92	V1.1b   usage
 *       9-nov-93       V1.2    times=
 *       7-may-02       minor code cleanup
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>  
#include <snapshot/body.h>
#include <snapshot/get_snap.c>

string defv[] = {
    "in=???\n              Input snapshot",
    "radii=0:2:0.1\n       Set of radii denoting edges of cylinders",
    "cos=0:4:1\n	   List of noj-zero cos(m.phi) terms",
    "sin=1:4:1\n	   List of non-zero sin(m.phi) terms",
    "xvar=x\n              X variable",
    "yvar=y\n              Y variable",
    "fvar=vy\n             Fourier Observable to be decomposed",
    "weight=1\n            Weight applied to observable",
    "amode=t\n             Display sin/cos amps or amp/phase if possible?",
    "times=all\n           Snapshots to select",
    "VERSION=1.2b\n        15-jul-04 PJT",
    NULL,
};

string usage = "Fourier coefficients of an N-body distribution";

#define TIMEFUZZ        0.0001  /* tolerance in time comparisons */

#define MAXRAD 513
#define MAXORDER 8

nemo_main()
{
    stream instr;
    string times, vstr;
    Body   *btab = NULL, *bp;
    int    i, n, nbody, bits, nrad, maxorder, tmpi[MAXORDER+1];
    real rad2[MAXRAD], tsnap;
    bool Qcos[MAXORDER+1], Qsin[MAXORDER+1], amode;
    rproc btrtrans(), xproc, yproc, fproc, wproc;

    times = getparam("times");
    nrad = nemoinpr(getparam("radii"),rad2,MAXRAD);     /* get radii */
    for (i=0; i<nrad; i++)
        rad2[i] = sqr(rad2[i]);             /* but actually save the square */

    for (i=0; i<=MAXORDER; i++) {       /* initially set all coefs to false */
        Qcos[i] = FALSE;
        Qsin[i] = FALSE;        /* Qsin[0] also set but never used though */
    }
    n = nemoinpi(getparam("cos"),tmpi,MAXORDER+1);  /* get true cos coefs */
    for (i=0; i<n; i++)
        if (tmpi[i]<0 || tmpi[i]>MAXORDER)
            warning("Illegal value %d for cos= skipped",tmpi[i]);
        else
            Qcos[tmpi[i]] = TRUE;
    n = nemoinpi(getparam("sin"),tmpi,MAXORDER+1);  /* get true sin coefs */
    for (i=0; i<n; i++)
        if (tmpi[i]<0 || tmpi[i]>MAXORDER)
            warning("Illegal value %d for sin= skipped",tmpi[i]);
        else
            Qsin[tmpi[i]] = TRUE;
    for (i=0, maxorder=-1; i<=MAXORDER; i++)
        if (Qsin[i] || Qcos[i]) maxorder=i;
    if (maxorder<0) error("No true sin or cos coefficients supplied");
    xproc = btrtrans(getparam("xvar"));
    yproc = btrtrans(getparam("yvar"));
    fproc = btrtrans(getparam("fvar"));
    wproc = btrtrans(getparam("weight"));
    amode = getbparam("amode");
        
    instr = stropen(getparam("in"), "r");           /* open input file */
    get_history(instr);                         /* get history */
    for (;;) {                          /* loop through snapshots */
        if (!get_tag_ok(instr, SnapShotTag))
                break;                           /* until done */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);
        if (!streq(times,"all") && !within(tsnap,times,0.0001))
            continue;                   /* skip work on this snapshot */
        if ((bits & MassBit) == 0 && (bits & PhaseSpaceBit) == 0) {
            dprintf (2,"Time= %f auto skipping ",tsnap);
            continue;       /* just skip - it maybe diagnostics */
        }
        dprintf (2,"Time= %f ",tsnap);
        snap_four(btab,nbody,tsnap,xproc,yproc,fproc,wproc,
                maxorder,Qcos,Qsin,rad2,nrad,amode);
    }
}


snap_four(btab,nbody,tsnap,xproc,yproc,fproc,wproc,
          maxorder,Qcos,Qsin,rad,nrad,amode)
Body *btab;                 /* pointer to snspshot with nbody Bodie's */
real rad[];                 /* radii for shells */
real tsnap;                 /* time of snapshot */
rproc xproc, yproc, fproc;  /* procedures to compute radius, angle and f */
rproc wproc;                /* weight factor per data point */
bool Qcos[], Qsin[];        /* designate if coef to be used */
int nbody, maxorder, nrad;     
bool amode;                 /* TRUE=amps only FALSE=amp+phase if all available */
{
    real   th, r2, v, rsum, vsum, cosk, amp, pha, radius, w;
    int    i,k,m,cnt,dim,ip;
    Body *bp;
    real mat[2*(MAXORDER+1)*(MAXORDER+1)],vec[2*(MAXORDER+1)];
    real sol[2*(MAXORDER+1)], a[2*(MAXORDER+1)+1];
    permanent bool first=TRUE;

    for (m=0, dim=0; m<=maxorder; m++)  /* count dimension of matrix needed */
        if (Qcos[m]) dim++;
    for (m=1; m<=maxorder; m++)
        if (Qsin[m]) dim++;
    if (dim==0) {
       warning("snap_four: dim=0");
       return;
    }
    dprintf(1,"snap_four: maxorder=%d dim=%d\n",maxorder,dim);
    if (!amode) {       /* check if OK to do phases and amplitudes */
        for (m=1; m<=maxorder; m++) {
            if (Qcos[m] && !Qsin[m]) amode=TRUE;
            if (!Qcos[m] && Qsin[m]) amode=TRUE;
        }
        if (amode) 
            warning("amode=f requested, but missing cos/sin terms");
    }
    if (first) {
        print_header(maxorder,Qcos,Qsin,amode);
        first = FALSE;
    }

    for (i=1; i<nrad; i++) {            /* foreach ring */
        lsq_zero(dim,mat,vec);          /* reset accum. matrix and vector */
        cnt = 0;                        /* count points in this ring */
        for (ip=0, bp=btab; ip<nbody; ip++, bp++) {  /* loop for all bodies */
            w = wproc(bp,tsnap,ip);
            r2 = sqr(xproc(bp,tsnap,ip)) + sqr(yproc(bp,tsnap,ip));
            if (r2<rad[i-1] || r2>rad[i])          /* if not in ring:  */
                continue;                           /* skip this particle */
            cnt++;
            th = atan2(yproc(bp,tsnap,ip) , xproc(bp,tsnap,ip));
            k=0;                            /* always count how many coefs */
            for(m=0; m<=maxorder; m++) {     /* cos(m.theta) */
                if (!Qcos[m]) continue;
                a[k++] = cos(m*th);
            }
            for(m=1; m<=maxorder; m++) {     /* sin(m.theta) */
                if (!Qsin[m]) continue;
                a[k++] = sin(m*th);
            }
            a[k] = fproc(bp,tsnap,ip);
            dprintf(1,"adding %d: r^2=%g th=%g, fvar=%g wt=%g\n",
			      cnt,r2,th*180/PI,a[k],w);
            if (k!=dim) error("snapfour: Counting error dim=%d k=%d",dim,k);
            lsq_accum(dim,mat,vec,a,w);  /* accumulate for LSQ normal matrix */
        } /* bp */
    
        radius = 0.5*(sqrt(rad[i-1])+sqrt(rad[i]));
	if (cnt<dim) {
            dprintf(0,"radius %g has %d points: skipped\n",radius,cnt);
            continue;
        }
        lsq_solve(dim,mat,vec,sol);
        printf("%g %d",radius,cnt);
	if (amode) {				/* only print amplitudes */
            for (k=0; k<dim; k++)
                printf(" %g",sol[k]);
        } else {			   /* figure out amp/phase stuff */
            k=0;        /* pointer which 'sol' has been printed */
            if (Qcos[0])		/* if offset wanted, print it now */
                printf(" %g",sol[k++]);
            for( ; k < (dim+1)/2; k++) {	/* go over all amp/phase */
                amp = sqrt(sqr(sol[k]) + sqr(sol[k+dim/2]));
                pha = atan2(sol[k+dim/2],sol[k]) * 180/PI;
                printf(" %g %g",amp,pha);
            }
	}
        printf("\n");
    } /* i */
}

print_header(maxorder,Qcos,Qsin,amode)
int maxorder;
bool Qcos[], Qsin[], amode;
{
    int m;
    
    dprintf(0,"<R> N ");
    if (amode) {
        for(m=0; m<=maxorder; m++)
            if (Qcos[m]) dprintf(0,"A%d ",m);
        for(m=1; m<=maxorder; m++)
            if (Qsin[m]) dprintf(0,"B%d ",m);
    } else {
        if (Qcos[0]) dprintf(0,"C0 ");
        for(m=1; m<=maxorder; m++)
            if (Qsin[m]) dprintf(0,"C%d P%d ",m,m);
    }
    dprintf(0,"\n");
}


#endif
