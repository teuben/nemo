/*
 * MKOP73:   set up a disk using original 1973 Ostriker-Peebles method
 *
 *	  -xxx-93	created as 'MKOP', cloned off mkexpdisk		JAN
 *      16-jul-93       tenergy         				JAN
 *	29-oct-93       optional output of phi/acc/aux/
 *                      added a bit of $real$ history
 *			fixed factor 10 error in vavg !!!
 *			added final velocity rescaling			PJT
 *	 3-nov-93	added zerocm=; implemented options=acc,phi also
 *                      implemented mode= to some extent
 *			Last bug (vel. rescaling) removed		PJT
 *	 4-nov-93	Renaming some parameters for formal release	PJT
 *			named program 'MKOP73'
 *	30-nov-95	Fixed system-specific function call (irint())	JAN
 *	27-mar-97  1.0d SINGLEPREC fixes, but not finished yet		pjt
 *
 */


#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>
#include <potential.h>
#include <moment.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {
    "out=???\n		  output snapshot file name ",
    "nbody=300\n	  usual number of particles ",
    "mdisk=300.0\n	  total disk mass",
    "rcut=1.00\n	  outer cutoff radius",
    "sigma=0.454\n	  scaling factor sigma from Toomre's Q parameter ",
    "nring=10\n           number of particles per ring;should divide nbody",
    "nvel=10\n            number of measurements used to determine circ. vel.",
    "nderiv=10\n          number of rings used to determine dv/dr",
    "seed=0\n		  usual random number seed [0=time since 1970]",
    "eps=0.05\n           softening for force calculation",
    "zerocm=t\n		  Re-center at center of mass?",
    "mode=1\n             creation mode: 1=OP73 2=theoretical rotation curve",
    "options=\n           Additional output (phi,acc,rot,sig)",
    "potname=op73\n       OP73 potential (5NEMO) [or any other]",
    "potpars=0,0\n	  .. with potential parameters",
    "potfile=\n		  .. and optional potential datafile name",    
    "headline=\n	  text headline for output",
    "VERSION=1.0d\n	  27-mar-97 PJT",
    NULL,
};

string usage = "Set up a Mestel disk using 1973 Ostriker & Peebles' method";

#define MAXRING   1000

local real rcut, mdisk, tsig, c;

local int  nbody,nring,nderiv,nvel;
local int  cmode;
local bool Qphi, Qacc, Qrot, Qsig;
local bool zerocm;
local int  mode;
local Body *disktab;

local proc potential;

local real rotcur(real);

extern double xrandom(double,double), grandom(double,double);

nemo_main()
{
    int seed;
    string options;
    bool scanopt();

#ifdef SINGLEPREC
	warning("this program probably doesn't work in SINGLEPREC");
#endif

    potential = get_potential(getparam("potname"),
                    getparam("potpars"), getparam("potfile"));
    rcut = getdparam("rcut");
    mdisk = getdparam("mdisk");
    tsig = getdparam("sigma");
    c = getdparam("eps");
    nbody = getiparam("nbody");
    nring = getiparam("nring");
    nvel  = getiparam("nvel");
    nderiv = getiparam("nderiv");
    set_xrandom(getiparam("seed"));
    zerocm = getbparam("zerocm");
    options = getparam("options");
    mode = getiparam("mode");
    Qphi = scanopt(options,"phi");
    Qacc = scanopt(options,"acc");
    Qrot = scanopt(options,"rot");	/* stores in Aux */
    Qsig = scanopt(options,"sig");	/* stores in Acc's */
    if (nbody%nring != 0)
        error("nbody(%d) must be multiple of nring(%d)",nbody, nring);
    if (nbody/nring > MAXRING)
        error("Too many rings (MAXRING=%d)",MAXRING);

    makedisk();
    writesnap(getparam("out"), getparam("headline"));
}


makedisk()
{
    Body *bp, *tp;
    int i,j,k, iring, izone, nzero=0, ndim=NDIM, ncheat=0;
    real mdsk_i, theta_i, vcir_i, omega, Aoort, kappa,t2,dr,dtheta;
    real mu, sig_r, sig_t, sig_z, vrad_i, veff_i, rbase, thetabase;
    real T0[MAXRING], T1[MAXRING], vavg[MAXRING], deriv[MAXRING], f;
    real nvsq,vsqavg,t=0.0,pot_i,dist,nu,*vorb_i, *rad_i;
    real dtemp,dt,dnder,dnbody,mpoint;
    vector acc_i, r_i, rb_i, dvec, dacc_i;
    double pos_d[NDIM], acc_d[NDIM], pot_d, time_d;
    Moment m;


    disktab = (Body *) allocate(nbody * sizeof(Body));
    rad_i   = (real *) allocate(nbody * sizeof(real));
    vorb_i  = (real *) allocate(nbody * sizeof(real));
    bp=disktab;
    for (i=0; i<MAXRING; i++)
        vavg[i] = T0[i] = T1[i] = 0.0;
    dnbody = (double) nbody;
    dnder = (double) nderiv;
    mpoint = mdisk / dnbody;

    /*CTEX
     *  \centerline{Initial Values in the Standard OP73 Model}
     *  \bigskip
     *     
     *  A particular model called model 1 is taken as the ``standard''
     *  to which we compare the results of varying any of the parameters.
     *  The initial values and other parameters for this standard model
     *  are described here.
     *
     *  The general parameters are: (i) particle number $N = 300$; (ii)
     *  disk radius $R = 1$; (iii) cutoff radius $c = 0.05$;
     *  (iv) timestep $\Delta t = 0.001$; and (v) halo $M_H = 0$.
     *
     *	The initial surface density of points $\Sigma(r)$ varies as
     *	$r^{-1}$. This is achieved by distributing the points uniformly
     *  in the radius interval $0 < r < R$, as follows. The disk is divided
     *  into $N/10$ rings, with width in radius $\Delta r = 10 R/N$ each,
     *  and into 10 equal pie-shaped radial segments ($\Delta\phi = 36^o$).
     *  One point is placed in each of the $N$ cells, with the radius
     *  $r$ and longitude $\phi$ randomly chosen over the range of arguments
     *  for the cell. The points are placed in a flat disk. The thickness
     *  of the disk is established by assigning random velocities normal
     *  to the disk.
     */
     
    dr= rcut * (double) nring / (double) nbody;     /* ring width */
    dtheta= TWO_PI / (double) nring;             /* segment width (radians) */
    for (i=0; i < nbody/nring; i++) {        /* set positions for all stars */
        rbase=(double) i * dr;
        for(j=0;j<nring;j++){
            Mass(bp) = mpoint;
            thetabase=(double) j * dtheta;
            rad_i[i] = rbase + xrandom(0.0,dr);
            theta_i = thetabase + xrandom(0.0, dtheta);

	    Pos(bp)[0] = rad_i[i] * sin(theta_i);	/* assign positions */
	    Pos(bp)[1] = rad_i[i] * cos(theta_i);
	    Pos(bp)[2] = 0.0;

            bp++;
        }
    }

    /*CTEX
     *  The first step in assigning initial velocities is to estimate
     *	the angular velocity needed to hold each particle in circular
     *	orbit. Each point is temporarely rotated by an angle $\phi_i$
     *  equal to a random fraction of $360^o$, and
     *  $\vec{r}_i . \vec{a}_i$ computed.
     *  The average value of this quantity yields the desired mean
     *  speed if the radial positions of the particles do not vary
     *  and if the angular positions are not correlated. The initial
     *  velocity of each particle is directed in the plane perpendicular
     *  to the radial vector, the magnitude being
     *  $<-\vec{r}_i.\vec{a}_i>^{1/2}$.
     */

    if (mode==1) {
      dprintf(1,"Estimating circular velocities (slow)\n");
      for (bp=disktab, i=0; i < nbody; bp++, i++) {  /* loop to find circ. velo.*/
	vsqavg=0.0;
	ABSV(rad_i[i],Pos(bp));			/* save rad_i array values */
	for (j=0; j<nvel;j++){			/* sample 'nvel' times */
	    if (j==0) {
                r_i[0] = pos_d[0] = Pos(bp)[0];
                r_i[1] = pos_d[1] = Pos(bp)[1];
                r_i[2] = pos_d[2] = Pos(bp)[2];
	    } else {
	    	dt = xrandom(0.0,TWO_PI); /* rotate by random angle */
	        r_i[0] = pos_d[0] = rad_i[i]*sin(dt);  	/* new sample pos */
      	        r_i[1] = pos_d[1] = rad_i[i]*cos(dt);
	        r_i[2] = pos_d[2] = 0.0;
            }
	    (*potential)(&ndim,pos_d,acc_d,&pot_d,&time_d);  /* include pot */
            SETV(acc_i,acc_d);
	    for (tp=disktab,k=0;k<nbody;tp++,k++){   /* sum over all others */
		if (k!=i){
		   DISTV(dist,Pos(tp),r_i);
		   dist=sqrt(dist*dist+c*c);
		   f=Mass(tp)/dist;
		   if (Qphi && j==0) pot_i -= f;
		   f /= dist*dist;
		   SUBV(dvec,Pos(tp),r_i);
		   MULVS(dacc_i,dvec,f);
		   SADDV(acc_i,dacc_i);
		}
	    }
	    DOTVP(nvsq,r_i,acc_i);
	    vsqavg -= nvsq;
	    if (Qphi && j==0) Phi(bp) = pot_i;
	    if (Qacc && j==0) SETV(Acc(bp),acc_i);
	}
	if (vsqavg > 0.0) {
	    vorb_i[i]=sqrt(vsqavg/nvel);
        } else {
            ncheat++;
      	    vorb_i[i]=0.0;
      	}

	if (Qrot) Aux(bp) = vorb_i[i];
      }
      if (ncheat > 0) {
	warning("%d/%d particles with too little orbital motion",ncheat,nbody);
	ncheat = 0;
      }
    } else {
      dprintf(1,"Using expected rotation curve to set circular velocities\n");
      for (bp=disktab, i=0; i < nbody; bp++, i++) {  /* loop to set circ. velo.*/
	ABSV(rad_i[i],Pos(bp));			/* save rad_i array values */
	vorb_i[i] = rotcur(rad_i[i]);
      }
    }


    /*CTEX
     *  The above procedure introduces a scatter of approximately
     *  20 percent in the circular velocity of particles in the
     *  same ring. We increase this scatter by adding a velocity
     *  dispersion designed to fit the Toomre (1964) criterion
     *  for stability against the development of small-scale
     *  irregularities. The added velocity components are drawn
     *  from random normal distributions with standard deviations
     *  $\sigma_r, \sigma_\theta, \sigma_z$. The standard deviation
     *  of the velocity added in the radial direction is the smaller
     *  of
     *  $$
     *     \sigma_r^T = \sigma N \nu^{-1}, \sigma_r^S = 0.4\nu       \eqno(7)
     *  $$
     *  where
     *  $$
     *          \nu \equiv v(1 + d \ln{v}/d \ln{r})^{1/2}
     *  $$
     *  Here $v(r)$ is the mean particle speed as a function of position
     *  obtained from the first step, and $N$ the number of mass
     *  points in the model. The first equation is based on Toomre's
     *  condition (1964) for stability against growth of irregularities
     *  of small scale.
     *
     *  In our standard integration we choose the constant
     *  $$
     *          \sigma = 0.454                                       \eqno(8)
     *  $$
     *  20 percent larger than the minimum for stability
     *  according to Toomre. The quantity $\sigma_r^T$ gets very large
     *  at small $r$, so the condition $\sigma_r \leq \sigma_r^T$
     *  defined by the second equation is added to assure that the 
     *  orbits are at least roughly circular.
     *
     *  The derivates of $v$ in equation (7) are obtained from $v$
     *  averaged in intervals $\Delta r = 0.1$, and $\sigma_r$ is
     *  evaluated for nine zones in radius. Finally, we reduce
     *  $\sigma_r$ in the outermost zone by a factor of 2 in order to
     *  reduce the expansion of the edges of the disk. The $\sigma_r$
     *  for the radial velocity added to any particle is the value for
     *  the zone in which the particle finds itself.
     *
     *  For $\sigma_\theta$ we take the equilibrium expression required
     *  for steady epicyclic motions:
     *  $$
     *       \sigma_\theta = {\sigma_r \over 2^{1/2}}
     *              (1 + d \ln{v}/d \ln{r})^{1/2}                   \eqno(9)
     *  $$
     *  The dispersion normal to the plane is taken to be
     *  $$
     *              \sigma_z = \sigma_\theta                        \eqno(10)
     *  $$
     */

    for (bp=disktab, i=0; i < nbody; bp++, i++) {       /* average in zones */
	izone = (int) (floor(rad_i[i]/rcut * dnder));
	vavg[izone] += vorb_i[i] / (dnbody/dnder);
    }
    
    dprintf(1,"# In %d zones: Mean velocity and it's derivative:\n",nderiv);
    for (i=0 ; i< nderiv ; i++) {
        if (i<nderiv-1)
    	    deriv[i]=(vavg[i+1]-vavg[i])*dnder/rcut;
    	else
    	    deriv[i]=0.0;
        dprintf(1,"%d %g %g\n",i+1,vavg[i],deriv[i]);
    }
     

    for (bp=disktab,i=0;i<nbody;bp++,i++){      /* set raw vel.disp */

	    if (rad_i[i]/rcut < 1.5*rcut/ dnder) {
		dtemp=deriv[0];
	    } else if (rad_i[i]/rcut > ( dnder-1.5)*rcut/dnder) {
		dtemp=deriv[nderiv-2];
	    } else
                dtemp=deriv[(int) (rad_i[i]/rcut* dnder - 1.0)];
            /* !!! dtemp might be < 0 */
            if (dtemp<-0.5) {   /* might even allow up to -1 ??? */
                dtemp= -0.5;   /* don't allow > kepler falloff */
                ncheat++;
            }
	    nu=vorb_i[i]*sqrt(1+rad_i[i]/vorb_i[i]*dtemp); /* smoother ??? */

	    if (tsig * mdisk/nu < 0.4 * nu) {
		sig_r = tsig * mdisk/nu;
	    } else
                sig_r = 0.4*nu;

            if (rad_i[i]/rcut > (dnder-1.0)*rcut/dnder)
                sig_r /= 2.0;            /* outmost zone: not as hot */
                
	    sig_t=sig_r/(sqrt(2.0)*vorb_i[i])*nu;
	    sig_z=sig_t;
	    vorb_i[i] += grandom(0.0,sig_t);
	    vrad_i = grandom(0.0,sig_r);
	    Vel(bp)[0]=(vrad_i * Pos(bp)[0] + vorb_i[i] * Pos(bp)[1]) / rad_i[i];
	    Vel(bp)[1]=(vrad_i * Pos(bp)[1] - vorb_i[i] * Pos(bp)[0]) / rad_i[i];
	    Vel(bp)[2] = grandom(0.0,sig_z);
	    if (Qsig) {
                Acc(bp)[0] = sig_r;
                Acc(bp)[1] = sig_t;
                Acc(bp)[2] = sig_z;
            }
    }
    if (ncheat>0) {
        warning("%d/%d orbits with apparent large local rotcur gradient",
                 ncheat/nbody);
        ncheat=0;
    }

    /*CTEX
     *  As a final step, the vector velocity $v_i$ is multiplied
     *  by a factor, separately determined for each of the original
     *  rings, such that the ring has the same total kinetic energy
     *  as before. The purpose of this last correction is to
     *  approximately reestablish the equilibrium between gravitational
     *  and centripetal accelarations that existed before the random
     *  energies were added.
     */

    for (bp=disktab,i=0;i<nbody;bp++,i++) {
        iring = (int) (floor(rad_i[i]/dr));
        T0[iring] += sqr(vorb_i[i]);
        ABSV(nvsq,Vel(bp));
        T1[iring] += sqr(nvsq);
    }
    dprintf(1,"Rescaling velocity factors per ring:\n");
    ini_moment(&m,4);
    for (i=0; i<nbody/nring; i++) {      /* compute scale factors per ring */
        f = sqrt(T1[i] / T0[i]);
        accum_moment(&m,f,1.0);
        dprintf(1,"%d %g %g %g\n",i+1,f,T0[i],T1[i]);
        T1[i] = 1.0/f;
    }
    dprintf(0,"Mean rescaling factor: %g +/- %g\n",
	mean_moment(&m), sigma_moment(&m));

    for (bp=disktab,i=0;i<nbody;bp++,i++) {     /* rescale all velocities */
        iring = (int) (floor(rad_i[i]/dr));
        SMULVS(Vel(bp),T1[iring]);
    }

    if (zerocm)
	centersnap(disktab, nbody);	/* OP73 don't mention this .... */
}

centersnap(btab, nb)
Body *btab;
int nb;
{
    real mtot;
    vector cmphase[2], tmp;
    Body *bp;

    mtot = 0.0;
    CLRV(cmphase[0]);
    CLRV(cmphase[1]);
    for (bp = btab; bp < btab + nb; bp++) {
	mtot = mtot + Mass(bp);
	MULVS(tmp, Phase(bp)[0], Mass(bp));
	ADDV(cmphase[0], cmphase[0], tmp);
	MULVS(tmp, Phase(bp)[1], Mass(bp));
	ADDV(cmphase[1], cmphase[1], tmp);
    }
    MULVS(cmphase[0], cmphase[0], 1.0/mtot);
    MULVS(cmphase[1], cmphase[1], 1.0/mtot);
    for (bp = btab; bp < btab + nb; bp++) {
	SUBV(Phase(bp)[0], Phase(bp)[0], cmphase[0]);
	SUBV(Phase(bp)[1], Phase(bp)[1], cmphase[1]);
    }
}

writesnap(name, headline)
string name;
string headline;
{
    stream outstr;
    real tzero = 0.0;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
	set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    if (Qacc || Qsig) bits |= AccelerationBit;
    if (Qphi) bits |= PotentialBit;
    if (Qrot) bits |= AuxBit;
    put_snap(outstr, &disktab, &nbody, &tzero, &bits);
    strclose(outstr);
}

local real rotcur(real r)
{
    local int ndim=3;
    local real pos[3] = {0,0,0};
    local real tsnap = 0.0;
    real acc[3], pot, vsq;

    /* external potential */
    
    pos[0] = r;
    (*potential)(&ndim,pos,acc,&pot,&tsnap);
    vsq = -r * acc[0];

    /* add theoretical mestel disk (this one is wrong though, just stub) */

    vsq += mdisk/rcut;

    return sqrt(vsq);
}
