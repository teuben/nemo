/* bodyext.h - ...
                orbitsegment, orbitsegmentptr,
                baseinfo, D1, D2, D3, AForce, AJerk,
	        AMass, Now, Old, polyinfo, APos, APos_now, T0, T1, T2, T3,
	        timesegments, ATimestep, trakinfo, AVel, AVel-now, aarseth0,
	        base_update, firstpass, get_force, integrate, main<TOOLBOX>, 
	        mkaarseth0, next_body, output, poly_update, pred_low_order,
	        pred_high_order, secondpass, start, thirdpass, timestep_update,
	        trak_update, trans_from_aarseth, trans_to_aarseth */

/*
 *  bodyext.h: for bodies in orbit-extrapolation representation (Aarseth type)
 *
 *      Feb. 1988  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

/*-----------------------------------------------------------------------------
 *  baseinfo  --  basic information about a particle to determine its dynamics
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  base_mass;	  /* mass of particle */
    real  base_t0;	  /* time of most recent update of particle position */
    real  base_pos[NDIM];	  /* position */
    real  base_vel[NDIM];	  /* velocity */
    } baseinfo;

/*-----------------------------------------------------------------------------
 *  trakinfo  --  extra information about a particle, for orbit extrapolation
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  trak_f[NDIM];	         /* force on particle */
    real  trak_fdot[NDIM];       /* jerk (time derivative of the force) */
    } trakinfo;

/*-----------------------------------------------------------------------------
 *  polyinfo  --  polynomial approximation to a particle orbit
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  poly_d1[NDIM];   /* coefficients of the 3rd-order force polynomial */
    real  poly_d2[NDIM];
    real  poly_d3[NDIM];
    real  poly_t1;	   /* time of second-most-recent update of particle */
    real  poly_t2;	   /*  "   "  third-   "     "      "    "    "     */
    real  poly_t3;	   /*  "   "  fourth-  "     "      "    "    "     */
    } polyinfo;

/*-----------------------------------------------------------------------------
 *  orbitsegment, orbitsegmentptr  --  full information about a particle
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    baseinfo  part_baseinfo;
    trakinfo  part_trakinfo;
    polyinfo  part_polyinfo;
    real  part_timestep;	/* time interval until next update required */
    real  part_pos_now[NDIM];	/* secondary position coordinates */
    real  part_vel_now[NDIM];	/* secondary velocity coordinates */
    } orbitsegment, *orbitsegmentptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual (sub)components from an orbitsegment
 *  note: in order to distinguish between similar macros in the main part of 
 *        newton0, a capital A is prefixed to some of the macros below, where
 *        A stands for Aarseth's method of orbit extrapolation.
 *-----------------------------------------------------------------------------
 */
#define  AMass(bodyptr)       ((bodyptr)->part_baseinfo.base_mass)
#define  APos_now(bodyptr)    (&(bodyptr)->part_pos_now[0])
#define  AVel_now(bodyptr)    (&(bodyptr)->part_vel_now[0])
#define  APos(bodyptr)        (&(bodyptr)->part_baseinfo.base_pos[0])
#define  AVel(bodyptr)        (&(bodyptr)->part_baseinfo.base_vel[0])
#define  AForce(bodyptr)      (&(bodyptr)->part_trakinfo.trak_f[0])
#define  AJerk(bodyptr)       (&(bodyptr)->part_trakinfo.trak_fdot[0])
#define  D1(bodyptr)         (&(bodyptr)->part_polyinfo.poly_d1[0])
#define  D2(bodyptr)         (&(bodyptr)->part_polyinfo.poly_d2[0])
#define  D3(bodyptr)         (&(bodyptr)->part_polyinfo.poly_d3[0])
#define  ATimestep(bodyptr)   ((bodyptr)->part_timestep)
#define  T0(bodyptr)         ((bodyptr)->part_baseinfo.base_t0)
#define  T1(bodyptr)         ((bodyptr)->part_polyinfo.poly_t1)
#define  T2(bodyptr)         ((bodyptr)->part_polyinfo.poly_t2)
#define  T3(bodyptr)         ((bodyptr)->part_polyinfo.poly_t3)

/*-----------------------------------------------------------------------------
 *  timesegments  --  contains old and new time intervals used in polynomials
 *-----------------------------------------------------------------------------
 */
typedef  struct
    {
    real  time_old[4];	       /*  old time intervals,
				*  in Aarseth's (1985), eq. (3) notation :
				*  _old[k] = t sub k prime = t sub 0 - t sub k
				*/
    real  time_now[4];         /*  new time intervals,
				*  in Aarseth's (1985), eq. (3) notation :
				*  _now[k] = t - t sub k
				*/
    } timesegments;

/*-----------------------------------------------------------------------------
 *  macros to extract individual (sub)components from a timesegments
 *-----------------------------------------------------------------------------
 */
#define  AOld(timeptr)      (&(timeptr)->time_old[0])
#define  ANow(timeptr)      (&(timeptr)->time_now[0])

/* end of: bodyext.c */
