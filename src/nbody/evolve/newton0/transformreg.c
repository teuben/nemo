/* transformreg.c - mk_reg_system, nonreg_to_reg, reg_to_nonreg, setup_reg */

/*
 *  transformreg.c:  for transforming to & from regularized coordinates
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static void noncom_to_com(bodyptr bodies, int npart, real com_pos[], real com_vel[]);
static void com_to_noncom(bodyptr bodies, int npart, real com_pos[], real com_vel[]);

/*-----------------------------------------------------------------------------
 *  setup_reg  --  creates the initial pair-representation in 4-dimensional
 *                 regularized coordinates of the standard newton0 system
 *                 (which has already received 4-dimensional vectors with 4th
 *                 components zero).
 *-----------------------------------------------------------------------------
 */
systptr  setup_reg(nonreg_sys)
systptr  nonreg_sys;
    {
    systptr  reg_sys;

    reg_sys = mk_reg_version(nonreg_sys);

    return(reg_sys);
    }

/*-----------------------------------------------------------------------------
 *  mk_reg_version  --  creates a new version of a system, in a
 *                      regularized representation:
 *                      translates from the standard newton0 3-dimensional 
 *                      representation, embedded in 4 dimensions by having the
 *                      fourth component identically zero in all vectors, to a
 *                      new presentation in regularized coordinates, using the
 *                      Heggie/Mikkola prescription.
 *                      Litt.: D.C. Heggie, Celestial Mechanics 10, 217 (1974);
 *                          S. Mikkola, Mon. Not. R. astr. Soc. 215, 171 (1985)
 *                      accepts: nonreg_sys:  a pointer to a newton0 system in
 *                                            4D form, but not regularized
 *                                            (4th components all vanishing);
 *                      returns: reg_sys: a pointer to a newly created 
 *                                        regularized system.
 *                      effect: initializes the regularized system as a 
 *                              transformation from the newton0 nonregularized
 *                              system, including memory allocation for the
 *                              mass matrix.
 *-----------------------------------------------------------------------------
 */
systptr  mk_reg_version(nonreg_sys)
systptr  nonreg_sys;
    {
    int  npart;
    int  nregpart;           /* number of regularized particles, i.e. number */
                             /* of particle pairs: npart(npart-1)/2          */
    systptr  reg_sys;

    npart = Nbody(nonreg_sys);
    nregpart = npart * (npart - 1) / 2;
    reg_sys = mk_system(nregpart);

    Massmatrix(reg_sys) = mk_massmatrix(npart);
    
    nonreg_to_reg(nonreg_sys, reg_sys);

    return(reg_sys);
    }

/*-----------------------------------------------------------------------------
 *  nonreg_to_reg  --  nonregularized  -->  regularized representation:
 *                     translates from the standard newton0 3-dimensional 
 *                     representation, embedded in 4 dimensions by having the
 *                     fourth component identically zero in all vectors, to a
 *                     new presentation in regularized coordinates, using the
 *                     Heggie/Mikkola prescription.
 *                     Litt.: D.C. Heggie, Celestial Mechanics 10, 217 (1974);
 *                          S. Mikkola, Mon. Not. R. astr. Soc. 215, 171 (1985)
 *                     accepts: nonreg_sys:  a pointer to a newton0 system in
 *                                           4D form, but not regularized (i.e.
 *                                           fourth components all vanishing);
 *                              reg_sys:  a pointer to the required amount
 *                                        and type of memory for the 
 *                                        corresponding regularized system.
 *                     effect: initializes the regularized system as a 
 *                             transformation from the newton0 nonregularized
 *                             system
 *                     NOTE: the memory for the regularized system, including
 *                           the mass matrix it points to, is assumed to be
 *                           provided; the easiest way to create this memory
 *                           is through an initial call to  mk_reg_version() ;
 *                           subsequent transformations can then be performed
 *                           through nonreg_to_reg()  and  reg_to_nonreg() .
 *                     note: com_to_noncom() transforms from to the original,
 *                           not necessarily center-of-mass, coordinate system.
 *                           combodies are created to preserve the old bodies;
 *                           at the end of the procedure  combodies  is freed.
 *-----------------------------------------------------------------------------
 */
void  nonreg_to_reg(nonreg_sys, reg_sys)
systptr  nonreg_sys;
systptr  reg_sys;
    {
    int  npart;
    real  npartinv;          /* 1/npart */
    real  q_k[NDIM];
    real  p_i[NDIM], p_j[NDIM], p_k[NDIM];      /* momentum vectors */
    real  two_p_k[NDIM];
    real  q_k_hat[NDIM][NDIM];       /* KS-matrix; see eq. 8 in Mikkola 1985 */
    real  q_k_hat_t[NDIM][NDIM];     /* matrix transpose of q_k_hat          */
    real  u_1;                       /* page 175, Mikkola 1985               */
    real  q_length;                  /*  | q_k |                             */
    real  *matrixelement;            /* element of the massmatrix            */
    bodyptr  bodies, combodies, regbodies;
    bodyptr  body_i, body_j, reg_k;
    bodyptr  body_i_prime, body_j_prime;

    Tnow(reg_sys) = Tnow(nonreg_sys);    /* real, unregularized time */

    npart = Nbody(nonreg_sys);
    npartinv = 1.0 / ((double) npart);

    bodies = Bodies(nonreg_sys);
    regbodies = Bodies(reg_sys);
/*
 * first separate the nonreg_sys center-of-mass motion from the internal
 * motion:
 */
    combodies = cp_bodies(bodies, npart);
    noncom_to_com(combodies, npart, Com_Pos(reg_sys),Com_Vel(reg_sys));

    reg_k = regbodies;
    for (body_i = combodies; body_i - combodies < npart-1; body_i++)
	for (body_j = body_i + 1; body_j - combodies < npart; body_j++)
	    {
	    SUBV(q_k, Pos(body_i), Pos(body_j));
            MULVS(p_i, Vel(body_i), Mass(body_i));
            MULVS(p_j, Vel(body_j), Mass(body_j));
	    SUBV(p_k, p_i, p_j);
	    INCMULVS(p_k, npartinv);

            PMass(reg_k) = Mass(body_i) * Mass(body_j);
	    
	    ABSV(q_length, q_k);
            u_1 = sqrt((q_length + ABS(q_k[0])) / 2.0);
	    if(q_k[0] >= 0.0)
		{
		PPos(reg_k)[0] = u_1;
		PPos(reg_k)[1] = q_k[1] / (2.0 * u_1);
		PPos(reg_k)[2] = q_k[2] / (2.0 * u_1);
	        PPos(reg_k)[3] = 0.0;
		}
	    else 
		{
		PPos(reg_k)[0] = q_k[1] / (2.0 * u_1);
		PPos(reg_k)[1] = u_1;
		PPos(reg_k)[2] = 0.0;
	        PPos(reg_k)[3] = q_k[2] / (2.0 * u_1);
		}

	    init_ks_matrix(q_k_hat, PPos(reg_k));
	    TRANM(q_k_hat_t, q_k_hat);
	    MULVS(two_p_k, p_k, 2.0);
	    MULMV(PMom(reg_k), q_k_hat_t, two_p_k);

	    reg_k++;
	    }

    matrixelement = Massmatrix(reg_sys);             /* first matrix element */
    for (body_i = combodies; body_i - combodies < npart-1; body_i++)
	for (body_j = body_i + 1; body_j - combodies < npart; body_j++)
            for (body_i_prime = combodies;
                           body_i_prime - combodies < npart-1; body_i_prime++)
	        for (body_j_prime = body_i_prime + 1;
                             body_j_prime - combodies < npart; body_j_prime++)
                    {
                    *matrixelement = 0.0;
		    if (body_i == body_i_prime)
			*matrixelement += 0.5 / Mass(body_i);
		    if (body_j == body_j_prime)
			*matrixelement += 0.5 / Mass(body_j);
		    if (body_i == body_j_prime)
			*matrixelement -= 0.5 / Mass(body_i);
		    if (body_j == body_i_prime)
			*matrixelement -= 0.5 / Mass(body_j);

                    matrixelement++;                  /* next matrix element */
		    }

    free(combodies);
    }

/*-----------------------------------------------------------------------------
 *  reg_to_nonreg  --  regularized  -->  nonregularized representation:
 *                     translates from a regularized representation back to the
 *                     the standard newton0 3-dimensional representation,
 *                     embedded in 4 dimensions, using the Heggie/Mikkola
 *                     prescription. In the absence of numerical
 *                     errors, the 4th components of positions and velocities
 *                     should all be identically zero. The extent to which they
 *                     are non-zero gives some measure of the numerical
 *                     accuracy.
 *                     Litt.: D.C. Heggie, Celestial Mechanics 10, 217 (1974);
 *                          S. Mikkola, Mon. Not. R. astr. Soc. 215, 171 (1985)
 *                     accepts: reg_sys:  a pointer to a regularized system
 *                                        in newton0 form;
 *                              nonreg_sys:  a pointer to a newton0 system in
 *                                           4D form, but not regularized.
 *                     effect: updates the nonregularized system as a 
 *                             transformation from the regularized system.
 *                     NOTE: the memory for the nonregularized system is 
 *                           assumed to be provided, and the masses
 *                           in nonreg_sys are assumed to still have their
 *                           correct values, left their from the time of an
 *                           earlier call either explicitly to  reg_to_nonreg()
 *                           or implicitly through  mk_reg_version() .
 *                           Having the old masses available is essential,
 *                           since they apear in Mikkola's equations on page
 *                           176, and are not available in the representation 
 *                           of the regularized system.
 *                     note: noncom_to_com() transforms back to the original,
 *                           not necessarily center-of-mass, coordinate system.
 *-----------------------------------------------------------------------------
 */
void  reg_to_nonreg(reg_sys, nonreg_sys)   /* systemptrs !!! no npart !!! */
systptr  reg_sys;
systptr  nonreg_sys;
    {
    int  npart;
    real  mtot;                      /* total mass */
    real  delta_r[NDIM];
    real  q_k[NDIM];
    real  p_k[NDIM];                 /* momentum vectors */
    real  q_k_hat[NDIM][NDIM];       /* KS-matrix; see eq. 8 in Mikkola 1985 */
    real  momentum_scalefactor;
    bodyptr  bodies, combodies, regbodies;
    bodyptr  body_i, body_j, reg_k;

    Tnow(nonreg_sys) = Tnow(reg_sys);    /* real, unregularized time */

    npart = Nbody(nonreg_sys);
    bodies = Bodies(nonreg_sys);
    regbodies = Bodies(reg_sys);

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	CLRV(Pos(body_i));
	CLRV(Vel(body_i));
	}

    reg_k = regbodies;
    for (body_i = bodies; body_i - bodies < npart-1; body_i++)
	for (body_j = body_i + 1; body_j - bodies < npart; body_j++)
	    {
	    init_ks_matrix(q_k_hat, PPos(reg_k));
            MULMV(q_k, q_k_hat, PPos(reg_k));
	    DOTVP(momentum_scalefactor, PPos(reg_k), PPos(reg_k));
            momentum_scalefactor = 1.0 / (2.0 * momentum_scalefactor);
            MULMV(p_k, q_k_hat, PMom(reg_k));
	    INCMULVS(p_k, momentum_scalefactor);

            MULVS(delta_r, q_k, Mass(body_j));
	    INCADDV(Pos(body_i), delta_r);
            MULVS(delta_r, q_k, Mass(body_i));
	    INCSUBV(Pos(body_j), delta_r);

	    INCADDV(Vel(body_i), p_k);
	    INCSUBV(Vel(body_j), p_k);

	    reg_k++;
	    }

    mtot = 0.0;
    for (body_i = bodies; body_i - bodies < npart; body_i++)
            mtot += Mass(body_i);

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	DIVVS(Pos(body_i), Pos(body_i), mtot);
	DIVVS(Vel(body_i), Vel(body_i), Mass(body_i));
	}

    com_to_noncom(bodies, npart, Com_Pos(reg_sys),Com_Vel(reg_sys));
    }


/*-----------------------------------------------------------------------------
 *  init_ks_matrix  --  initializes a Kustaanheimo-Stiefel matrix from a vector
 *                      accepts: ks_matrix: Kustaanheimo-Stiefel matrix;
 *                                     vec: a vector which may either be a
 *                                          position or a momentum vector;
 *-----------------------------------------------------------------------------
 */
void  init_ks_matrix(ks_matrix, vec)
real  ks_matrix[NDIM][NDIM];
real  vec[NDIM];
    {
    ks_matrix[0][0] = vec[0];
    ks_matrix[0][1] = -vec[1];
    ks_matrix[0][2] = -vec[2];
    ks_matrix[0][3] = vec[3];

    ks_matrix[1][0] = vec[1];
    ks_matrix[1][1] = vec[0];
    ks_matrix[1][2] = -vec[3];
    ks_matrix[1][3] = -vec[2];

    ks_matrix[2][0] = vec[2];
    ks_matrix[2][1] = vec[3];
    ks_matrix[2][2] = vec[0];
    ks_matrix[2][3] = vec[1];

    ks_matrix[3][0] = vec[3];
    ks_matrix[3][1] = -vec[2];
    ks_matrix[3][2] = vec[1];
    ks_matrix[3][3] = -vec[0];
    }

/*-----------------------------------------------------------------------------
 *  mk_massmatrix  --  allocate memory for the regularized mass matrix which
 *                     contains the coefficients of the products of velocity 
 *                     components in the Hamiltonian.
 *                     Litt.: eq. (6) in Mikkola (1985).
 *-----------------------------------------------------------------------------
 */
realptr  mk_massmatrix(npart)
int  npart;
    {
    int  nregpart;           /* number of regularized particles, i.e. number */
                             /* of particle pairs: npart(npart-1)/2          */
    realptr  massmatrix;

    nregpart = npart * (npart - 1) / 2;

    massmatrix = (realptr) malloc((unsigned)nregpart * nregpart * sizeof(real));
    if (massmatrix == NULL)
	error("mk_massmatrix: not enough memory left\n");

    return(massmatrix);
    }

/*-----------------------------------------------------------------------------
 *  noncom_to_com  --  transforms  bodies  to the center-of-mass coordinate
 *                     system. The position and velocity of the center of mass
 *                     motion are stored in  com_pos  and  com_vel .
 *-----------------------------------------------------------------------------
 */
local void  noncom_to_com(bodies, npart, com_pos, com_vel)
bodyptr  bodies;
int  npart;
real  com_pos[];
real  com_vel[];
    {
    bodyptr  body_i;
    real  total_mass;
    real  mpi[NDIM], mvi[NDIM];

    total_mass = 0.0;
    CLRV(com_pos);
    CLRV(com_vel);

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	total_mass += Mass(body_i);
	MULVS(mpi, Pos(body_i), Mass(body_i));
	INCADDV(com_pos, mpi);                        /* here: mass-weighted */
	MULVS(mvi, Vel(body_i), Mass(body_i));
	INCADDV(com_vel, mvi);                        /* here: mass-weighted */
	}

    INCDIVVS(com_pos, total_mass);         /* here: the real c.o.m. position */
    INCDIVVS(com_vel, total_mass);         /* here: the real c.o.m. velocity */

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	INCSUBV(Pos(body_i), com_pos);
	INCSUBV(Vel(body_i), com_vel);
	}
    }

/*-----------------------------------------------------------------------------
 *  com_to_noncom  --  transforms  bodies  from the center-of-mass coordinate
 *                     system. The position and velocity of the center of mass
 *                     motion are stored in  com_pos  and  com_vel .
 *-----------------------------------------------------------------------------
 */
local void  com_to_noncom(bodies, npart, com_pos, com_vel)
bodyptr  bodies;
int  npart;
real  com_pos[];
real  com_vel[];
    {
    bodyptr  body_i;

    for (body_i = bodies; body_i - bodies < npart; body_i++)
	{
	INCADDV(Pos(body_i), com_pos);
	INCADDV(Vel(body_i), com_vel);
	}
    }

/* endof: transformreg.c */
