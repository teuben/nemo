/* differentiatereg.c - BIGNUMBER, get_massmatrix_ij 
                        heggie_mikkola_equations_of_motion */

/*
 *  differentiatereg.c:  contains helper functions for regularized integration
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static real get_massmatrix_ij(realptr massmatrix, int nregpart, int i, int j);

#define  BIGNUMBER  1.0e20

/*-----------------------------------------------------------------------------
 *  heggie_mikkola_equations_of_motion  --  force calculation for regularized
 *                                          coordinates
 *  NOTE:
 *       this giant procedure works, but is not very elegant;
 *       future modularization plans are on the back burner.
 *-----------------------------------------------------------------------------
 */

void  heggie_mikkola_equations_of_motion(reg_sys, specs)
systptr  reg_sys;
specptr  specs;
    {
    int  nregpart;
    real  q_k[NDIM];
    real  q_k_hat[NDIM][NDIM];       /* KS-matrix; see eq. 8 in Mikkola 1985 */
    real  q_k_hat_t[NDIM][NDIM];     /* its transpose */
    real  p_k_hat[NDIM][NDIM];       /* KS-matrix; see eq. 8 in Mikkola 1985 */
    real  p_k_hat_t[NDIM][NDIM];     /* its transpose */
    real  scalefactor;
    real  negpot_energy;          /* minus potential energy */
    real  kin_energy;
    real  t_p_k[NDIM];
    real  t_q_k[NDIM];
    real  u_q_k[NDIM];
    real  vec1[NDIM], vec2[NDIM];    /* tmp storage for intermediate results */
    real  g_hamiltonian;
    real  g_t;
    real  g_u;
    realptr  d_k;
    realptr  the_d_k;
    real  d_a_k[NDIM];
    real  pair_min;                  /* minimum pair separation */
    cbodyptr  one_p_k;                 /* momentum vector */
    cbodyptr  all_p_k;
    cbodyptr  one_a_k;
    cbodyptr  all_a_k;
    cbodyptr  one_a_star_k;
    cbodyptr  all_a_star_k;
    bodyptr  regbodies;
    bodyptr  reg_k;

    regbodies = Bodies(reg_sys);
    nregpart = Nbody(reg_sys);
    all_p_k = mk_cbodies(nregpart);
    all_a_k = mk_cbodies(nregpart);
    all_a_star_k = mk_cbodies(nregpart);
    the_d_k = (realptr) malloc((unsigned)nregpart * sizeof(real));
    if (the_d_k == NULL)
	error("heggie_mikkola_equations_of_motion: no memory for d_k array\n");

    pair_min = BIGNUMBER;
    negpot_energy = 0.0;
    for (reg_k = regbodies, one_p_k = all_p_k; reg_k - regbodies < nregpart;
	 reg_k++, one_p_k++)
        {
        init_ks_matrix(q_k_hat, PPos(reg_k));
        DOTVP(scalefactor, PPos(reg_k), PPos(reg_k));
	pair_min = MIN(pair_min, scalefactor);
        negpot_energy += PMass(reg_k) / scalefactor;  /* potential interlude */
        scalefactor = 1.0 / (2.0 * scalefactor);
        MULMV(Config(one_p_k), q_k_hat, PMom(reg_k));
        INCMULVS(Config(one_p_k), scalefactor);
	}
    SPECmin_pair(specs) = MIN(SPECmin_pair(specs), pair_min);

    for (one_a_k = all_a_k; one_a_k - all_a_k < nregpart; one_a_k++)
        {
        CLRV(Config(one_a_k));
        for (one_p_k = all_p_k; one_p_k - all_p_k < nregpart; one_p_k++)
	    {
	    MULVS(d_a_k, Config(one_p_k),get_massmatrix_ij(Massmatrix(reg_sys),
                  nregpart, one_a_k - all_a_k, one_p_k - all_p_k));

	    INCADDV(Config(one_a_k), d_a_k);
	    }
	}

    for (one_a_k = all_a_k, one_a_star_k = all_a_star_k;
         one_a_k - all_a_k < nregpart; one_a_k++, one_a_star_k++)
        {
	Config(one_a_star_k)[0] =  Config(one_a_k)[0];
	Config(one_a_star_k)[1] =  Config(one_a_k)[1];
	Config(one_a_star_k)[2] =  Config(one_a_k)[2];
	Config(one_a_star_k)[3] = -Config(one_a_k)[3];
	}

    kin_energy = 0.0;
    for (d_k = the_d_k, one_a_k = all_a_k, one_p_k = all_p_k;
         d_k - the_d_k < nregpart; d_k++, one_a_k++, one_p_k++)
        {
        DOTVP(*d_k, Config(one_a_k), Config(one_p_k));
	kin_energy += *d_k;
        }

    for (reg_k = regbodies, one_p_k = all_p_k, d_k = the_d_k,
           one_a_k = all_a_k, one_a_star_k = all_a_star_k;
         reg_k - regbodies < nregpart;
	 reg_k++, one_p_k++, d_k++, one_a_k++, one_a_star_k++)
        {
        init_ks_matrix(q_k_hat, PPos(reg_k));
        TRANM(q_k_hat_t, q_k_hat);
        init_ks_matrix(p_k_hat, PMom(reg_k));
        TRANM(p_k_hat_t, p_k_hat);

        DOTVP(scalefactor, PPos(reg_k), PPos(reg_k));
        scalefactor = 1.0 / scalefactor;

	MULMV(t_p_k, q_k_hat_t, Config(one_a_k));
	INCMULVS(t_p_k, scalefactor);

	MULMV(vec1, p_k_hat_t, Config(one_a_star_k));
	MULVS(vec2, PPos(reg_k), *d_k);
	INCMULVS(vec2, 4.0);
        SUBV(t_q_k, vec1, vec2);
	INCMULVS(t_q_k, scalefactor);
	
        scalefactor = -2.0 * scalefactor * scalefactor;
	MULVS(u_q_k, PPos(reg_k), PMass(reg_k));
	INCMULVS(u_q_k, scalefactor);

        Reglagrangian(reg_sys) = kin_energy + negpot_energy;
        Reghamiltonian(reg_sys) = kin_energy - negpot_energy;
	g_hamiltonian = (Reghamiltonian(reg_sys) - Regenergy(reg_sys)) /
                                                        Reglagrangian(reg_sys);
	g_t = (1.0 - g_hamiltonian) / Reglagrangian(reg_sys);
	g_u = -(1.0 + g_hamiltonian) / Reglagrangian(reg_sys);

        MULVS(PdPos_ds(reg_k), t_p_k, g_t);
        MULVS(vec1, t_q_k, -g_t);
        MULVS(vec2, u_q_k, g_u);
        SUBV(PdMom_ds(reg_k), vec1, vec2);
	}

    free(the_d_k);
    free(all_p_k);
    free(all_a_k);
    free(all_a_star_k);
    }

/*-----------------------------------------------------------------------------
 *  get_massmatrix_ij  --  return the [i][j] element of the mass matrix.
 *                         This procedures hides the bookkeeping, caused by
 *                         fact that the mass matrix is really stored as a
 *                         one-dimensional array (since at compilation time
 *                         the size of the matrix is not yet known).
 *-----------------------------------------------------------------------------
 */
local real  get_massmatrix_ij(massmatrix, nregpart, i, j)
realptr  massmatrix;
int  nregpart;
int  i;
int  j;
    {
    return( massmatrix[i * nregpart + j] );
    }

/* endof: differentiatereg.c */
