/* -*- C -*-
 *******************************************************************************
 *                                                                             *
 * falcON_C.h                                                                  *
 *                                                                             *
 * C header file                                                               *
 *                                                                             *
 * Copyright Walter Dehnen, 2000-2004                                          *
 * e-mail:   walter.dehnen@astro.le.ac.uk                                      *
 * address:  Department of Physics and Astronomy, University of Leicester      *
 *           University Road, Leicester LE1 7RH, United Kingdom                *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 * falcON = Force ALgorithm with Complexity O(N)                               *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 * header file for C users.                                                    *
 * (C++ and FORTRAN users, see files falcON.h and falcON.f, respectively)      *
 *                                                                             *
 * C routines  implementing the code described by Dehnen (2000,2002). These    *
 * routines call the original C++ functions, declared in file falcON.h.        *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 * See file readme.falcON.html for some guidelines on how to install, #include *
 * and link this file and the code declared in it. You have the choice between *
 * various options determining some properties of the code.                    *
 *                                                                             *
 *******************************************************************************
 */
#ifndef falcON_included_falcONC_h             /* ensure this file is seen     */
#define falcON_included_falcONC_h 1           /* once only by the compiler    */

#ifdef __cplusplus                            /* C++ users should better use  */
extern "C" {                                  /* falcON.h.                    */
#ifndef falcONC_cc
#warning "falcON_C.h is for C. with C++ better use falcON.h"
#endif
#endif                                        /* this is for savety only      */

#if defined(falcON_SINGLE_DOUBLE) || defined(falcON_DOUBLE_DOUBLE)
#define INPUT_TYPE double                     /*   define input type = double */
#else                                         /*                              */
#define INPUT_TYPE float                      /*   define input type = float  */
#endif

#ifndef falcON_NDIM
#  define falcON_NDIM 3
#else
#  define falcON_NDIM 2
#endif
/*******************************************************************************
 *                                                                             *
 *  CONTENTS                                                                   *
 *  ========                                                                   *
 *                                                                             *
 *  1  Initialisation and clearing up                                          *
 *     1.1  Before using the code                                              *
 *     1.2  After using the code                                               *
 *  2  Meaning of the body flags                                               *
 *  3  Force approximation and related                                         *
 *     3.1  Generating a tree structure                                        *
 *     3.2  Approximating the accelerations                                    *
 *     3.3  A crude estimation of the mass- and number density                 *
 *  4  Search for and counting of neighbours and collision partners            *
 *     4.1  Neighbour search (SPH support)                                     *
 *     4.2  Collision partner search (for sticky particles)                    *
 *     4.3  SPH Neighbour counting                                             *
 *  5  Other features                                                          *
 *  6  Known bugs and problems                                                 *
 *     6.1  Test-bodies are not possible                                       *
 *     6.2  Bodies at identical positions                                      *
 *  References                                                                 *
 *  A  Routines for FORTRAN support only                                       *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 *  1 INITIALISATION AND CLEARING UP                                           *
 *  ================================                                           *
 *                                                                             *
 *  1.1 BEFORE USING THE CODE                                                  *
 *  -------------------------                                                  *
 *  In order to use the code, one has first to initialize it. This is done     *
 *  using the following routine, which actually does no computation at all.    *
 */
void  falcON_initialize(const int*,               /* flags                    */
			const INPUT_TYPE*,        /* masses                   */
			const INPUT_TYPE**,       /* positions                */
#ifdef falcON_INDI
			const INPUT_TYPE*,        /* eps                      */
#endif
			      INPUT_TYPE**,       /* accelerations            */
		              INPUT_TYPE*,        /* potentials               */
		              INPUT_TYPE*,        /* densities                */
			const int,                /* N   = # bodies           */
			const INPUT_TYPE,         /* eps = softening length   */
			const INPUT_TYPE,         /* theta= opening angle     */
			const int,                /* softening kernel         */
			const INPUT_TYPE);        /* Newton's constant G      */
/*                                                                             *
 * The first 11[10] arguments specify the sink and source properties of the    *
 * bodies. Each body has the sink properties: position (x,y,z), mass,          *
 * softening length (for the case of individual softening lengths), and flag   *
 * (see section 2 below). The source properties are: acceleration(ax,ay,az),   *
 * potential and mass-density. If individual softening is enabled (via the     *
 * preprocessor flag "falcON_INDI", see file make.defs) but global softening   *
 * is used, a NULL pointer may be used instead of an array holding eps_i. The  *
 * same applies to the arrays for density and potential: if a NULL pointer is  *
 * given, they will never be updated (but possibly computed).                  *
 *                                                                             *
 * The last 4 arguments are:                                                   *
 *                                                                             *
 * N:                  size of arrays == number of bodies.                     *
 * EPS:   global overall softening length, unless eps_i are given.             *
 * THETA: if(THETA>0): theta = theta(M) with theta_min = THETA (Dehnen 2002)   *
 *        if(THETA<0): theta = THETA = const                                   *
 *        RECOMMENDED: THETA = 0.5 - 0.6                                       *
 * KERN:  0,1,2,3:     P_n kernels (P_0=Plummer) (Dehnen & Teuben, 2002)       *
 *                                                                             *
 * With the functions                                                          *
 */
void falcON_resetsoftening(const INPUT_TYPE, /* fixed/maximum softening length*/
			   const int);       /* softening kernel              */
void falcON_resetopening  (const INPUT_TYPE);/* opening angle                 */
/*
 * it is possible to change, after initialisation, the softening length and    *
 * kernel as well as the opening criterion.                                    *
 *                                                                             *
 * If falcON_initialize() is called more than once, only the last              *
 * initialisation applies, the older ones will be deleted. That is, you can    *
 * have only one tree at a time -- use the C++ version if you need more.       *
 *                                                                             *
 * 1.2 AFTER USING THE CODE                                                    *
 * ------------------------                                                    *
 * After using the code, one should delete allocated memory by calling         *
 */
void falcON_clearup();                       /*                               */
/*
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * 2 MEANING OF THE BODY FLAGS                                                 *
 * ===========================                                                 *
 *                                                                             *
 * The body flags are integers with the following meaning:                     *
 *                                                                             *
 * bit  value     meaning                                                      *
 * --------------------------------------------------------------              *
 *   1      1     this body is active, i.e. wants update                       *
 *   2      2     don't load this body into the tree, ignore it                *
 *   3      4     this body is a SPH particle                                  *
 *   4      8     this body is a sticky particle                               *
 * i>4   2^(i-1)  not used                                                     *
 *                                                                             *
 *                                                                             *
 * The default, flag=0, represents a plain body that is inactive, but still    *
 * source of gravity.                                                          *
 * The flag is obtained by setting the bits or, equivalently, adding the       *
 * values. The flags are not used by falcON_initialize() and will unfold their *
 * effect only when The functions in section 3 and 4 below are called.         *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * 3 FORCE APPROXIMATION AND RELATED                                           *
 * =================================                                           *
 *                                                                             *
 * 3.1 Generating a tree structure                                             *
 * -------------------------------                                             *
 * In order to establish a hierarchical tree structure from all bodies that    *
 * are not flagged to be ignored, and subsequently pre-compute the centers     *
 * of mass and sizes of the cells, use                                         *
 */
void falcON_grow(                            /*  grow new tree structure      */
		 const int);                 /*  Ncrit: critical cell size    */
/*
 * which creates a tree from scratch. Cells containing Ncrit or less bodies    *
 * will NOT be splitted. Experiments showed that for a full force calculation  *
 * Ncrit = 6-8 is about optimal and results in a few % decrease in CPU time    *
 * consumption (compared to Ncrit=1) and about 20% reduction in memory.        *
 *                                                                             *
 *                                                                             *
 * You may instead also re-use an old tree structure and only re-compute the   *
 * centers of mass and sizes of the cells (since the bodies have moved), by    *
 * using                                                                       *
 */
void falcON_reuse();
/*
 * In this case, the logical linkage between cells and bodies is preserved     *
 * and no new tree structure is established. However, the code accounts for    *
 * the fact that the bodies might have moved. If the bodies have moved a lot,  *
 * the sizes of the cells will be much larger than the physical size of the    *
 * associated boxes, and the tree traversal will be very inefficient.          *
 * However, if the bodies have moved only little, the force compuation is      *
 * hardly slowed down and re-using an old tree saves the costs for             *
 * establishing a new one. Note that reuse() does not allow for any change in  *
 * the flags indicating whether or not a body shall be ignored: changing this  *
 * flag will have an effect only at the next call of grow().                   *
 *                    +----------------------------+                           *
 *                    | USE THIS OPTION CAREFULLY! |                           *
 *                    +----------------------------+                           *
 * You have been warned!                                                       *
 *                                                                             *
 *                                                                             *
 * 3.2 Approximating the Accelerations                                         *
 * -----------------------------------                                         *
 * Once a tree structure is established, you can compute the forces by         *
 */
void falcON_approx_grav();
/*
 * after a call to falcON_grow()  or falcON_reuse().                           *
 *                                                                             *
 * falcON_approx_grav() implements the pre-computation of the quadrupoles, as  *
 * well as the interaction and evaluation phase. See src/exe/C/TestGravC.c     *
 * for an example application.                                                 *
 */
#ifdef falcON_ADAP
/* For individual adaptive softening the routine                               *
 *                                                                            */
  void falcON_adjust_epsi_and_approx_grav(INPUT_TYPE,  /* I: Nsoft            */
					  int,         /* I: Nref             */
					  INPUT_TYPE); /* I: fac              */
/*
 * can do more for you before the forces are actually computed:                *
 *                                                                             *
 * If Nsoft [2nd argument] is non-zero, it estimates for each active particle  *
 * the local number density (using the number density of the smallest cell     *
 * containing that particle with not less than Nref [3rd argument] bodies).    *
 * If fac [4th arg] is zero, the bodies softening lengths are set such that,   *
 * based on the estimated number density, their eps-spheres contain Nsoft      *
 * bodies, but eps <= EPS [global parameter].                                  *
 * If fac [4th arg] is non zero, Eps is computed in the same way, and the      *
 * new softening is set to                                                     *
 *        eps_new = Eps^2 / eps_old,                                           *
 * with the restriction  eps_new in [eps_old/fac, eps_old*fac]                 *
 *                                                                             *
 * If eps_i is adjusted in this way, it will be copied back to the bodies.     *
 */
#endif
/*
 *                                                                             *
 * 3.3 A crude Estimation of the Mass- and Number-Density                      *
 * ------------------------------------------------------                      *
 * There is also the possibility to obtain a rough estimate of the mass- or    *
 * number density of bodies in the neighbourhood of every body flagged being   *
 * active, via                                                                 *
 */
void falcON_estimate_rho  (                  /* estimate mass density         */
			   const int);       /* I:  critical cell size        */
void falcON_estimate_n    (                  /* estimate number density       */
			   const int);       /* I:  critical cell size        */
/*
 * These estimates are simply the mean density within the smallest cells with  *
 * more then NX (1st arg) bodies. Note that when using test bodies (bodies     *
 * with zero or tiny mass), this guess for the mass density can have terrible  *
 * errors. Moreover, when the tree has not been grow()n but simply reuse()ed,  *
 * these estimate will not change. Be careful using these functions.           *
 * YOU HAVE BEEN WARNED.                                                       *
 *                                                                             *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * 4 SEARCH FOR AND COUNTING OF COLLISION PARTNERS                             *
 * ===============================================                             *
 *                                                                             *
 * Once a tree structure is established (i.e. after a call to falcON_grow()    *
 * or falcON_reuse()), you can also use it to create interaction lists for SPH *
 * and sticky particles via the routine                                        *
 */
void falcON_iactionlist(      int*,          /* list of indices: 1st of pair  */
		              int*,          /* list of indices: 2nd of pair  */
			const int,           /* physical size of list         */
			      int*,          /* actual size of list           */
			const INPUT_TYPE*,   /* array with body sizes         */
			const bool,          /* use Max(h_i,h_j) OR h_i+h_j ? */
			const INPUT_TYPE,    /* time step tau                 */
			const INPUT_TYPE**); /* arrays with V                 */
/*
 * In case of overflow, i.e. if the number of pairs found exceeds the size     *
 * (3rd arg) of the list (1st&2nd args), a warning is issued to stderr, but    *
 * the search is not truncated, rather pairs are no longer copied into the     *
 * interaction list. In this case, the value returned for the actual size      *
 * (4th arg) exceeds the maximum size (3rd arg), but reflects the actual       *
 * number of interactions found.                                               *
 * It is the responsibility of the user to ensure that the arrays for sizes    *
 * and velocity components of the bodies are properly allocated (it is         *
 * sufficient to have entries for all bodies flagged as sticky/SPH.            *
 *                                                                             *
 *                                                                             *
 * 4.1 SPH support: neighbour or interaction partner search                    *
 * --------------------------------------------------------                    *
 * In order to make a list of all pairs {i,j} of indices for which             *
 *                                                                             *
 *      (1) both flags indicate SPH particles,                                 *
 * and  (2) at least one flagged being active,                                 *
 * and  (3)     | x_i - x_j | < max(size_i,size_j)     IF Max==true            *
 *          OR  | x_i - x_j | < size_i + size_j        IF Max==false           *
 *                                                                             *
 * use falcON_iactionlist() with negative time step (7th arg) and provide an   *
 * array with body sizes (5th arg) but give NULL pointers for the arrays with  *
 * velocity components (last args). Note that the 6th argument determines      *
 * whether you search for neihbours or interaction partners.                   *
 *                                                                             *
 *                                                                             *
 * 4.1 Sticky-particle support                                                 *
 * ---------------------------                                                 *
 * In order to make a list of all pairs {i,j} of indices for which             *
 *                                                                             *
 *      (1) both flags indicate sticky particles,                              *
 * and  (2) at least one is flagged being active,                              *
 * and  (3) | (x_i+t*v_i)-(x_j+t*v_j) | < size_i+size_j  with t in [0,tau],    *
 *                                                                             *
 * use flacON_iactionlist() with the time step (6th arg) >= 0 and provide      *
 * arrays with sizes (5th arg) and velocity components (last args) of the      *
 * bodies.                                                                     *
 *                                                                             *
 *                                                                             *
 * See file src/C/TestPairC.c for an example application for both supports.    *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * 5 OTHER FEATURES                                                            *
 * ================                                                            *
 *                                                                             *
 * There are many further routines by which the user can obtain information    *
 * about the behaviour of the code.                                            *
 */
INPUT_TYPE falcON_root_center(const int i);  /* ith (0,1,2) comp. of          */
					     /*	center of root                */

INPUT_TYPE falcON_root_radius();             /* radius of root                */

INPUT_TYPE falcON_current_eps();             /* softening length              */

int  falcON_current_kernel();                /* kernel: 0/1/2/3/              */

int  falcON_softening();                     /* softening:                    */
 					     /* 0/1 == fixed/individual       */

int  falcON_No_cells();                      /* number of cells used          */

int  falcON_depth();                         /* max depth of tree             */

void falcON_stats();                         /* writes some statistics        */
					     /*	to standard output            */
/*
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * 6 KNOWN BUGS AND PROBLEMS                                                   *
 * =========================                                                   *
 *                                                                             *
 * 6.1 Test-bodies are not possible                                            *
 * --------------------------------                                            *
 * A body that is loaded into the tree but has zero mass, will not acquire     *
 * any acceleration. This is because the code computes first the force = mass  *
 * times acceleration (it is symmetric and hence better suited for             *
 * computation of mutual interactions) and then divides by the mass to obtain  *
 * the acceleration. The only possible work-around this problem is to set      *
 * the mass of potential test bodies much smaller than the masses of source    *
 * bodies. However, this must be done such that the gravity of test bodies     *
 * is everywhere neglible compared to that of the source bodies.               *
 * Note, however, that this work-around is wasteful: it computes the forces    *
 * generated by the test bodies. (You have been warned!)                       *
 *                                                                             *
 * 6.2 Bodies at identical positions                                           *
 * ---------------------------------                                           *
 * The code cannot cope with more than Ncrit bodies at the same position       *
 * (within the floating point accuracy). This situation will lead to an        *
 * infinitely deep tree, i.e. the maximum allowed tree depth will be exceeded  *
 * and the code aborts with an error message.                                  *
 *                                                                             *
 *******************************************************************************
 *                                                                             *
 *                                                                             *
 * REFERENCES                                                                  *
 * ==========                                                                  *
 *                                                                             *
 * Dehnen W., 2000, ApJ, 536, L39                                              *
 * Dehnen W., 2001, MNRAS, 324, 273                                            *
 * Dehnen W., 2002, J. Comp. Phys, 179, 27                                     *
 * Dehnen W. & Teuben P.J., 2002, in preparation                               *
 *                                                                             *
 *******************************************************************************
 */
#undef INPUT_TYPE                             /* not to be used elsewhere     */
/******************************************************************************/
#ifdef __cplusplus
}
#endif
#endif                                        /* falcON_included_falcONC_h    */
