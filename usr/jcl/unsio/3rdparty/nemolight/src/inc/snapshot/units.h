/*
 * UNITS.H: definition file for various systems of units.
 */

/*
 * Names for different systems of units in which to represent galaxy models.
 * All systems are defined by requiring that the gravitational constant G = 1,
 * and that the total mass M = 1. The remaining freedom lies in the length
 * scale, for which the natural choice is less obvious than the above two.
 * In fact, for different applications, different systems of units are
 * convenient. We list here all systems presently defined, and their
 * corresponding definition of length scale.
 */

#define  VIRIAL  001     /*      total energy E = -1/4 . 
                          *  In these units the virial radius r_v = 1 ,
			  *  where r_v is defined as the harmonic mean particle
			  *  separation (Heggie and Mathieu).
                          *  With this choice, a system of equal-mass
			  *  particles has a root-mean-square particle
			  *                2  1/2     -1/2
			  *  velocity   < v  >    =  2     .
                          *  In terms of the half-mass radius r_h, typically
                          *  r_h = 0.8 r_v , with the coefficient being
			  *  generally in the range 0.76 - 0.98 (Spitzer).
                          *  Litt.: Heggie, D.C. and Mathieu, R.D. 1986, in
			  *           "The Use of Supercomputers in Stellar
			  *           Dynamics", eds. P. Hut and S. McMillan
			  *           (Springer, Lecture Notes in Physics 267),
			  *           p. 233.
			  *         Spitzer, L. 1987, Dynamical Evolution of 
			  *           globular clusters (Princeton Univ.
			  *           Press), Sect. 1.2a .
                          *  Example: for a Plummer model,
			  *                      ___________
			  *             16      /  2/3      |
			  *      r  = ------ \ /  2    - 1    r  = 1.3011 r
			  *       v    3.pi   V                h           h
			  *
                          *  Example: for a constant-density sphere,
			  *
			  *              5
			  *      r  = -------- r  = 1.0499 r
			  *       v       2/3   h           h
			  *            3.2
			  */

#define  RMSVEL  002     /*      total energy E = -1/2 . 
                          *  With this choice, a system of equal-mass
			  *  particles has a root-mean-square particle
			  *                2  1/2
			  *  velocity   < v  >    = 1 .
			  *  The unit length in this system is r_{rv} = 2 r_v ,
			  *  where r_v is defined as the harmonic mean particle
			  *  separation (Heggie and Mathieu). In terms of the
                          *  half-mass radius r_h, typically  r_h = 0.4 r_{rv},
			  *  with the coefficient being generally in the range
			  *  0.38 - 0.49 (Spitzer).
                          *  Litt.: Heggie, D.C. and Mathieu, R.D. 1986, in
			  *           "The Use of Supercomputers in Stellar
			  *           Dynamics", eds. P. Hut and S. McMillan
			  *           (Springer, Lecture Notes in Physics 267),
			  *           p. 233.
                          *  Example: for a Plummer model,
			  *                      ___________
			  *             32      /  2/3      |
			  *     r   = ------ \ /  2    - 1    r  = 2.6022 r
			  *      rv    3.pi   V                h           h
			  *
                          *  Example: for a constant-density sphere,
			  * 
			  *               1/3
			  *            5.2
			  *     r   = ------  r  = 2.0999 r
			  *      rv     3      h           h
			  */

#define  STRUCTURAL  003 /*      structural length r_0 = 1 .
			  *  The structural length is a unit of length chosen
			  *  so as to simplify one of the basic expression,
			  *  e.g. the density (in a Plummer model) or Poisson's
			  *  equation (for a King model). It is therefore not
			  *  a uniquely defined expression, and care has to be
			  *  taken to check the precise definition for each
			  *  model separately.
			  *  Litt.: Binney, J.J. and Tremaine, S.D. 1987,
                          *         Galactic Astronomy, Vol. 2: Galactic
                          *         Dynamics (Princeton Univ. Press), Ch. 4 .
                          *  Example: for a Plummer model,
			  *                ___________
			  *               /  2/3      |
			  *      r  =  \ /  2    - 1    r  = 0.7664 r
			  *       0     V                h           h
			  */

#define  HALFMASS   004  /*      half-mass radius r_h = 1 .
			  *  This radius contains half the mass.
                          */

#define HALFPROJMASS 005 /*      half-mass radius r_{hP} = 1 .
			  *  This radius includes half the integrated
			  *  projected surface mass density.
			  *  Litt.: Spitzer, L. 1987, Dynamical Evolution of 
			  *           globular clusters (Princeton Univ.
			  *           Press), Sect. 1.2b .
                          *  Example: for a Plummer model,
			  *                   ___________
			  *                  /  2/3      |
			  *    r   = r  = \ /  2    - 1    r  = 0.7664 r
			  *     hP    0    V                h           h
                          *
                          *  Example: for a constant-density sphere (same!),
			  *                   ___________
			  *                  /  2/3      |
			  *    r   = r  = \ /  2    - 1    r  = 0.7664 r
			  *     hP    0    V                h           h
                          */

#define  TOTALMASS  006  /*      total-mass radius r_t = 1 .
			  *  This is the radius of the whole system, the point
			  *  at which the density drops to zero. 
                          *  If the density is non-zero out to infinity, this
			  *  system of units does not apply.
                          *  Example: for a constant-density sphere,
			  *           
			  *                    1/3
			  *            r   =  2    r  = 1.2599 r
			  *             t           h           h
                          */

#define  PROJCORE   007  /*      core-mass radius r_c = 1 .
                          *  The core radius is defined as the radius at which
			  *  the projected surface density drops off to half
			  *  the central value (King). This observationally
			  *  motivated definition applies only for models with
			  *  finite central surface density. For example, a
			  *  model with a central density similar to a singular
			  *  isothermal sphere will have r_c = 0 .
                          *  Note, however, that a divergent density can still
			  *  lead to a well-defined finite core radius, as in
			  *  the case of a de Vaucouleurs model.
                          *  Example: for a Plummer model,
			  *              ___________________
			  *             /  2/3       1/2    |
			  *      r = \ / (2   - 1) (2   - 1)  r  = 0.4933 r
			  *       c   V                        h           h
			  *           
                          *  Example: for a constant-density sphere,
			  *                 
			  *                  -2/3  1/2
			  *           r  =  2     3    r  =  1.0911 r
			  *            c                h            h
                          */

/*
 * A2UNIT: Integer-valued function to convert from ascii to unit number.
 */

int a2unit();

/*
 * _UNIT_TABLE: body of table with ascii names for units, used by a2unit().
 * NOTE: the numbering of units must correspond to the order in this table!
 */

#define _UNIT_TABLE {   \
    "VIRIAL",		\
    "RMSVEL",		\
    "STRUCTURAL",	\
    "HALFMASS",		\
    "HALFPROJMASS",	\
    "TOTALMASS",	\
    "PROJCORE",		\
    NULL,		\
}
