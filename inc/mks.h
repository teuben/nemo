/* 
 *  Some fundamental physical constants - mks units
 *
 *  See also:
 *  "The Fundamental Physical Constants" by
 *      E. Richard Cohen and Barry N. Taylor
 *          (PHYICS TODAY, August 1989, Austus 1994)
 *  https://physics.nist.gov/
 */


#define c_MKS 299792458.0000    /* speed of light [m/s] */
#define k_MKS 1.380649e-23      /* Boltzmann  [J/K] */
#define h_MKS 6.62607015e-34    /* Planck [J-s] */
#define G_MKS 6.67430e-11       /* Newton's Gravitational Constant "G" [N m^2/kg^2] 6.674 30 */
#define g_MKS 9.80665           /* standard acceleration of gravity */

#define M_SOLAR   1.989e30	/* Mass of Sun [kg] */
#define L_SOLAR   3.827e26 	/* Luminosity of Sun [W] */
#define R_SOLAR	  6.96e8	/* Radius of Sun [m] */
#define M_EARTH   5.974e24	/* Mass of earth [kg] */
/* #define AU	  1.49597870691e11     Astronomical Unit [m]  */
#define AU        1.49597870700e11  /* AU as per IAU 2012 definition */
#define PC        3.0856e16     /* Parsec [m] */

#define HI_MHz    1420.40575177 /* 21cm (21.10611405413 cm) line  of HI */

/*
 *  1au = exactly 149597870700 metres, in agreement with the value adopted in IAU 2009 Resolution B2
 *  1pc = 648000/pi AU (IAU 2015 Resolution B2) = 3.085677581491367e+16
 */
