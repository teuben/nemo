/*
 * QUADFIELD.H: definitions for quadrupole force tables.  Values for
 * internal and external moments Qnk and Pnk are tabulated for radii
 * running from 0 to the radius of the outermost body plus epsilon.
 * with MQTAB-2 values stored at the radii of actual bodies.
 */

#define MQTAB  (2 + 64)			/* maximum entries in table */

typedef struct {
    int    nqtab;				/* entries in tables */
    real   radtab[MQTAB];			/* tabulated radii */
    real   Q00tab[MQTAB];			/* tabulated moments */
    real   P00tab[MQTAB];
    real   Q10tab[MQTAB];
    real   P10tab[MQTAB];
    vector Q11tab[MQTAB];
    vector P11tab[MQTAB];
    matrix Q22tab[MQTAB];
    matrix P22tab[MQTAB];
} quadfield;
