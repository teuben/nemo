/*  multiscat.h - */

/*
 *  multiscat.h: header file for series of gravitational scattering experiments
 *		(c) 1989  Piet Hut  Princeton, NJ, USA
 */

#include  "scat.h"

/******  #define  BUFF_LENGTH    128  *******  this is now in file  diag.h   */
#define  MAX_DIVERSITY  256     /* maximum number of qualitatively different */
                                /* different types of outcome of a series of */
                                /* scattering experiments.                   */
/*-----------------------------------------------------------------------------
 *  multiclass  --  
 *-----------------------------------------------------------------------------
 */
typedef struct
    {
    char  hier_string[MAX_DIVERSITY][BUFF_LENGTH];
    int  subtotal[MAX_DIVERSITY];
    int  nscattype;
    } multiclass, *multiclassptr;

/*-----------------------------------------------------------------------------
 *  macros to extract individual (sub)components from a multi_class; e.g.:
 *
 *  Scattype(ptr)[3]       is the string which contains the 4th type of all
 *                         qualitatively different outcomes.
 *  Scatsubtotal(ptr)[3]   gives the total number of scattering experiments
 *                         which had an outcome of the 4th type.
 *  Nscattype(ptr)         gives the number of qualitatively different
 *                         outcomes encountered during the whole series of
 *                         scattering experiments.
 *-----------------------------------------------------------------------------
 */
#define  Scattype(ptr)       ((ptr)->hier_string)         /* type: stringptr */
#define  Scatsubtotal(ptr)   ((ptr)->subtotal)            /* type: intptr    */
#define  Nscattype(ptr)      ((ptr)->nscattype)           /* type: int       */

/* endof: multiscat.h */

