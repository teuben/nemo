/* potential0.c - get_potential, get_pattern, get_inipotential */

/*------------------------------------------------------------------------------
 *  POTENTIAL:  
 *	23-feb-93	Created - for standalone potential(5NEMO) usage
 *			although one has to recompile a program with
 *			this module for each potential to be used,
 *			the parameter passing (potpars,potfile) still
 *			works as advertised.
 *	11-oct-93	added in local_omega
 *------------------------------------------------------------------------------
 */

#include  <stdinc.h>
#include  <getparam.h>
#include  <loadobj.h>
#include "potential.h"


#define MAXPAR 64

local real local_par[MAXPAR];   /* NOTE: first par is reserved for pattern speed*/
local int  local_npar=0;        /* actual used number of par's                  */
local real local_omega=0;	/* pattern speed                             */
local proc l_potential=NULL;    /* actual storage of pointer to external worker */
local proc l_inipotential=NULL; /* actual storage of pointer to external inits  */

extern void inipotential(), potential();   /* fixed compilation, no loadobj */

/*-----------------------------------------------------------------------------
 *  get_potential --  returns the pointer ptr to the function which carries out
 *                the calculation of potential and accelerations.
 *-----------------------------------------------------------------------------
 */
proc  get_potential(potname, potpars, potfile)
string  potname;
string  potpars;
string  potfile;
{
    dprintf(0,"[Fixed potential0 potname=%s",potname);
    if (potname == NULL || *potname == NULL)	/* if no name provided */
        return NULL;				/* return no potential */
    local_npar = nemoinpd(potpars,local_par,MAXPAR);
    if (local_npar>MAXPAR)
        error ("get_potential: potential has too many parameters (%d)",
                            local_npar);
    if (local_npar<0)
        warning("get_potential: parsing error in: %s",potpars);

    if (local_npar > 0 && local_par[0] != 0.0) { /* aid multiple potentials */
        local_omega = local_par[0];
        dprintf(1,"get_potential: setting local_omega = %g\n",local_omega);
    }

    inipotential(&local_npar, local_par, potfile);
    if (local_npar > 0 && local_par[0] != local_omega) {
    	local_omega = local_par[0];
    	dprintf(1,"get_potential: modified omega=%g\n",local_omega);
    }

    l_potential = potential;
    return l_potential;
}

/*-----------------------------------------------------------------------------
 *  get_inipotential --  returns the pointer ptr to the last inipotential
 *          function which initializes the potential
 *-----------------------------------------------------------------------------
 */
proc get_inipotential()
{
    return l_inipotential;
}

/*-----------------------------------------------------------------------------
 *  get_pattern --  return the last set pattern speed
 *		    note that a 0.0 parameter does *not* reset it !!!
 *-----------------------------------------------------------------------------
 */

real get_pattern()
{
    return local_omega;
}

/* endof: potential.c */
