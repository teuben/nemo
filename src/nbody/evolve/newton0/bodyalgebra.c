/*
 *  bodyalgebra.c:  for algebraic operations on bodies in newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  "newton0.h"

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |   bodyalgebra.c :                                                   |  */
/*  |                                                                     |  */
/*  |           This file contains procedures for operations on xbodies.  |  */
/*  |                                                                     |  */
/*  |           xbodies stands for: bodies, ebodies, mbodies, pbodies,    |  */
/*  |                               and cbodies.                          |  */
/*  |                                                                     |  */
/*  |           operations stands for: creating, clearing, copying,       |  */
/*  |                                  multiplying, adding, subtracting,  |  */
/*  |                                  scaling, incrementing and          |  */
/*  |                                  decrementing.                      |  */
/*  |                                                                     |  */
/*  |           these operations preserve the type of xbody;              |  */
/*  |           other operations which mix different types of x           |  */
/*  |           reside in the file  bodyconversion.c .                    |  */
/*  |                                                                     |  */
/*  |           meaning of operations:                                    |  */
/*  |                                                                     |  */
/*  |           multiplying: take a xbody, multiply all components        |  */
/*  |                        with a constant scale factor, and            |  */
/*  |                        assign the resulting xbody to a newly        |  */
/*  |                        created xbody whose pointer is returned;     |  */
/*  |                        the old body remains unchanged.              |  */
/*  |           scaling: as multiplying, but now the results are          |  */
/*  |                    stored in the old body, and nothing is           |  */
/*  |                    returned explicitly by the procedure.            |  */
/*  |           adding: add two xbodies component-wise, assigning         |  */
/*  |                   the result to a new xbodies a pointer to          |  */
/*  |                   which is returned; the old xbodies remain         |  */
/*  |                   unchanged.                                        |  */
/*  |           incrementing: as adding, but now the results are          |  */
/*  |                         stored in the first of the two bodiess;     |  */
/*  |                         the other bodies remains unchanged,         |  */
/*  |                         and nothing is returned.                    |  */
/*  |           subtracting: analog of adding.                            |  */
/*  |           decrementing: analog of incrementing.                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |     TABLE OF CONTENTS:                                              |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART I                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on bodies:             |  */
/*  |                                                                     |  */
/*  |               bodyptr  mk_bodies(              npart )              |  */
/*  |               void     clear_bodies(     bod , npart )              |  */
/*  |               bodyptr  cp_bodies(        bod , npart )              |  */
/*  |               bodyptr  mul_bodies(       bod , npart, factor )      |  */
/*  |               bodyptr  add_bodies( bod1, bod2, npart )              |  */
/*  |               bodyptr  sub_bodies( bod1, bod2, npart )              |  */
/*  |               void     scale_bodies(     bod , npart, factor )      |  */
/*  |               void     inc_bodies( bod ,dbod , npart )              |  */
/*  |               void     dec_bodies( bod ,dbod , npart )              |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART II                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on ebodies:            |  */
/*  |                                                                     |  */
/*  |               ebodyptr  mk_ebodies(                npart )          |  */
/*  |               void      clear_ebodies(      ebod , npart )          |  */
/*  |               ebodyptr  cp_ebodies(         ebod , npart )          |  */
/*  |               ebodyptr  mul_ebodies(        ebod , npart, factor )  |  */
/*  |               ebodyptr  add_ebodies( ebod1, ebod2, npart )          |  */
/*  |               ebodyptr  sub_ebodies( ebod1, ebod2, npart )          |  */
/*  |               void      scale_ebodies(      ebod , npart, factor )  |  */
/*  |               void      inc_ebodies( ebod, debod , npart )          |  */
/*  |               void      dec_ebodies( ebod, debod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART III                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on mbodies:            |  */
/*  |                                                                     |  */
/*  |               mbodyptr  mk_mbodies(                npart )          |  */
/*  |               void      clear_mbodies(      mbod , npart )          |  */
/*  |               mbodyptr  cp_mbodies(         mbod , npart )          |  */
/*  |               mbodyptr  mul_mbodies(        mbod , npart, factor )  |  */
/*  |               mbodyptr  add_mbodies( mbod1, mbod2, npart )          |  */
/*  |               mbodyptr  sub_mbodies( mbod1, mbod2, npart )          |  */
/*  |               void      scale_mbodies(      mbod , npart, factor )  |  */
/*  |               void      inc_mbodies( mbod, dmbod , npart )          |  */
/*  |               void      dec_mbodies( mbod, dmbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                   (continued)       |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |     TABLE OF CONTENTS (continued):                                  |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART IV                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on pbodies:            |  */
/*  |                                                                     |  */
/*  |               pbodyptr  mk_pbodies(                npart )          |  */
/*  |               void      clear_pbodies(      pbod , npart )          |  */
/*  |               pbodyptr  cp_pbodies(         pbod , npart )          |  */
/*  |               pbodyptr  mul_pbodies(        pbod , npart, factor )  |  */
/*  |               pbodyptr  add_pbodies( pbod1, pbod2, npart )          |  */
/*  |               pbodyptr  sub_pbodies( pbod1, pbod2, npart )          |  */
/*  |               void      scale_pbodies(      pbod , npart, factor )  |  */
/*  |               void      inc_pbodies( pbod, dpbod , npart )          |  */
/*  |               void      dec_pbodies( pbod, dpbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART V                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on cbodies:            |  */
/*  |                                                                     |  */
/*  |               cbodyptr  mk_cbodies(                npart )          |  */
/*  |               void      clear_cbodies(      cbod , npart )          |  */
/*  |               cbodyptr  cp_cbodies(         cbod , npart )          |  */
/*  |               cbodyptr  mul_cbodies(        cbod , npart, factor )  |  */
/*  |               cbodyptr  add_cbodies( cbod1, cbod2, npart )          |  */
/*  |               cbodyptr  sub_cbodies( cbod1, cbod2, npart )          |  */
/*  |               void      scale_cbodies(      cbod , npart, factor )  |  */
/*  |               void      inc_cbodies( cbod, dcbod , npart )          |  */
/*  |               void      dec_cbodies( cbod, dcbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_bodies  --  allocates memory for bodies in standard newton0
 *                 form, and properly initializes the types of all bodies.
 *                 accepts: npart: the number of particles in the bodies.
 *                 returns: new_bod: a pointer to the new bodies.
 *-----------------------------------------------------------------------------
 */
bodyptr  mk_bodies(npart)
int  npart;
    {
    bodyptr  new_bod;
    bodyptr  new_i;      /* points to an individual particle of new_bod */

    new_bod = (bodyptr) malloc((unsigned)npart * sizeof(body));
    if (new_bod == NULL)
	error("mk_bodies: not enough memory left for a %d-particle bodies\n",
                                                                       npart);
    for (new_i = new_bod; new_i - new_bod < npart; new_i++)
        PartType(new_i) = BODY;

    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  clear_bodies  --  sets all values to zero in a bodies in standard newton0
 *                    form.
 *                    accepts: bod: pointer to bodies;
 *                           npart: the number of particles in that bodies.
 *-----------------------------------------------------------------------------
 */
void  clear_bodies(bod, npart)
bodyptr  bod;
int  npart;
    {
    bodyptr  body_i;      /* points to an individual particle of bod */
    
    for (body_i = bod; body_i - bod < npart; body_i++)
        {
        CLRV(Pos(body_i));
        CLRV(Vel(body_i));
	Mass(body_i) = 0.0;
	Pot(body_i) = 0.0;
        CLRV(Acc(body_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  cp_bodies  --  makes a copy of a bodies in standard newton0 form.
 *                 accepts: old_bod: pointer to bodies;
 *                            npart: the number of particles in that bodies.
 *                 returns: new_bod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
bodyptr  cp_bodies(old_bod, npart)
bodyptr  old_bod;
int  npart;
    {
    bodyptr  new_bod;
    bodyptr  new_i;      /* points to an individual particle of new_bod */
    bodyptr  old_i;      /* points to an individual particle of old_bod */
    
    new_bod = mk_bodies(npart);
    
    for (old_i = old_bod, new_i = new_bod; old_i - old_bod < npart;
                                                              old_i++, new_i++)
        {
        SETV(Pos(new_i), Pos(old_i));
        SETV(Vel(new_i), Vel(old_i));
	Mass(new_i) = Mass(old_i);
	Pot(new_i) = Pot(old_i);
        SETV(Acc(new_i), Acc(old_i));
        }
    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  mul_bodies  --  form a new bodies by multiplying an old bodies by a scalar.
 *                  accepts: old_bod: pointer to bodies in standard
 *                                    newton0 form.
 *                             npart: the number of particles in that bodies;
 *                      scale_factor: with which every entry in "old_bod"
 *                                    multiplied.
 *                  returns: new_bod: pointer to the new scaled copy.
 *                  note: EVERYTHING is scaled, even the masses, potentials and
 *                        accelerations; if this is not what you want, then use
 *                        body subset representations such as pbody or cbody.
 *-----------------------------------------------------------------------------
 */
bodyptr  mul_bodies(old_bod, npart, scale_factor)
bodyptr  old_bod;
int  npart;
real  scale_factor;
    {
    bodyptr  new_bod;
    bodyptr  new_i;      /* points to an individual particle of new_bod */
    bodyptr  old_i;      /* points to an individual particle of old_bod */
    
    new_bod = mk_bodies(npart);
    
    for (old_i = old_bod, new_i = new_bod; old_i - old_bod < npart;
                                                              old_i++, new_i++)
        {
        MULVS(Pos(new_i), Pos(old_i), scale_factor);
        MULVS(Vel(new_i), Vel(old_i), scale_factor);
	Mass(new_i) = scale_factor * Mass(old_i);
	Pot(new_i) = scale_factor * Pot(old_i);
        MULVS(Acc(new_i), Acc(old_i), scale_factor);
        }
    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  add_bodies  --  adds two bodiess of bodies together
 *                  accepts: bod1: pointer to bodies in standard newton0 form
 *                           bod2: pointer to second such bodies;
 *                           npart: the number of particles in each bodies.
 *                  returns: new_bod: pointer to new bodies which contains
 *                                    the sum of both old ones.
 *                  note: EVERYTHING is added, even the masses, potentials and
 *                        accelerations; if this is not what you want, then
 *                        either put these equal to zero in one of the bodiess
 *                        so that the sum inherits the values of the other
 *                        bodies, or use body subset representations such as
 *                        pbody or cbody.
 *-----------------------------------------------------------------------------
 */
bodyptr  add_bodies(bod1, bod2, npart)
bodyptr  bod1;
bodyptr  bod2;
int  npart;
    {
    bodyptr  new_bod;
    bodyptr  new_i, old1_i, old2_i;
    
    new_bod = mk_bodies(npart);
    
    for (new_i = new_bod, old1_i= bod1, old2_i = bod2;
                          new_i - new_bod < npart; new_i++, old1_i++, old2_i++)
        {
        ADDV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        ADDV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) + Mass(old2_i);
        Pot(new_i) = Pot(old1_i) + Pot(old2_i);
        ADDV(Acc(new_i), Acc(old1_i), Acc(old2_i));
	}    
    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  sub_bodies  --  subtracts two bodiess of bodies
 *                  accepts: bod1: pointer to bodies in standard newton0 form
 *                           bod2: pointer to second such bodies;
 *                           npart: the number of particles in each bodies.
 *                  returns: new_bod: pointer to new bodies which contains
 *                                    the difference of both old ones.
 *                  note: EVERYTHING is subtracted, even the masses, potentials
 *                        and accelerations; if this is not what you want, then
 *                        either put these equal to zero in one of the bodiess
 *                        so that the difference inherits the values of the
 *                        other bodies (modulo minus sign), or use body subset
 *                        representations such as pbody or cbody.
 *            disclaimer: the management is not responsible for such repulsive
 *                        creations as are formed by subtracting large masses
 *                        from small ones ...
 *-----------------------------------------------------------------------------
 */
bodyptr  sub_bodies(bod1, bod2, npart)
bodyptr  bod1;
bodyptr  bod2;
int  npart;
    {
    bodyptr  new_bod;
    bodyptr  new_i, old1_i, old2_i;
    
    new_bod = mk_bodies(npart);
    
    for (new_i = new_bod, old1_i= bod1, old2_i = bod2;
                          new_i - new_bod < npart; new_i++, old1_i++, old2_i++)
        {
        SUBV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        SUBV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) - Mass(old2_i);
        Pot(new_i) = Pot(old1_i) - Pot(old2_i);
        SUBV(Acc(new_i), Acc(old1_i), Acc(old2_i));
	}    
    return(new_bod);
    }

/*-----------------------------------------------------------------------------
 *  scale_bodies  --  multiplies an existing bodies by a scalar.
 *                    accepts: old_bod: pointer to bodies in
 *                                      standard newton0 form.
 *                               npart: the number of particles in that bodies;
 *                        scale_factor: with which every entry in "old_bod"
 *                                      multiplied.
 *                    note: EVERYTHING is scaled, even the masses, potentials
 *                          and accelerations; if this is not what you want,
 *                          then use body subset representations such as pbody
 *                          or cbody.
 *-----------------------------------------------------------------------------
 */
void  scale_bodies(bod, npart, scale_factor)
bodyptr  bod;
int  npart;
real  scale_factor;
    {
    bodyptr  body_i;      /* points to an individual particle of bod */
    
    for (body_i = bod; body_i - bod < npart; body_i++)
        {
        INCMULVS(Pos(body_i), scale_factor);
        INCMULVS(Vel(body_i), scale_factor);
	Mass(body_i) *= scale_factor;
	Pot(body_i) *= scale_factor;
        INCMULVS(Acc(body_i), scale_factor);
        }
    }

/*-----------------------------------------------------------------------------
 *  inc_bodies  --  increments a bodies by adding an other bodies
 *                  accepts: bod: pointer to bodies in standard newton0 form
 *                           d_bod: pointer to the bodies with which  bod
 *                                 will be incremented.
 *                           npart: the number of particles in each bodies.
 *                  note: EVERYTHING is added, even the masses, potentials and
 *                        accelerations; if this is not what you want, then
 *                        either put these equal to zero in one of the bodiess
 *                        so that the sum inherits the values of the other
 *                        bodies, or use body subset representations such as
 *                        pbody or cbody.
 *-----------------------------------------------------------------------------
 */
void  inc_bodies(bod, d_bod, npart)
bodyptr  bod;
bodyptr  d_bod;
int  npart;
    {
    bodyptr  body_i, d_body_i;
    
    for (body_i= bod, d_body_i = d_bod; body_i - bod < npart; 
                                                          body_i++, d_body_i++)
        {
        INCADDV(Pos(body_i), Pos(d_body_i));
        INCADDV(Vel(body_i), Vel(d_body_i));
        Mass(body_i) += Mass(d_body_i);
        Pot(body_i) += Pot(d_body_i);
        INCADDV(Acc(body_i), Acc(d_body_i));
	}    
    }

/*-----------------------------------------------------------------------------
 *  dec_bodies  --  decrements a bodies by subtracting an other bodies
 *                  accepts: bod: pointer to bodies in standard newton0 form
 *                           d_bod: pointer to the bodies with which  bod
 *                                 will be decremented.
 *                           npart: the number of particles in each bodies.
 *                  note: EVERYTHING is subtracted, even the masses, potentials
 *                        and accelerations; if this is not what you want, then
 *                        either put these equal to zero in one of the bodiess
 *                        so that the difference inherits the values of the
 *                        other bodies (modulo minus sign), or use body subset
 *                        representations such as pbody or cbody.
 *            disclaimer: the management is not responsible for such repulsive
 *                        creations as are formed by subtracting large masses
 *                        from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_bodies(bod, d_bod, npart)
bodyptr  bod;
bodyptr  d_bod;
int  npart;
    {
    bodyptr  body_i, d_body_i;
    
    for (body_i= bod, d_body_i = d_bod; body_i - bod < npart; 
                                                          body_i++, d_body_i++)
        {
        INCSUBV(Pos(body_i), Pos(d_body_i));
        INCSUBV(Vel(body_i), Vel(d_body_i));
        Mass(body_i) -= Mass(d_body_i);
        Pot(body_i) -= Pot(d_body_i);
        INCSUBV(Acc(body_i), Acc(d_body_i));
	}    
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART II                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on ebodies:            |  */
/*  |                                                                     |  */
/*  |               ebodyptr  mk_ebodies(                npart )          |  */
/*  |               void      clear_ebodies(      ebod , npart )          |  */
/*  |               ebodyptr  cp_ebodies(         ebod , npart )          |  */
/*  |               ebodyptr  mul_ebodies(        ebod , npart, factor )  |  */
/*  |               ebodyptr  add_ebodies( ebod1, ebod2, npart )          |  */
/*  |               ebodyptr  sub_ebodies( ebod1, ebod2, npart )          |  */
/*  |               void      scale_ebodies(      ebod , npart, factor )  |  */
/*  |               void      inc_ebodies( ebod, debod , npart )          |  */
/*  |               void      dec_ebodies( ebod, debod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_ebodies  --  allocates memory for ebodies in standard
 *                  newton0 form, and properly initializes the types of 
 *                  all bodies.
 *                  accepts: npart: the number of particles in the ebodies.
 *                  returns: new_ebod: a pointer to the new ebodies.
 *-----------------------------------------------------------------------------
 */
ebodyptr  mk_ebodies(npart)
int  npart;
    {
    ebodyptr  new_ebod;
    ebodyptr  new_e_i;      /* points to an individual particle of new_ebod */

    new_ebod = (ebodyptr) malloc((unsigned)npart * sizeof(ebody));
    if (new_ebod == NULL)
	error("mk_ebodies: not enough memory left for a %d-particle ebodies\n",
                                                                       npart);
    for (new_e_i = new_ebod; new_e_i - new_ebod < npart; new_e_i++)
        PartType(new_e_i) = EBODY;

    return(new_ebod);
    }

/*-----------------------------------------------------------------------------
 *  clear_ebodies  --  sets all values to zero in an ebodies in standard
 *                     newton0 form.
 *                     accepts: ebod: pointer to ebodies;
 *                             npart: the number of particles in that ebodies.
 *-----------------------------------------------------------------------------
 */
void  clear_ebodies(ebod, npart)
ebodyptr  ebod;
int  npart;
    {
    ebodyptr  ebody_i;      /* points to an individual particle of ebod */
    
    for (ebody_i = ebod; ebody_i - ebod < npart; ebody_i++)
        {
        CLRV(Pos(ebody_i));
        CLRV(Vel(ebody_i));
	Mass(ebody_i) = 0.0;
	Pot(ebody_i) = 0.0;
        }
    }

/*-----------------------------------------------------------------------------
 *  cp_ebodies  --  makes a copy of an ebodies in standard newton0 form.
 *                  accepts: old_ebod: pointer to ebodies;
 *                              npart: the number of particles in that ebodies.
 *                  returns: new_ebod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
ebodyptr  cp_ebodies(old_ebod, npart)
ebodyptr  old_ebod;
int  npart;
    {
    ebodyptr  new_ebod;
    ebodyptr  new_i;      /* points to an individual particle of new_ebod */
    ebodyptr  old_i;      /* points to an individual particle of old_ebod */
    
    new_ebod = mk_ebodies(npart);
    
    for (old_i = old_ebod, new_i = new_ebod; old_i - old_ebod < npart;
                                                              old_i++, new_i++)
        {
        SETV(Pos(new_i), Pos(old_i));
        SETV(Vel(new_i), Vel(old_i));
	Mass(new_i) = Mass(old_i);
	Pot(new_i) = Pot(old_i);
        }
    return(new_ebod);
    }

/*-----------------------------------------------------------------------------
 *  mul_ebodies  --  form a new ebodies by multiplying an old ebodies by a
 *                   scalar.
 *                   accepts: old_ebod: pointer to many-ebody ebodies in
 *                                      standard newton0 form.
 *                               npart: the number of particles in that
 *                                      ebodies;
 *                        scale_factor: with which every entry in "old_ebod"
 *                                      multiplied.
 *                   returns: new_ebod: pointer to the new scaled copy.
 *                   note: EVERYTHING is scaled, even the masses and
 *                         potentials; if this is not what you want, then use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *-----------------------------------------------------------------------------
 */
ebodyptr  mul_ebodies(old_ebod, npart, scale_factor)
ebodyptr  old_ebod;
int  npart;
real  scale_factor;
    {
    ebodyptr  new_ebod;
    ebodyptr  new_i;      /* points to an individual particle of new_ebod */
    ebodyptr  old_i;      /* points to an individual particle of old_ebod */
    
    new_ebod = mk_ebodies(npart);
    
    for (old_i = old_ebod, new_i = new_ebod; old_i - old_ebod < npart;
                                                              old_i++, new_i++)
        {
        MULVS(Pos(new_i), Pos(old_i), scale_factor);
        MULVS(Vel(new_i), Vel(old_i), scale_factor);
	Mass(new_i) = scale_factor * Mass(old_i);
	Pot(new_i) = scale_factor * Pot(old_i);
        }
    return(new_ebod);
    }

/*-----------------------------------------------------------------------------
 *  add_ebodies  --  adds two ebodiess of ebodies together
 *                   accepts: ebod1: pointer to ebodies in standard newton0
 *                                   form;
 *                            ebod2: pointer to second such ebodies;
 *                            npart: the number of particles in each ebodies.
 *                   returns: new_ebod: pointer to new ebodies which contains
 *                                     the sum of both old ones.
 *                   note: EVERYTHING is added, even the masses and potentials;
 *                         if this is not what you want, then either put these
 *                         equal to zero in one of the ebodiess so that the sum
 *                         inherits the values of the other ebodies, or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *-----------------------------------------------------------------------------
 */
ebodyptr  add_ebodies(ebod1, ebod2, npart)
ebodyptr  ebod1;
ebodyptr  ebod2;
int  npart;
    {
    ebodyptr  new_ebod;
    ebodyptr  new_i, old1_i, old2_i;
    
    new_ebod = mk_ebodies(npart);
    
    for (new_i = new_ebod, old1_i= ebod1, old2_i = ebod2;
                         new_i - new_ebod < npart; new_i++, old1_i++, old2_i++)
        {
        ADDV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        ADDV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) + Mass(old2_i);
        Pot(new_i) = Pot(old1_i) + Pot(old2_i);
	}    
    return(new_ebod);
    }

/*-----------------------------------------------------------------------------
 *  sub_ebodies  --  subtracts two ebodiess of ebodies
 *                   accepts: ebod1: pointer to ebodies in standard newton0
 *                                   form;
 *                            ebod2: pointer to second such ebodies;
 *                            npart: the number of particles in each ebodies.
 *                   returns: new_ebod: pointer to new ebodies which contains
 *                                      the difference of both old ones.
 *                   note: EVERYTHING is subtracted, even the masses and
 *                         potentials; if this is not what you want, then
 *                         either put these equal to zero in one of the
 *                         ebodiess so that the difference inherits the values
 *                         of the other ebodies (modulo minus sign), or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
ebodyptr  sub_ebodies(ebod1, ebod2, npart)
ebodyptr  ebod1;
ebodyptr  ebod2;
int  npart;
    {
    ebodyptr  new_ebod;
    ebodyptr  new_i, old1_i, old2_i;
    
    new_ebod = mk_ebodies(npart);
    
    for (new_i = new_ebod, old1_i= ebod1, old2_i = ebod2;
                         new_i - new_ebod < npart; new_i++, old1_i++, old2_i++)
        {
        SUBV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        SUBV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) - Mass(old2_i);
        Pot(new_i) = Pot(old1_i) - Pot(old2_i);
	}    
    return(new_ebod);
    }

/*-----------------------------------------------------------------------------
 *  scale_ebodies  --  multiplies an existing ebodies by a scalar.
 *                     accepts: old_ebod: pointer to ebodies in
 *                                        standard newton0 form.
 *                                 npart: the number of particles in that
 *                                        ebodies;
 *                          scale_factor: with which every entry in "old_ebod"
 *                                        multiplied.
 *                     note: EVERYTHING is scaled, even the masses and
 *                           potentials; if this is not what you want, then use
 *                           other body subset representations such as pbody or
 *                           cbody.
 *-----------------------------------------------------------------------------
 */
void  scale_ebodies(ebod, npart, scale_factor)
ebodyptr  ebod;
int  npart;
real  scale_factor;
    {
    ebodyptr  ebody_i;      /* points to an individual particle of ebod */
    
    for (ebody_i = ebod; ebody_i - ebod < npart; ebody_i++)
        {
        INCMULVS(Pos(ebody_i), scale_factor);
        INCMULVS(Vel(ebody_i), scale_factor);
	Mass(ebody_i) *= scale_factor;
	Pot(ebody_i) *= scale_factor;
        }
    }

/*-----------------------------------------------------------------------------
 *  inc_ebodies  --  increments an ebodies by adding an other ebodies
 *                   accepts: ebod: pointer to ebodies in standard newton0
 *                                   form;
 *                            d_ebod: pointer to the ebodies with which  ebod
 *                                   will be incremented.
 *                            npart: the number of particles in each ebodies.
 *                   note: EVERYTHING is added, even the masses and potentials;
 *                         if this is not what you want, then either put these
 *                         equal to zero in one of the ebodiess so that the sum
 *                         inherits the values of the other ebodies, or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *-----------------------------------------------------------------------------
 */
void  inc_ebodies(ebod, d_ebod, npart)
ebodyptr  ebod;
ebodyptr  d_ebod;
int  npart;
    {
    ebodyptr  ebody_i, d_ebody_i;
    
    for (ebody_i= ebod, d_ebody_i = d_ebod; ebody_i - ebod < npart; 
                                                        ebody_i++, d_ebody_i++)
        {
        INCADDV(Pos(ebody_i), Pos(d_ebody_i));
        INCADDV(Vel(ebody_i), Vel(d_ebody_i));
        Mass(ebody_i) += Mass(d_ebody_i);
        Pot(ebody_i) += Pot(d_ebody_i);
	}    
    }

/*-----------------------------------------------------------------------------
 *  dec_ebodies  --  decrements an ebodies by subtracting an other ebodies
 *                   accepts: ebod: pointer to ebodies in standard newton0
 *                                   form;
 *                            d_ebod: pointer to the ebodies with which  ebod
 *                                   will be decremented.
 *                            npart: the number of particles in each ebodies.
 *                   note: EVERYTHING is subtracted, even the masses and
 *                         potentials; if this is not what you want, then
 *                         either put these equal to zero in one of the
 *                         ebodiess so that the difference inherits the values
 *                         of the other ebodies (modulo minus sign), or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_ebodies(ebod, d_ebod, npart)
ebodyptr  ebod;
ebodyptr  d_ebod;
int  npart;
    {
    ebodyptr  ebody_i, d_ebody_i;
    
    for (ebody_i= ebod, d_ebody_i = d_ebod; ebody_i - ebod < npart; 
                                                        ebody_i++, d_ebody_i++)
        {
        INCSUBV(Pos(ebody_i), Pos(d_ebody_i));
        INCSUBV(Vel(ebody_i), Vel(d_ebody_i));
        Mass(ebody_i) -= Mass(d_ebody_i);
        Pot(ebody_i) -= Pot(d_ebody_i);
	}    
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART III                               |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on mbodies.            |  */
/*  |                                                                     |  */
/*  |               mbodyptr  mk_mbodies(                npart )          |  */
/*  |               void      clear_mbodies(      mbod , npart )          |  */
/*  |               mbodyptr  cp_mbodies(         mbod , npart )          |  */
/*  |               mbodyptr  mul_mbodies(        mbod , npart, factor )  |  */
/*  |               mbodyptr  add_mbodies( mbod1, mbod2, npart )          |  */
/*  |               mbodyptr  sub_mbodies( mbod1, mbod2, npart )          |  */
/*  |               void      scale_mbodies(      mbod , npart, factor )  |  */
/*  |               void      inc_mbodies( mbod, dmbod , npart )          |  */
/*  |               void      dec_mbodies( mbod, dmbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_mbodies  --  allocates memory for mbodies in standard
 *                  newton0 form, and properly initializes the types of 
 *                  all bodies.
 *                  accepts: npart: the number of particles in the mbodies.
 *                  returns: new_mbod: a pointer to the new mbodies.
 *-----------------------------------------------------------------------------
 */
mbodyptr  mk_mbodies(npart)
int  npart;
    {
    mbodyptr  new_mbod;
    mbodyptr  new_m_i;      /* points to an individual particle of new_mbod */

    new_mbod = (mbodyptr) malloc((unsigned)npart * sizeof(mbody));
    if (new_mbod == NULL)
	error("mk_mbodies: not enough memory left for a %d-particle mbodies\n",
                                                                       npart);
    for (new_m_i = new_mbod; new_m_i - new_mbod < npart; new_m_i++)
        PartType(new_m_i) = MBODY;

    return(new_mbod);
    }

/*-----------------------------------------------------------------------------
 *  clear_mbodies  --  sets all values to zero in a mbodies in standard newton0
 *                     form.
 *                     accepts: mbod: pointer to mbodies;
 *                             npart: the number of particles in that mbodies.
 *-----------------------------------------------------------------------------
 */
void  clear_mbodies(mbod, npart)
mbodyptr  mbod;
int  npart;
    {
    mbodyptr  mbody_i;      /* points to an individual particle of mbod */
    
    for (mbody_i = mbod; mbody_i - mbod < npart; mbody_i++)
        {
        CLRV(Pos(mbody_i));
        CLRV(Vel(mbody_i));
	Mass(mbody_i) = 0.0;
        }
    }

/*-----------------------------------------------------------------------------
 *  cp_mbodies  --  makes a copy of a mbodies in standard newton0 form.
 *                  accepts: old_mbod: pointer to mbodies;
 *                              npart: the number of particles in that mbodies.
 *                  returns: new_mbod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
mbodyptr  cp_mbodies(old_mbod, npart)
mbodyptr  old_mbod;
int  npart;
    {
    mbodyptr  new_mbod;
    mbodyptr  new_i;      /* points to an individual particle of new_mbod */
    mbodyptr  old_i;      /* points to an individual particle of old_mbod */
    
    new_mbod = mk_mbodies(npart);
    
    for (old_i = old_mbod, new_i = new_mbod; old_i - old_mbod < npart;
                                                              old_i++, new_i++)
        {
        SETV(Pos(new_i), Pos(old_i));
        SETV(Vel(new_i), Vel(old_i));
	Mass(new_i) = Mass(old_i);
        }
    return(new_mbod);
    }

/*-----------------------------------------------------------------------------
 *  mul_mbodies  --  form a new mbodies by multiplying an old mbodies by a
 *                   scalar.
 *                   accepts: old_mbod: pointer to many-mbody mbodies in
 *                                      standard newton0 form.
 *                               npart: the number of particles in that
 *                                      mbodies;
 *                        scale_factor: with which every entry in "old_mbod"
 *                                      multiplied.
 *                   returns: new_mbod: pointer to the new scaled copy.
 *                   note: EVERYTHING is scaled, even the masses;
 *                         if this is not what you want, then use other body
 *                         subset representations such as pbody or cbody.
 *-----------------------------------------------------------------------------
 */
mbodyptr  mul_mbodies(old_mbod, npart, scale_factor)
mbodyptr  old_mbod;
int  npart;
real  scale_factor;
    {
    mbodyptr  new_mbod;
    mbodyptr  new_i;      /* points to an individual particle of new_mbod */
    mbodyptr  old_i;      /* points to an individual particle of old_mbod */
    
    new_mbod = mk_mbodies(npart);
    
    for (old_i = old_mbod, new_i = new_mbod; old_i - old_mbod < npart;
                                                              old_i++, new_i++)
        {
        MULVS(Pos(new_i), Pos(old_i), scale_factor);
        MULVS(Vel(new_i), Vel(old_i), scale_factor);
	Mass(new_i) = scale_factor * Mass(old_i);
        }
    return(new_mbod);
    }

/*-----------------------------------------------------------------------------
 *  add_mbodies  --  adds two mbodiess of mbodies together
 *                   accepts: mbod1: pointer to mbodies in standard newton0
 *                                   form;
 *                            mbod2: pointer to second such mbodies;
 *                            npart: the number of particles in each mbodies.
 *                   returns: new_mbod: pointer to new mbodies which contains
 *                                     the sum of both old ones.
 *                   note: EVERYTHING is added, even the masses;
 *                         if this is not what you want, then either put these
 *                         equal to zero in one of the mbodiess so that the sum
 *                         inherits the values of the other mbodies, or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *-----------------------------------------------------------------------------
 */
mbodyptr  add_mbodies(mbod1, mbod2, npart)
mbodyptr  mbod1;
mbodyptr  mbod2;
int  npart;
    {
    mbodyptr  new_mbod;
    mbodyptr  new_i, old1_i, old2_i;
    
    new_mbod = mk_mbodies(npart);
    
    for (new_i = new_mbod, old1_i= mbod1, old2_i = mbod2;
                         new_i - new_mbod < npart; new_i++, old1_i++, old2_i++)
        {
        ADDV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        ADDV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) + Mass(old2_i);
	}    
    return(new_mbod);
    }

/*-----------------------------------------------------------------------------
 *  sub_mbodies  --  subtracts two mbodiess of mbodies
 *                   accepts: mbod1: pointer to mbodies in standard newton0
 *                                   form;
 *                            mbod2: pointer to second such mbodies;
 *                            npart: the number of particles in each mbodies.
 *                   returns: new_mbod: pointer to new mbodies which contains
 *                                      the difference of both old ones.
 *                   note: EVERYTHING is subtracted, even the masses; if this
 *                         is not what you want, then either put these equal to
 *                         zero in one of the mbodiess so that the difference
 *                         inherits the values of the other mbodies (modulo
 *                         minus sign), or use other body subset
 *                         representations such as pbody or cbody.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
mbodyptr  sub_mbodies(mbod1, mbod2, npart)
mbodyptr  mbod1;
mbodyptr  mbod2;
int  npart;
    {
    mbodyptr  new_mbod;
    mbodyptr  new_i, old1_i, old2_i;
    
    new_mbod = mk_mbodies(npart);
    
    for (new_i = new_mbod, old1_i= mbod1, old2_i = mbod2;
                         new_i - new_mbod < npart; new_i++, old1_i++, old2_i++)
        {
        SUBV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        SUBV(Vel(new_i), Vel(old1_i), Vel(old2_i));
        Mass(new_i) = Mass(old1_i) - Mass(old2_i);
	}    
    return(new_mbod);
    }

/*-----------------------------------------------------------------------------
 *  scale_mbodies  --  multiplies an existing mbodies by a scalar.
 *                     accepts: old_mbod: pointer to mbodies in
 *                                       standard newton0 form.
 *                                 npart: the number of particles in that
 *                                        mbodies;
 *                        scale_factor: with which every entry in "old_mbod"
 *                                      multiplied.
 *                  note: EVERYTHING is scaled, even the masses;
 *                        if this is not what you want, then use other body
 *                        subset representations such as pbody or cbody.
 *-----------------------------------------------------------------------------
 */
void  scale_mbodies(mbod, npart, scale_factor)
mbodyptr  mbod;
int  npart;
real  scale_factor;
    {
    mbodyptr  mbody_i;      /* points to an individual particle of mbod */
    
    for (mbody_i = mbod; mbody_i - mbod < npart; mbody_i++)
        {
        INCMULVS(Pos(mbody_i), scale_factor);
        INCMULVS(Vel(mbody_i), scale_factor);
	Mass(mbody_i) *= scale_factor;
        }
    }

/*-----------------------------------------------------------------------------
 *  inc_mbodies  --  increments a mbodies by adding an other mbodies
 *                   accepts: mbod: pointer to mbodies in standard newton0
 *                                   form;
 *                            d_mbod: pointer to the mbodies with which  mbod
 *                                   will be incremented.
 *                            npart: the number of particles in each mbodies.
 *                   note: EVERYTHING is added, even the masses;
 *                         if this is not what you want, then either put these
 *                         equal to zero in one of the mbodiess so that the sum
 *                         inherits the values of the other mbodies, or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *-----------------------------------------------------------------------------
 */
void  inc_mbodies(mbod, d_mbod, npart)
mbodyptr  mbod;
mbodyptr  d_mbod;
int  npart;
    {
    mbodyptr  mbody_i, d_mbody_i;
    
    for (mbody_i= mbod, d_mbody_i = d_mbod; mbody_i - mbod < npart; 
                                                        mbody_i++, d_mbody_i++)
        {
        INCADDV(Pos(mbody_i), Pos(d_mbody_i));
        INCADDV(Vel(mbody_i), Vel(d_mbody_i));
        Mass(mbody_i) += Mass(d_mbody_i);
	}    
    }

/*-----------------------------------------------------------------------------
 *  dec_mbodies  --  decrements a mbodies by subtracting an other mbodies
 *                   accepts: mbod: pointer to mbodies in standard newton0
 *                                   form;
 *                            d_mbod: pointer to the mbodies with which  mbod
 *                                   will be decremented.
 *                            npart: the number of particles in each mbodies.
 *                   note: EVERYTHING is subtracted, even the masses;
 *                         if this is not what you want, then either put these
 *                         equal to zero in one of the mbodiess so that the sum
 *                         inherits the values of the other mbodies, or use
 *                         other body subset representations such as pbody or
 *                         cbody.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_mbodies(mbod, d_mbod, npart)
mbodyptr  mbod;
mbodyptr  d_mbod;
int  npart;
    {
    mbodyptr  mbody_i, d_mbody_i;
    
    for (mbody_i= mbod, d_mbody_i = d_mbod; mbody_i - mbod < npart; 
                                                        mbody_i++, d_mbody_i++)
        {
        INCSUBV(Pos(mbody_i), Pos(d_mbody_i));
        INCSUBV(Vel(mbody_i), Vel(d_mbody_i));
        Mass(mbody_i) -= Mass(d_mbody_i);
	}    
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART IV                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on pbodies.            |  */
/*  |                                                                     |  */
/*  |               pbodyptr  mk_pbodies(                npart )          |  */
/*  |               void      clear_pbodies(      pbod , npart )          |  */
/*  |               pbodyptr  cp_pbodies(         pbod , npart )          |  */
/*  |               pbodyptr  mul_pbodies(        pbod , npart, factor )  |  */
/*  |               pbodyptr  add_pbodies( pbod1, pbod2, npart )          |  */
/*  |               pbodyptr  sub_pbodies( pbod1, pbod2, npart )          |  */
/*  |               void      scale_pbodies(      pbod , npart, factor )  |  */
/*  |               void      inc_pbodies( pbod, dpbod , npart )          |  */
/*  |               void      dec_pbodies( pbod, dpbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_pbodies  --  allocates memory for pbodies in standard
 *                  newton0 form, and properly initializes the types of 
 *                  all bodies.
 *                  accepts: npart: the number of particles in the pbodies.
 *                  returns: new_pbod: a pointer to the new pbodies.
 *-----------------------------------------------------------------------------
 */
pbodyptr  mk_pbodies(npart)
int  npart;
    {
    pbodyptr  new_pbod;
    pbodyptr  new_p_i;      /* points to an individual particle of new_pbod */

    new_pbod = (pbodyptr) malloc((unsigned)npart * sizeof(pbody));
    if (new_pbod == NULL)
	error("mk_pbodies: not enough memory left for a %d-particle pbodies\n",
                                                                       npart);
    for (new_p_i = new_pbod; new_p_i - new_pbod < npart; new_p_i++)
        PartType(new_p_i) = PBODY;

    return(new_pbod);
    }

/*-----------------------------------------------------------------------------
 *  clear_pbodies  --  sets all values to zero in a pbodies in standard newton0
 *                     form.
 *                     accepts: pbod: pointer to pbodies;
 *                             npart: the number of particles in that pbodies.
 *-----------------------------------------------------------------------------
 */
void  clear_pbodies(pbod, npart)
pbodyptr  pbod;
int  npart;
    {
    pbodyptr  pbody_i;      /* points to an individual particle of pbod */
    
    for (pbody_i = pbod; pbody_i - pbod < npart; pbody_i++)
        {
        CLRV(Pos(pbody_i));
        CLRV(Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  cp_pbodies  --  makes a copy of a pbodies in standard newton0 form.
 *                  accepts: old_pbod: pointer to pbodies;
 *                              npart: the number of particles in that pbodies.
 *                  returns: new_pbod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
pbodyptr  cp_pbodies(old_pbod, npart)
pbodyptr  old_pbod;
int  npart;
    {
    pbodyptr  new_pbod;
    pbodyptr  new_i;      /* points to an individual particle of new_pbod */
    pbodyptr  old_i;      /* points to an individual particle of old_pbod */
    
    new_pbod = mk_pbodies(npart);
    
    for (old_i = old_pbod, new_i = new_pbod; old_i - old_pbod < npart;
                                                              old_i++, new_i++)
        {
        SETV(Pos(new_i), Pos(old_i));
        SETV(Vel(new_i), Vel(old_i));
        }
    return(new_pbod);
    }

/*-----------------------------------------------------------------------------
 *  mul_pbodies  --  form a new pbodies by multiplying an old pbodies by a
 *                   scalar.
 *                   accepts: old_pbod: pointer to many-pbody pbodies in
 *                                      standard newton0 form.
 *                               npart: the number of particles in that
 *                                      pbodies;
 *                        scale_factor: with which every entry in "old_pbod"
 *                                      multiplied.
 *                   returns: new_pbod: pointer to the new scaled copy.
 *-----------------------------------------------------------------------------
 */
pbodyptr  mul_pbodies(old_pbod, npart, scale_factor)
pbodyptr  old_pbod;
int  npart;
real  scale_factor;
    {
    pbodyptr  new_pbod;
    pbodyptr  new_i;      /* points to an individual particle of new_pbod */
    pbodyptr  old_i;      /* points to an individual particle of old_pbod */
    
    new_pbod = mk_pbodies(npart);
    
    for (old_i = old_pbod, new_i = new_pbod; old_i - old_pbod < npart;
                                                              old_i++, new_i++)
        {
        MULVS(Pos(new_i), Pos(old_i), scale_factor);
        MULVS(Vel(new_i), Vel(old_i), scale_factor);
        }
    return(new_pbod);
    }

/*-----------------------------------------------------------------------------
 *  add_pbodies  --  adds two pbodiess of mbodies together
 *                   accepts: pbod1: pointer to pbodies in standard newton0
 *                                   form;
 *                            pbod2: pointer to second such pbodies;
 *                            npart: the number of particles in each pbodies.
 *                   returns: new_pbod: pointer to new pbodies which contains
 *                                     the sum of both old ones.
 *-----------------------------------------------------------------------------
 */
pbodyptr  add_pbodies(pbod1, pbod2, npart)
pbodyptr  pbod1;
pbodyptr  pbod2;
int  npart;
    {
    pbodyptr  new_pbod;
    pbodyptr  new_i, old1_i, old2_i;
    
    new_pbod = mk_pbodies(npart);
    
    for (new_i = new_pbod, old1_i= pbod1, old2_i = pbod2;
                         new_i - new_pbod < npart; new_i++, old1_i++, old2_i++)
        {
        ADDV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        ADDV(Vel(new_i), Vel(old1_i), Vel(old2_i));
	}    
    return(new_pbod);
    }

/*-----------------------------------------------------------------------------
 *  sub_pbodies  --  subtracts two pbodiess of mbodies
 *                   accepts: pbod1: pointer to pbodies in standard newton0
 *                                   form;
 *                            pbod2: pointer to second such pbodies;
 *                            npart: the number of particles in each pbodies.
 *                   returns: new_pbod: pointer to new pbodies which contains
 *                                      the difference of both old ones.
 *-----------------------------------------------------------------------------
 */
pbodyptr  sub_pbodies(pbod1, pbod2, npart)
pbodyptr  pbod1;
pbodyptr  pbod2;
int  npart;
    {
    pbodyptr  new_pbod;
    pbodyptr  new_i, old1_i, old2_i;
    
    new_pbod = mk_pbodies(npart);
    
    for (new_i = new_pbod, old1_i= pbod1, old2_i = pbod2;
                         new_i - new_pbod < npart; new_i++, old1_i++, old2_i++)
        {
        SUBV(Pos(new_i), Pos(old1_i), Pos(old2_i));
        SUBV(Vel(new_i), Vel(old1_i), Vel(old2_i));
	}    
    return(new_pbod);
    }

/*-----------------------------------------------------------------------------
 *  scale_pbodies  --  multiplies an existing pbodies by a scalar.
 *                     accepts: old_pbod: pointer to pbodies in
 *                                       standard newton0 form.
 *                                 npart: the number of particles in that
 *                                        pbodies;
 *                        scale_factor: with which every entry in "old_pbod"
 *                                      multiplied.
 *-----------------------------------------------------------------------------
 */
void  scale_pbodies(pbod, npart, scale_factor)
pbodyptr  pbod;
int  npart;
real  scale_factor;
    {
    pbodyptr  pbody_i;      /* points to an individual particle of pbod */
    
    for (pbody_i = pbod; pbody_i - pbod < npart; pbody_i++)
        {
        INCMULVS(Pos(pbody_i), scale_factor);
        INCMULVS(Vel(pbody_i), scale_factor);
        }
    }

/*-----------------------------------------------------------------------------
 *  inc_pbodies  --  increments a pbodies by adding an other pbodies
 *                   accepts: pbod: pointer to pbodies in standard newton0
 *                                   form;
 *                            d_pbod: pointer to the pbodies with which  pbod
 *                                   will be incremented.
 *                            npart: the number of particles in each pbodies.
 *-----------------------------------------------------------------------------
 */
void  inc_pbodies(pbod, d_pbod, npart)
pbodyptr  pbod;
pbodyptr  d_pbod;
int  npart;
    {
    pbodyptr  pbody_i, d_pbody_i;
    
    for (pbody_i= pbod, d_pbody_i = d_pbod; pbody_i - pbod < npart; 
                                                        pbody_i++, d_pbody_i++)
        {
        INCADDV(Pos(pbody_i), Pos(d_pbody_i));
        INCADDV(Vel(pbody_i), Vel(d_pbody_i));
	}    
    }

/*-----------------------------------------------------------------------------
 *  dec_pbodies  --  decrements a pbodies by subtracting an other pbodies
 *                   accepts: pbod: pointer to pbodies in standard newton0
 *                                   form;
 *                            d_pbod: pointer to the pbodies with which  pbod
 *                                   will be decremented.
 *                            npart: the number of particles in each pbodies.
 *-----------------------------------------------------------------------------
 */
void  dec_pbodies(pbod, d_pbod, npart)
pbodyptr  pbod;
pbodyptr  d_pbod;
int  npart;
    {
    pbodyptr  pbody_i, d_pbody_i;
    
    for (pbody_i= pbod, d_pbody_i = d_pbod; pbody_i - pbod < npart; 
                                                        pbody_i++, d_pbody_i++)
        {
        INCSUBV(Pos(pbody_i), Pos(d_pbody_i));
        INCSUBV(Vel(pbody_i), Vel(d_pbody_i));
	}    
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART V                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on cbodies.            |  */
/*  |                                                                     |  */
/*  |               cbodyptr  mk_cbodies(                npart )          |  */
/*  |               void      clear_cbodies(      cbod , npart )          |  */
/*  |               cbodyptr  cp_cbodies(         cbod , npart )          |  */
/*  |               cbodyptr  mul_cbodies(        cbod , npart, factor )  |  */
/*  |               cbodyptr  add_cbodies( cbod1, cbod2, npart )          |  */
/*  |               cbodyptr  sub_cbodies( cbod1, cbod2, npart )          |  */
/*  |               void      scale_cbodies(      cbod , npart, factor )  |  */
/*  |               void      inc_cbodies( cbod, dcbod , npart )          |  */
/*  |               void      dec_cbodies( cbod, dcbod , npart )          |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_cbodies  --  allocates memory for cbodies in standard
 *                  newton0 form, and properly initializes the types of 
 *                  all bodies.
 *                  accepts: npart: the number of particles in the cbodies.
 *                  returns: new_cbod: a pointer to the new cbodies.
 *-----------------------------------------------------------------------------
 */
cbodyptr  mk_cbodies(npart)
int  npart;
    {
    cbodyptr  new_cbod;
    cbodyptr  new_c_i;      /* points to an individual particle of new_cbod */

    new_cbod = (cbodyptr) malloc((unsigned)npart * sizeof(cbody));
    if (new_cbod == NULL)
	error("mk_cbodies: not enough memory left for a %d-particle cbodies\n",
                                                                       npart);
    for (new_c_i = new_cbod; new_c_i - new_cbod < npart; new_c_i++)
        PartType(new_c_i) = CBODY;

    return(new_cbod);
    }

/*-----------------------------------------------------------------------------
 *  clear_cbodies  --  sets all values to zero in a cbodies in standard newton0
 *                     form.
 *                     accepts: cbod: pointer to cbodies;
 *                             npart: the number of particles in that cbodies.
 *-----------------------------------------------------------------------------
 */
void  clear_cbodies(cbod, npart)
cbodyptr  cbod;
int  npart;
    {
    cbodyptr  cbody_i;      /* points to an individual particle of cbod */
    
    for (cbody_i = cbod; cbody_i - cbod < npart; cbody_i++)
        CLRV(Config(cbody_i));
    }

/*-----------------------------------------------------------------------------
 *  cp_cbodies  --  makes a copy of a cbodies in standard newton0 form.
 *                  accepts: old_cbod: pointer to cbodies;
 *                              npart: the number of particles in that cbodies.
 *                  returns: new_cbod: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
cbodyptr  cp_cbodies(old_cbod, npart)
cbodyptr  old_cbod;
int  npart;
    {
    cbodyptr  new_cbod;
    cbodyptr  new_i;      /* points to an individual particle of new_cbod */
    cbodyptr  old_i;      /* points to an individual particle of old_cbod */
    
    new_cbod = mk_cbodies(npart);
    
    for (old_i = old_cbod, new_i = new_cbod; old_i - old_cbod < npart;
                                                              old_i++, new_i++)
        SETV(Config(new_i), Config(old_i));

    return(new_cbod);
    }

/*-----------------------------------------------------------------------------
 *  mul_cbodies  --  form a new cbodies by multiplying an old cbodies by a
 *                   scalar.
 *                   accepts: old_cbod: pointer to many-cbody cbodies in
 *                                      standard newton0 form.
 *                               npart: the number of particles in that
 *                                      cbodies;
 *                        scale_factor: with which every entry in "old_cbod"
 *                                      multiplied.
 *                   returns: new_cbod: pointer to the new scaled copy.
 *-----------------------------------------------------------------------------
 */
cbodyptr  mul_cbodies(old_cbod, npart, scale_factor)
cbodyptr  old_cbod;
int  npart;
real  scale_factor;
    {
    cbodyptr  new_cbod;
    cbodyptr  new_i;      /* points to an individual particle of new_cbod */
    cbodyptr  old_i;      /* points to an individual particle of old_cbod */
    
    new_cbod = mk_cbodies(npart);
    
    for (old_i = old_cbod, new_i = new_cbod; old_i - old_cbod < npart;
                                                              old_i++, new_i++)
        MULVS(Config(new_i), Config(old_i), scale_factor);

    return(new_cbod);
    }

/*-----------------------------------------------------------------------------
 *  add_cbodies  --  adds two cbodiess of mbodies together
 *                   accepts: cbod1: pointer to cbodies in standard newton0
 *                                   form;
 *                            cbod2: pointer to second such cbodies;
 *                            npart: the number of particles in each cbodies.
 *                   returns: new_cbod: pointer to new cbodies which contains
 *                                     the sum of both old ones.
 *-----------------------------------------------------------------------------
 */
cbodyptr  add_cbodies(cbod1, cbod2, npart)
cbodyptr  cbod1;
cbodyptr  cbod2;
int  npart;
    {
    cbodyptr  new_cbod;
    cbodyptr  new_i, old1_i, old2_i;
    
    new_cbod = mk_cbodies(npart);
    
    for (new_i = new_cbod, old1_i= cbod1, old2_i = cbod2;
                         new_i - new_cbod < npart; new_i++, old1_i++, old2_i++)
        ADDV(Config(new_i), Config(old1_i), Config(old2_i));

    return(new_cbod);
    }

/*-----------------------------------------------------------------------------
 *  sub_cbodies  --  subtracts two cbodiess of mbodies
 *                   accepts: cbod1: pointer to cbodies in standard newton0
 *                                   form;
 *                            cbod2: pointer to second such cbodies;
 *                            npart: the number of particles in each cbodies.
 *                   returns: new_cbod: pointer to new cbodies which contains
 *                                      the difference of both old ones.
 *-----------------------------------------------------------------------------
 */
cbodyptr  sub_cbodies(cbod1, cbod2, npart)
cbodyptr  cbod1;
cbodyptr  cbod2;
int  npart;
    {
    cbodyptr  new_cbod;
    cbodyptr  new_i, old1_i, old2_i;
    
    new_cbod = mk_cbodies(npart);
    
    for (new_i = new_cbod, old1_i= cbod1, old2_i = cbod2;
                         new_i - new_cbod < npart; new_i++, old1_i++, old2_i++)
        SUBV(Config(new_i), Config(old1_i), Config(old2_i));

    return(new_cbod);
    }

/*-----------------------------------------------------------------------------
 *  scale_cbodies  --  multiplies an existing cbodies by a scalar.
 *                     accepts: old_cbod: pointer to cbodies in
 *                                       standard newton0 form.
 *                                 npart: the number of particles in that
 *                                        cbodies;
 *                        scale_factor: with which every entry in "old_cbod"
 *                                      multiplied.
 *-----------------------------------------------------------------------------
 */
void  scale_cbodies(cbod, npart, scale_factor)
cbodyptr  cbod;
int  npart;
real  scale_factor;
    {
    cbodyptr  cbody_i;      /* points to an individual particle of cbod */
    
    for (cbody_i = cbod; cbody_i - cbod < npart; cbody_i++)
        INCMULVS(Config(cbody_i), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_cbodies  --  increments a cbodies by adding an other cbodies
 *                   accepts: cbod: pointer to cbodies in standard newton0
 *                                   form;
 *                            d_cbod: pointer to the cbodies with which  cbod
 *                                   will be incremented.
 *                            npart: the number of particles in each cbodies.
 *-----------------------------------------------------------------------------
 */
void  inc_cbodies(cbod, d_cbod, npart)
cbodyptr  cbod;
cbodyptr  d_cbod;
int  npart;
    {
    cbodyptr  cbody_i, d_cbody_i;
    
    for (cbody_i= cbod, d_cbody_i = d_cbod; cbody_i - cbod < npart; 
                                                        cbody_i++, d_cbody_i++)
        INCADDV(Config(cbody_i), Config(d_cbody_i));
    }

/*-----------------------------------------------------------------------------
 *  dec_cbodies  --  decrements a cbodies by subtracting an other cbodies
 *                   accepts: cbod: pointer to cbodies in standard newton0
 *                                   form;
 *                            d_cbod: pointer to the cbodies with which  cbod
 *                                   will be decremented.
 *                            npart: the number of particles in each cbodies.
 *-----------------------------------------------------------------------------
 */
void  dec_cbodies(cbod, d_cbod, npart)
cbodyptr  cbod;
cbodyptr  d_cbod;
int  npart;
    {
    cbodyptr  cbody_i, d_cbody_i;
    
    for (cbody_i= cbod, d_cbody_i = d_cbod; cbody_i - cbod < npart; 
                                                        cbody_i++, d_cbody_i++)
        INCSUBV(Config(cbody_i), Config(d_cbody_i));
    }

/* endof: bodyalgebra.c */
