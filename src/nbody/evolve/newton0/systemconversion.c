/* systemconversion.c - */

/*
 *  systemconversion.c:  for converting different type of systems
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  "newton0.h"

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |   systemconversion.c :                                              |  */
/*  |                                                                     |  */
/*  |           This file contains procedures for operations on xsystems. |  */
/*  |                                                                     |  */
/*  |           xsystems stands for: systems, esystems, msystems,         |  */
/*  |                                psystems and csystems.               |  */
/*  |                                                                     |  */
/*  |           operations stands for: selecting, constructing            |  */
/*  |                                                                     |  */
/*  |           these operations are of mixed type, involving two         |  */
/*  |           basic types:                                              |  */
/*  |                                                                     |  */
/*  |           1) selection of smaller xsystems as separate entities     |  */
/*  |           extracted from ysystems (e.g. selecting accelerations     |  */
/*  |           from bodies);                                             |  */
/*  |                                                                     |  */
/*  |           2) construction of larger ysystems, either a) from        |  */
/*  |           smaller xsystems entirely, or b) from replacing smaller   |  */
/*  |           xsystems inside existing ysystems, or c) from adding      |  */
/*  |           smaller xsystems incrementally to existing ysystems       |  */
/*  |           (e.g. constructing psystems from positions and            |  */
/*  |           velocities and implementing an integration step           |  */
/*  |           by adding psystems to systems).                           |  */
/*  |                                                                     |  */
/*  |           summary: this file contains the following four            |  */
/*  |                    types of procedures:                             |  */
/*  |                                                                     |  */
/*  |                  PART I     :  selectors                            |  */
/*  |                                                                     |  */
/*  |                  PART II    :  constructors                         |  */
/*  |                    PART IIa :    pure constructors                  |  */
/*  |                    PART IIb :    renovators                         |  */
/*  |                    PART IIc :    annexers                           |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART I                                 |  */
/*  |                                                                     |  */
/*  |                  contains selection procedures                      |  */
/*  |                                                                     |  */
/*  |           the following selectors are offered:                      |  */
/*  |                                                                     |  */
/*  |           starting with a system, we can                            |  */
/*  |                select:          using:  returning a system of type: |  */
/*  |                                                                     |  */
/*  |             positions           sel_pos()       csystem             |  */
/*  |             velocities          sel_vel()       csystem             |  */
/*  |             accelerations       sel_acc()       csystem             |  */
/*  |             phase space part    sel_phase()     psystem             |  */
/*  |            (d/dt) phase space   sel_phasedot()  psystem             |  */
/*  |             its msyst part      sel_msys()      msystem             |  */
/*  |             its esyst part      sel_esys()      esystem             |  */
/*  |                                                                     |  */
/*  |           starting with a esystem, we can                           |  */
/*  |                select:          using:  returning a system of type: |  */
/*  |                                                                     |  */
/*  |             positions           sel_epos()      csystem             |  */
/*  |             velocities          sel_evel()      csystem             |  */
/*  |             phase space part    sel_ephase()    psystem             |  */
/*  |             its msyst part      sel_emsys()     msystem             |  */
/*  |                                                                     |  */
/*  |           starting with a msystem, we can                           |  */
/*  |                select:          using:  returning a system of type: |  */
/*  |                                                                     |  */
/*  |             positions           sel_mpos()      csystem             |  */
/*  |             velocities          sel_mvel()      csystem             |  */
/*  |             phase space part    sel_mphase()    psystem             |  */
/*  |                                                                     |  */
/*  |           starting with a psystem, we can                           |  */
/*  |                select:          using:  returning a system of type: |  */
/*  |                                                                     |  */
/*  |             positions           sel_ppos()      csystem             |  */
/*  |             velocities          sel_pvel()      csystem             |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  sel_emsys()  can be remembered       |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "select e's m part",                        |  */
/*  |                 just as you read  sel_mpos()  as                    |  */
/*  |                         "select m's positions".                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  sel_pos  --  take a newton0 system in standard form, and select the
 *               positions which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_pos(sys)
systptr  sys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(sys);
    CBodies(csys) = sel_b_pos(Bodies(sys), Nbody(sys));
    Nbody(csys) = Nbody(sys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_vel  --  take a newton0 system in standard form, and select the
 *               velocities which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_vel(sys)
systptr  sys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(sys);
    CBodies(csys) = sel_b_vel(Bodies(sys), Nbody(sys));
    Nbody(csys) = Nbody(sys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_acc  --  take a newton0 system in standard form, and select the
 *               accelerations which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_acc(sys)
systptr  sys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(sys);
    CBodies(csys) = sel_b_acc(Bodies(sys), Nbody(sys));
    Nbody(csys) = Nbody(sys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_phase  --  take a newton0 system in standard form, and select the
 *                 phase space part which is returned in the form of a psystem .
 *-----------------------------------------------------------------------------
 */
psystptr  sel_phase(sys)
systptr  sys;
    {
    psystptr  psys;

    psys = mk_empty_psystem();

    Tnow(psys) = Tnow(sys);
    PBodies(psys) = sel_b_phase(Bodies(sys), Nbody(sys));
    Nbody(psys) = Nbody(sys);

    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_phasedot  --  take a newton0 system in standard form, and select the
 *                    time derivative of the phase space part (velocities and
 *                    accelerations) which is returned in the form of a  pbody.
 *                    note: the exception for the REGULARIZATION option is not
 *                          very elegant; but it works for the moment.
 *-----------------------------------------------------------------------------
 */
psystptr  sel_phasedot(sys)
systptr  sys;
    {
    psystptr  psys;

    psys = mk_empty_psystem();     /* empty, since the next line creates the */
                                   /* bodies (this used to be a bug, when I  */
                                   /* used mk_psystem; repaired Oct.24 1987) */
    PBodies(psys) = sel_b_phasedot(Bodies(sys), Nbody(sys));
    Nbody(psys) = Nbody(sys);

#ifndef REGULARIZATION
    Tnow(psys) = 1.0;                                          /* dt/dt */
#else
    Tnow(psys) = 1.0 / Reglagrangian(sys);            /* pseudo time advance */

    MULVS(Com_Pos(psys), Com_Vel(sys), Tnow(psys));  /* c.o.m.: dr/dt = v*dt */
    CLRV(Com_Vel(psys));                     /* center-of-mass: dv/dt = 0    */

#endif

    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_msys  --  take a newton0 system in standard form, and select the
 *                msystem part which is returned in the form of a  msystem .
 *-----------------------------------------------------------------------------
 */
msystptr  sel_msys(sys)
systptr  sys;
    {
    msystptr  msys;

    msys = mk_empty_msystem();

    Tnow(msys) = Tnow(sys);
    MBodies(msys) = sel_b_msys(Bodies(sys), Nbody(sys));
    Nbody(msys) = Nbody(sys);

    return(msys);
    }

/*-----------------------------------------------------------------------------
 *  sel_esys  --  take a newton0 system in standard form, and select the
 *                esystem part which is returned in the form of a  esystem .
 *-----------------------------------------------------------------------------
 */
esystptr  sel_esys(sys)
systptr  sys;
    {
    esystptr  esys;

    esys = mk_empty_esystem();

    Tnow(esys) = Tnow(sys);
    EBodies(esys) = sel_b_esys(Bodies(sys), Nbody(sys));
    Nbody(esys) = Nbody(sys);

    return(esys);
    }

/*-----------------------------------------------------------------------------
 *  sel_epos  --  take a newton0 esystem in standard form, and select the
 *                positions which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_epos(esys)
esystptr  esys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(esys);
    CBodies(csys) = sel_b_epos(EBodies(esys), Nbody(esys));
    Nbody(csys) = Nbody(esys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_evel  --  take a newton0 esystem in standard form, and select the
 *                velocities which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_evel(esys)
esystptr  esys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(esys);
    CBodies(csys) = sel_b_evel(EBodies(esys), Nbody(esys));
    Nbody(csys) = Nbody(esys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_ephase  --  take a newton0 esystem in standard form, and select the
 *                  phase space part which is returned in the form of a  pbody.
 *-----------------------------------------------------------------------------
 */
psystptr  sel_ephase(esys)
esystptr  esys;
    {
    psystptr  psys;

    psys = mk_empty_psystem();

    Tnow(psys) = Tnow(esys);
    PBodies(psys) = sel_b_ephase(EBodies(esys), Nbody(esys));
    Nbody(psys) = Nbody(esys);

    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_emsys  --  take a newton0 esystem in standard form, and select the
 *                 msystem part which is returned in the form of a  msystem .
 *-----------------------------------------------------------------------------
 */
msystptr  sel_emsys(esys)
esystptr  esys;
    {
    msystptr  msys;

    msys = mk_empty_msystem();

    Tnow(msys) = Tnow(esys);
    MBodies(msys) = sel_b_emsys(EBodies(esys), Nbody(esys));
    Nbody(msys) = Nbody(esys);

    return(msys);
    }

/*-----------------------------------------------------------------------------
 *  sel_mpos  --  take a newton0 msystem in standard form, and select the
 *                positions which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_mpos(msys)
msystptr  msys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(msys);
    CBodies(csys) = sel_b_mpos(MBodies(msys), Nbody(msys));
    Nbody(csys) = Nbody(msys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_mvel  --  take a newton0 msystem in standard form, and select the
 *                velocities which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_mvel(msys)
msystptr  msys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(msys);
    CBodies(csys) = sel_b_mvel(MBodies(msys), Nbody(msys));
    Nbody(csys) = Nbody(msys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_mphase  --  take a newton0 msystem in standard form, and select the
 *                  phase space part which is returned in the form of a  psys.
 *-----------------------------------------------------------------------------
 */
psystptr  sel_mphase(msys)
msystptr  msys;
    {
    psystptr  psys;

    psys = mk_empty_psystem();

    Tnow(psys) = Tnow(msys);
    PBodies(psys) = sel_b_mphase(MBodies(msys), Nbody(msys));
    Nbody(psys) = Nbody(msys);

    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_ppos  --  take a newton0 psystem in standard form, and select the
 *                positions which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_ppos(psys)
psystptr  psys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(psys);
    CBodies(csys) = sel_b_ppos(PBodies(psys), Nbody(psys));
    Nbody(csys) = Nbody(psys);

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_pvel  --  take a newton0 psystem in standard form, and select the
 *                velocities which are returned in the form of a  csystem .
 *-----------------------------------------------------------------------------
 */
csystptr  sel_pvel(psys)
psystptr  psys;
    {
    csystptr  csys;

    csys = mk_empty_csystem();

    Tnow(csys) = Tnow(psys);
    CBodies(csys) = sel_b_pvel(PBodies(psys), Nbody(psys));
    Nbody(csys) = Nbody(psys);

    return(csys);
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                            PART II                                  |  */
/*  |                                                                     |  */
/*  |                  contains construction procedures                   |  */
/*  |                                                                     |  */
/*  |           the following pure constructors are provided,             |  */
/*  |           pure in the sense that a xbody and a ybody are            |  */
/*  |           composed to form a new zbody:                             |  */
/*  |                                                                     |  */
/*  |               constructors: these take a xbody and a ybody          |  */
/*  |                             and return the combination as           |  */
/*  |                             a zbody; only two are possible:         |  */
/*  |                                                                     |  */
/*  |               cons_pos_vel()   returns the combined   pbody         |  */
/*  |               cons_e_acc()     returns the combined   body          |  */
/*  |                                                                     |  */
/*  |           all other combinations of xsystems and ysystems           |  */
/*  |           lead to either an overdetermination or an                 |  */
/*  |           underdetermination in a resulting zbody;                  |  */
/*  |           two natural solutions are the following classes           |  */
/*  |           of impure constructors:                                   |  */
/*  |                                                                     |  */
/*  |               renovators: these take a xbody to replace             |  */
/*  |                           the corresponding xbody part of           |  */
/*  |                           a larger ybody;                           |  */
/*  |           and                                                       |  */
/*  |               annexers:   these take a xbody to increment           |  */
/*  |                           the corresponding xbody part of           |  */
/*  |                           a larger ybody.                           |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART IIa                                |  */
/*  |                                                                     |  */
/*  |                  contains pure constructors                         |  */
/*  |                                                                     |  */
/*  |           starting with two csystems, we can use                    |  */
/*  |                                                                     |  */
/*  |               cons_pos_vel()   to construct a   psystem;            |  */
/*  |                                                                     |  */
/*  |           starting with a cbody and a ebody, we can use             |  */
/*  |                                                                     |  */
/*  |               cons_e_acc()     to construct a    system;            |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  cons_pos_vel  --  constructs a phase space psystem from a position csystem
 *                    and a velocity csystem.
 *                    accepts: csys_pos: a csystem which contains the
 *                                       positions of all particles;
 *                             csys_vel: a csystem which contains the
 *                                       velocities of all particles;
 *                    returns: psys: the psystem resulting from the combination
 *                                   of the two csystems.
 *-----------------------------------------------------------------------------
 */
psystptr  cons_pos_vel(csys_pos, csys_vel)
csystptr  csys_pos;
csystptr  csys_vel;
    {
    psystptr  psys;
    
    psys = mk_empty_psystem();

    if (Nbody(csys_vel) == Nbody(csys_pos))
        Nbody(psys) = Nbody(csys_pos);
    else
	error("cons_pos_vel: the pos N = %d and vel N = %d differ\n", 
                                             Nbody(csys_pos), Nbody(csys_vel));
    if (coincide(Tnow(csys_pos), Tnow(csys_vel)))
	Tnow(psys) = Tnow(csys_pos);
    else
	error("cons_pos_vel: the pos time %g and vel time %g differ\n", 
                                               Tnow(csys_pos), Tnow(csys_vel));

    PBodies(psys) = cons_b_pos_vel(CBodies(csys_pos), CBodies(csys_vel),
                                                              Nbody(csys_pos));

    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  cons_e_acc  --  constructs a system from an acceleration csystem
 *                  and an esystem.
 *                  accepts: esys: a csystem which contains the complete
 *                                 evergy information of all particles;
 *                           csys: a csystem which contains the
 *                                 accelerations of all particles;
 *                  returns: sys: the system resulting from the combination
 *                                of the two csystems, in standard newton0 form
 *-----------------------------------------------------------------------------
 */
systptr  cons_e_acc(esys, csys)
esystptr  esys;
csystptr  csys;
    {
    systptr  sys;
    
    sys = mk_empty_system();
    
    if (Nbody(esys) == Nbody(csys))
        Nbody(sys) = Nbody(esys);
    else
	error("cons_e_acc: the esys N = %d and csys N = %d differ\n", 
                                                     Nbody(esys), Nbody(csys));
    if (coincide(Tnow(esys), Tnow(csys)))
	Tnow(sys) = Tnow(esys);
    else
	error("cons_e_acc: the esys time %g and csys time %g differ\n", 
                                                       Tnow(esys), Tnow(csys));

    Bodies(sys) = cons_b_e_acc(EBodies(esys), CBodies(csys), Nbody(esys));

    return(sys);
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART IIb                                |  */
/*  |                                                                     |  */
/*  |                       contains renovators                           |  */
/*  |                                                                     |  */
/*  |           the following renovators are offered:                     |  */
/*  |                                                                     |  */
/*  |           starting with a body, we can                              |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_pos()                       |  */
/*  |               velocities          renov_vel()                       |  */
/*  |               accelerations       renov_acc()                       |  */
/*  |               its pbody part      renov_phase()                     |  */
/*  |               its mbody part      renov_msys()                      |  */
/*  |               its ebody part      renov_esys()                      |  */
/*  |                                                                     |  */
/*  |           starting with a ebody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_epos()                      |  */
/*  |               velocities          renov_evel()                      |  */
/*  |               its pbody part      renov_ephase()                    |  */
/*  |               its mbody part      renov_emsys()                     |  */
/*  |                                                                     |  */
/*  |           starting with a mbody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_mpos()                      |  */
/*  |               velocities          renov_mvel()                      |  */
/*  |               its pbody part      renov_mphase()                    |  */
/*  |                                                                     |  */
/*  |           starting with a pbody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_ppos()                      |  */
/*  |               velocities          renov_pvel()                      |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  renov_emsys()  can be remembered     |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "renovate e's m part",                      |  */
/*  |                 just as you read  renov_mpos()  as                  |  */
/*  |                         "renovate m's positions".                   |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  renov_pos  --  renovates a system "sys" by replacing its positions by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_pos(sys, csys)
systptr  sys;
csystptr  csys;
    {
    renov_b_pos(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_vel  --  renovates a system "sys" by replacing its velocities by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_vel(sys, csys)
systptr  sys;
csystptr  csys;
    {
    renov_b_vel(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_acc  --  renovates a system "sys" by replacing its accelerations by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_acc(sys, csys)
systptr  sys;
csystptr  csys;
    {
    renov_b_acc(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_phase  --  renovates a system "sys" by replacing its phase space part
 *                   by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_phase(sys, psys)
systptr  sys;
psystptr  psys;
    {
    renov_b_phase(Bodies(sys), PBodies(psys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_msys  --  renovates a system "sys" by replacing its phase space part
 *                  and its masses by that found in "msys".
 *-----------------------------------------------------------------------------
 */
void  renov_msys(sys, msys)
systptr  sys;
msystptr  msys;
    {
    renov_b_msys(Bodies(sys), MBodies(msys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_esys  --  renovates a system "sys" by replacing its phase space part,
 *                  masses and potentials by that found in "esys".
 *-----------------------------------------------------------------------------
 */
void  renov_esys(sys, esys)
systptr  sys;
esystptr  esys;
    {
    renov_b_esys(Bodies(sys), EBodies(esys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  renov_epos  --  renovates a esystem "esys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_epos(esys, csys)
esystptr  esys;
csystptr  csys;
    {
    renov_b_epos(EBodies(esys), CBodies(csys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  renov_evel  --  renovates a esystem "esys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_evel(esys, csys)
esystptr  esys;
csystptr  csys;
    {
    renov_b_evel(EBodies(esys), CBodies(csys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  renov_ephase  --  renovates a esystem "esys" by replacing its phase space
 *                    part by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_ephase(esys, psys)
esystptr  esys;
psystptr  psys;
    {
    renov_b_ephase(EBodies(esys), PBodies(psys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  renov_emsys  --  renovates a esystem "esys" by replacing its phase space
 *                   part and its masses by that found in "msys".
 *-----------------------------------------------------------------------------
 */
void  renov_emsys(esys, msys)
esystptr  esys;
msystptr  msys;
    {
    renov_b_emsys(EBodies(esys), MBodies(msys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  renov_mpos  --  renovates a msystem "msys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_mpos(msys, csys)
msystptr  msys;
csystptr  csys;
    {
    renov_b_mpos(MBodies(msys), CBodies(csys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  renov_mvel  --  renovates a msystem "msys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_mvel(msys, csys)
msystptr  msys;
csystptr  csys;
    {
    renov_b_mvel(MBodies(msys), CBodies(csys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  renov_mphase  --  renovates a msystem "msys" by replacing its phase space
 *                    part by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_mphase(msys, psys)
msystptr  msys;
psystptr  psys;
    {
    renov_b_mphase(MBodies(msys), PBodies(psys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  renov_ppos  --  renovates a psystem "psys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_ppos(psys, csys)
psystptr  psys;
csystptr  csys;
    {
    renov_b_ppos(PBodies(psys), CBodies(csys), Nbody(psys));
    }

/*-----------------------------------------------------------------------------
 *  renov_pvel  --  renovates a psystem "psys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_pvel(psys, csys)
psystptr  psys;
csystptr  csys;
    {
    renov_b_pvel(PBodies(psys), CBodies(csys), Nbody(psys));
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART IIc                                |  */
/*  |                                                                     |  */
/*  |                        contains annexers                            |  */
/*  |                                                                     |  */
/*  |           the following annexers are offered:                       |  */
/*  |                                                                     |  */
/*  |           starting with a body, we can                              |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_pos()                       |  */
/*  |               velocities          annex_vel()                       |  */
/*  |               accelerations       annex_acc()                       |  */
/*  |               its pbody part      annex_phase()                     |  */
/*  |               its mbody part      annex_msys()                      |  */
/*  |               its ebody part      annex_esys()                      |  */
/*  |                                                                     |  */
/*  |           starting with a ebody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_epos()                      |  */
/*  |               velocities          annex_evel()                      |  */
/*  |               its pbody part      annex_ephase()                    |  */
/*  |               its mbody part      annex_emsys()                     |  */
/*  |                                                                     |  */
/*  |           starting with a mbody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_mpos()                      |  */
/*  |               velocities          annex_mvel()                      |  */
/*  |               its pbody part      annex_mphase()                    |  */
/*  |                                                                     |  */
/*  |           starting with a pbody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_ppos()                      |  */
/*  |               velocities          annex_pvel()                      |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  annex_emsys()  can be remembered     |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "annex e's m part",                         |  */
/*  |                 just as you read  annex_mpos()  as                  |  */
/*  |                         "annex m's positions".                      |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  annex_pos  --  annexes a csystem "csys" by adding its positions by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_pos(sys, csys)
systptr  sys;
csystptr  csys;
    {
    annex_b_pos(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_vel  --  annexes a csystem "csys" by adding its velocities by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_vel(sys, csys)
systptr  sys;
csystptr  csys;
    {
    annex_b_vel(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_acc  --  annexes a csystem "csys" by adding its accelerations by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_acc(sys, csys)
systptr  sys;
csystptr  csys;
    {
    annex_b_acc(Bodies(sys), CBodies(csys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_phase  --  annexes a psystem "psys" by adding its phase space part
 *                   by that found in the system "sys".
 *                   note: the exception for the REGULARIZATION option is not
 *                         very elegant; but it works for the moment.
 *-----------------------------------------------------------------------------
 */
void  annex_phase(sys, psys)
systptr  sys;
psystptr  psys;
    {
    Tnow(sys) += Tnow(psys);                               /*  t = t + dt  */

#ifdef REGULARIZATION
    INCADDV(Com_Pos(sys), Com_Pos(psys));
    INCADDV(Com_Vel(sys), Com_Vel(psys));
#endif

    annex_b_phase(Bodies(sys), PBodies(psys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_msys  --  annexes a msystem "msys" by adding its phase space part
 *                  and its masses by that found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_msys(sys, msys)
systptr  sys;
msystptr  msys;
    {
    annex_b_msys(Bodies(sys), MBodies(msys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_esys  --  annexes a esystem "esys" by adding its phase space part,
 *                  masses and potentials by that found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_esys(sys, esys)
systptr  sys;
esystptr  esys;
    {
    annex_b_esys(Bodies(sys), EBodies(esys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  annex_epos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_epos(esys, csys)
esystptr  esys;
csystptr  csys;
    {
    annex_b_epos(EBodies(esys), CBodies(csys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  annex_evel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_evel(esys, csys)
esystptr  esys;
csystptr  csys;
    {
    annex_b_evel(EBodies(esys), CBodies(csys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  annex_ephase  --  annexes a psystem "psys" by adding its phase space
 *                    part by that found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_ephase(esys, psys)
esystptr  esys;
psystptr  psys;
    {
    annex_b_ephase(EBodies(esys), PBodies(psys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  annex_emsys  --  annexes a msystem "msys" by adding its phase space
 *                   part and its masses by that found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_emsys(esys, msys)
esystptr  esys;
msystptr  msys;
    {
    annex_b_emsys(EBodies(esys), MBodies(msys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  annex_mpos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_mpos(msys, csys)
msystptr  msys;
csystptr  csys;
    {
    annex_b_mpos(MBodies(msys), CBodies(csys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  annex_mvel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_mvel(msys, csys)
msystptr  msys;
csystptr  csys;
    {
    annex_b_mvel(MBodies(msys), CBodies(csys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  annex_mphase  --  annexes a psystem "psys" by adding its phase space
 *                    part by that found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_mphase(msys, psys)
msystptr  msys;
psystptr  psys;
    {
    annex_b_mphase(MBodies(msys), PBodies(psys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  annex_ppos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the psystem "psys".
 *-----------------------------------------------------------------------------
 */
void  annex_ppos(psys, csys)
psystptr  psys;
csystptr  csys;
    {
    annex_b_ppos(PBodies(psys), CBodies(csys), Nbody(psys));
    }

/*-----------------------------------------------------------------------------
 *  annex_pvel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the psystem "psys".
 *-----------------------------------------------------------------------------
 */
void  annex_pvel(psys, csys)
psystptr  psys;
csystptr  csys;
    {
    annex_b_pvel(PBodies(psys), CBodies(csys), Nbody(psys));
    }

/* endof: systemconversion.c */
