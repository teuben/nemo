/*
 *  bodyconversion.c:  for converting different type of bodies into each other
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
#include  "newton0.h"

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |   bodyconversion.c :                                                |  */
/*  |                                                                     |  */
/*  |           This file contains procedures for operations on xbodies.  |  */
/*  |                                                                     |  */
/*  |           xbodies stands for: bodies, ebodies, mbodies, pbodies,    |  */
/*  |                               and cbodies.                          |  */
/*  |                                                                     |  */
/*  |           operations stands for: selecting, constructing            |  */
/*  |                                                                     |  */
/*  |           these operations are of mixed type, involving two         |  */
/*  |           basic types:                                              |  */
/*  |                                                                     |  */
/*  |           1) selection of smaller xbodies as separate entities      |  */
/*  |           extracted from ybodies (e.g. selecting accelerations      |  */
/*  |           from bodies);                                             |  */
/*  |                                                                     |  */
/*  |           2) construction of larger ybodies, either a) from         |  */
/*  |           smaller xbodies entirely, or b) from replacing smaller    |  */
/*  |           xbodies inside existing ybodies, or c) from adding        |  */
/*  |           smaller xbodies incrementally to existing ybodies         |  */
/*  |           (e.g. constructing pbodies from positions and             |  */
/*  |           velocities and implementing an integration step           |  */
/*  |           by adding pbodies to bodies).                             |  */
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
/*  |       NOTE: the terminology used has to be ironed out;              |  */
/*  |             for example, the _b_ in the function names is a         |  */
/*  |             simple (but inelegant) trick to avoid clashing          |  */
/*  |             with similar names of functions which convert           |  */
/*  |             systems rather than bodies (the higher level            |  */
/*  |             equivalent functions which reside in the file           |  */
/*  |             systemconversion.c).                                    |  */
/*  |             the reason for the inelegant naming is historical:      |  */
/*  |             at first bodies were introduced as the basic entities,  |  */
/*  |             and a group of bodies was called a system; later        |  */
/*  |             the current system idea was added as a higher-level     |  */
/*  |             abstraction. some of the names of descriptions in       |  */
/*  |             the present file will still use the old terminology;    |  */
/*  |             some day this will all be cleaned up and ironed out.    |  */
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
/*  |           starting with a body, we can                              |  */
/*  |                select:          using:   returning a body of type:  |  */
/*  |                                                                     |  */
/*  |             positions           sel_b_pos()       cbody             |  */
/*  |             velocities          sel_b_vel()       cbody             |  */
/*  |             accelerations       sel_b_acc()       cbody             |  */
/*  |             phase space part    sel_b_phase()     pbody             |  */
/*  |            (d/dt) phase space   sel_b_phasedot()  pbody             |  */
/*  |             its mbody part      sel_b_msys()      mbody             |  */
/*  |             its ebody part      sel_b_esys()      ebody             |  */
/*  |                                                                     |  */
/*  |           starting with a ebody, we can                             |  */
/*  |                select:          using:   returning a body of type:  |  */
/*  |                                                                     |  */
/*  |             positions           sel_b_epos()      cbody             |  */
/*  |             velocities          sel_b_evel()      cbody             |  */
/*  |             phase space part    sel_b_ephase()    pbody             |  */
/*  |             its mbody part      sel_b_emsys()     mbody             |  */
/*  |                                                                     |  */
/*  |           starting with a mbody, we can                             |  */
/*  |                select:          using:   returning a body of type:  |  */
/*  |                                                                     |  */
/*  |             positions           sel_b_mpos()      cbody             |  */
/*  |             velocities          sel_b_mvel()      cbody             |  */
/*  |             phase space part    sel_b_mphase()    pbody             |  */
/*  |                                                                     |  */
/*  |           starting with a pbody, we can                             |  */
/*  |                select:          using:   returning a body of type:  |  */
/*  |                                                                     |  */
/*  |             positions           sel_b_ppos()      cbody             |  */
/*  |             velocities          sel_b_pvel()      cbody             |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  sel_b_emsys()  can be remembered     |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "select e's m part",                        |  */
/*  |                 just as you read  sel_b_mpos()  as                  |  */
/*  |                         "select m's positions".                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  sel_b_pos  --  take a newton0 system in standard form, and select the
 *               positions which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_pos(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  csys;
    cbodyptr  pos_i;

    csys = mk_cbodies(npart);

    for (pos_i = csys, body_i = sys; body_i - sys < npart; pos_i++, body_i++) 
	SETV(Config(pos_i), Pos(body_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_vel  --  take a newton0 system in standard form, and select the
 *               velocities which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_vel(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  csys;
    cbodyptr  vel_i;

    csys = mk_cbodies(npart);

    for (vel_i = csys, body_i = sys; body_i - sys < npart; vel_i++, body_i++) 
	SETV(Config(vel_i), Vel(body_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_acc  --  take a newton0 system in standard form, and select the
 *               accelerations which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_acc(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  csys;
    cbodyptr  acc_i;

    csys = mk_cbodies(npart);

    for (acc_i = csys, body_i = sys; body_i - sys < npart; acc_i++, body_i++) 
	SETV(Config(acc_i), Acc(body_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_phase  --  take a newton0 system in standard form, and select the
 *                 phase space part which is returned in the form of a  pbody .
 *-----------------------------------------------------------------------------
 */
pbodyptr  sel_b_phase(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    pbodyptr  psys;
    pbodyptr  p_i;

    psys = mk_pbodies(npart);

    for (p_i = psys, body_i = sys; body_i - sys < npart; p_i++, body_i++) 
        {
	SETV(Pos(p_i), Pos(body_i));
	SETV(Vel(p_i), Vel(body_i));
        }
    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_phasedot  --  take a newton0 system in standard form, and select the
 *                    time derivative of the phase space part (velocities and
 *                    accelerations) which is returned in the form of a  pbody.
 *-----------------------------------------------------------------------------
 */
pbodyptr  sel_b_phasedot(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    pbodyptr  psys;
    pbodyptr  p_i;

    psys = mk_pbodies(npart);

    for (p_i = psys, body_i = sys; body_i - sys < npart; p_i++, body_i++) 
        {
#ifndef REGULARIZATION
	SETV(Pos(p_i), Vel(body_i));
	SETV(Vel(p_i), Acc(body_i));
#else
        SETV(PPos(p_i), PdPos_ds(body_i));
        SETV(PMom(p_i), PdMom_ds(body_i));
#endif
        }
    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_msys  --  take a newton0 system in standard form, and select the
 *                msystem part which is returned in the form of a  mbody .
 *-----------------------------------------------------------------------------
 */
mbodyptr  sel_b_msys(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    mbodyptr  msys;
    mbodyptr  m_i;

    msys = mk_mbodies(npart);

    for (m_i = msys, body_i = sys; body_i - sys < npart; m_i++, body_i++) 
        {
	SETV(Pos(m_i), Pos(body_i));
	SETV(Vel(m_i), Vel(body_i));
        Mass(m_i) = Mass(body_i);
        }
    return(msys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_esys  --  take a newton0 system in standard form, and select the
 *                esystem part which is returned in the form of a  ebody .
 *-----------------------------------------------------------------------------
 */
ebodyptr  sel_b_esys(sys, npart)
bodyptr  sys;
int  npart;
    {
    bodyptr  body_i;
    ebodyptr  esys;
    ebodyptr  e_i;

    esys = mk_ebodies(npart);

    for (e_i = esys, body_i = sys; body_i - sys < npart; e_i++, body_i++) 
        {
	SETV(Pos(e_i), Pos(body_i));
	SETV(Vel(e_i), Vel(body_i));
        Mass(e_i) = Mass(body_i);
        }
    return(esys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_epos  --  take a newton0 esystem in standard form, and select the
 *                positions which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_epos(esys, npart)
ebodyptr  esys;
int  npart;
    {
    ebodyptr  e_i;
    cbodyptr  csys;
    cbodyptr  pos_i;

    csys = mk_cbodies(npart);

    for (pos_i = csys, e_i = esys; e_i - esys < npart; pos_i++, e_i++) 
	SETV(Config(pos_i), Pos(e_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_evel  --  take a newton0 esystem in standard form, and select the
 *                velocities which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_evel(esys, npart)
ebodyptr  esys;
int  npart;
    {
    ebodyptr  e_i;
    cbodyptr  csys;
    cbodyptr  vel_i;

    csys = mk_cbodies(npart);

    for (vel_i = csys, e_i = esys; e_i - esys < npart; vel_i++, e_i++) 
	SETV(Config(vel_i), Vel(e_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_ephase  --  take a newton0 esystem in standard form, and select the
 *                  phase space part which is returned in the form of a  pbody.
 *-----------------------------------------------------------------------------
 */
pbodyptr  sel_b_ephase(esys, npart)
ebodyptr  esys;
int  npart;
    {
    ebodyptr  e_i;
    pbodyptr  psys;
    pbodyptr  p_i;

    psys = mk_pbodies(npart);

    for (p_i = psys, e_i = esys; e_i - esys < npart; p_i++, e_i++)
        {
	SETV(Pos(p_i), Pos(e_i));
	SETV(Vel(p_i), Vel(e_i));
        }
    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_emsys  --  take a newton0 esystem in standard form, and select the
 *                 msystem part which is returned in the form of a  mbody .
 *-----------------------------------------------------------------------------
 */
mbodyptr  sel_b_emsys(esys, npart)
ebodyptr  esys;
int  npart;
    {
    ebodyptr  e_i;
    mbodyptr  msys;
    mbodyptr  m_i;

    msys = mk_mbodies(npart);

    for (m_i = msys, e_i = esys; e_i - esys < npart; m_i++, e_i++)
        {
	SETV(Pos(m_i), Pos(e_i));
	SETV(Vel(m_i), Vel(e_i));
        Mass(m_i) = Mass(e_i);
        }
    return(msys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_mpos  --  take a newton0 msystem in standard form, and select the
 *                positions which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_mpos(msys, npart)
mbodyptr  msys;
int  npart;
    {
    mbodyptr  m_i;
    cbodyptr  csys;
    cbodyptr  pos_i;

    csys = mk_cbodies(npart);

    for (pos_i = csys, m_i = msys; m_i - msys < npart; pos_i++, m_i++)
	SETV(Config(pos_i), Pos(m_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_mvel  --  take a newton0 msystem in standard form, and select the
 *                velocities which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_mvel(msys, npart)
mbodyptr  msys;
int  npart;
    {
    mbodyptr  m_i;
    cbodyptr  csys;
    cbodyptr  vel_i;

    csys = mk_cbodies(npart);

    for (vel_i = csys, m_i = msys; m_i - msys < npart; vel_i++, m_i++) 
	SETV(Config(vel_i), Vel(m_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_mphase  --  take a newton0 msystem in standard form, and select the
 *                  phase space part which is returned in the form of a  pbody.
 *-----------------------------------------------------------------------------
 */
pbodyptr  sel_b_mphase(msys, npart)
mbodyptr  msys;
int  npart;
    {
    mbodyptr  m_i;
    pbodyptr  psys;
    pbodyptr  p_i;

    psys = mk_pbodies(npart);

    for (p_i = psys, m_i = msys; m_i - msys < npart; p_i++, m_i++)
        {
	SETV(Pos(p_i), Pos(m_i));
	SETV(Vel(p_i), Vel(m_i));
        }
    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_ppos  --  take a newton0 psystem in standard form, and select the
 *                positions which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_ppos(psys, npart)
pbodyptr  psys;
int  npart;
    {
    pbodyptr  p_i;
    cbodyptr  csys;
    cbodyptr  pos_i;

    csys = mk_cbodies(npart);

    for (pos_i = csys, p_i = psys; p_i - psys < npart; pos_i++, p_i++) 
	SETV(Config(pos_i), Pos(p_i));

    return(csys);
    }

/*-----------------------------------------------------------------------------
 *  sel_b_pvel  --  take a newton0 psystem in standard form, and select the
 *                velocities which are returned in the form of a  cbody .
 *-----------------------------------------------------------------------------
 */
cbodyptr  sel_b_pvel(psys, npart)
pbodyptr  psys;
int  npart;
    {
    pbodyptr  p_i;
    cbodyptr  csys;
    cbodyptr  vel_i;

    csys = mk_cbodies(npart);

    for (vel_i = csys, p_i = psys; p_i - psys < npart; vel_i++, p_i++) 
	SETV(Config(vel_i), Vel(p_i));

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
/*  |               cons_b_pos_vel()   returns the combined   pbody       |  */
/*  |               cons_b_e_acc()     returns the combined   body        |  */
/*  |                                                                     |  */
/*  |           all other combinations of xbodies and ybodies             |  */
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
/*  |           starting with two cbodies, we can use                     |  */
/*  |                                                                     |  */
/*  |               cons_b_pos_vel()   to construct a   pbody;            |  */
/*  |                                                                     |  */
/*  |           starting with a cbody and a ebody, we can use             |  */
/*  |                                                                     |  */
/*  |               cons_b_e_acc()     to construct a    body;            |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  cons_b_pos_vel  --  constructs a phase space psystem from a position 
 *                      csystem and a velocity csystem.
 *                      accepts: csys_pos: a csystem which contains the
 *                                         positions of all particles;
 *                               csys_vel: a csystem which contains the
 *                                         velocities of all particles;
 *                      returns: psys: the psystem resulting from the 
 *                                     combination of the two csystems.
 *-----------------------------------------------------------------------------
 */
pbodyptr  cons_b_pos_vel(csys_pos, csys_vel, npart)
cbodyptr  csys_pos;
cbodyptr  csys_vel;
int  npart;
    {
    cbodyptr  pos_i, vel_i;
    pbodyptr  psys;
    pbodyptr  pbody_i;
    
    psys = mk_pbodies(npart);
    
    for (pbody_i = psys, pos_i= csys_pos, vel_i = csys_vel;
                           pbody_i - psys < npart; pbody_i++, pos_i++, vel_i++)
        {
        SETV(Pos(pbody_i), Config(pos_i));
        SETV(Vel(pbody_i), Config(vel_i));
	}    
    return(psys);
    }

/*-----------------------------------------------------------------------------
 *  cons_b_e_acc  --  constructs a system from an acceleration csystem
 *                  and an esystem.
 *                  accepts: esys: a csystem which contains the complete
 *                                 evergy information of all particles;
 *                           csys: a csystem which contains the
 *                                 accelerations of all particles;
 *                  returns: sys: the system resulting from the combination
 *                                of the two csystems, in standard newton0 form
 *-----------------------------------------------------------------------------
 */
bodyptr  cons_b_e_acc(esys, csys, npart)
ebodyptr  esys;
cbodyptr  csys;
int  npart;
    {
    ebodyptr  ebody_i;
    cbodyptr  acc_i;
    bodyptr  sys;
    bodyptr  body_i;
    
    sys = mk_bodies(npart);
    
    for (body_i = sys, ebody_i= esys, acc_i = csys;
                            body_i - sys < npart; body_i++, ebody_i++, acc_i++)
        {
        SETV(Pos(body_i), Pos(ebody_i));
        SETV(Vel(body_i), Vel(ebody_i));
        Mass(body_i) = Mass(ebody_i);
        Pot(body_i) = Pot(ebody_i);
        SETV(Acc(body_i), Config(acc_i));
	}    
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
/*  |               positions           renov_b_pos()                     |  */
/*  |               velocities          renov_b_vel()                     |  */
/*  |               accelerations       renov_b_acc()                     |  */
/*  |               its pbody part      renov_b_phase()                   |  */
/*  |               its mbody part      renov_b_msys()                    |  */
/*  |               its ebody part      renov_b_esys()                    |  */
/*  |                                                                     |  */
/*  |           starting with a ebody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_b_epos()                    |  */
/*  |               velocities          renov_b_evel()                    |  */
/*  |               its pbody part      renov_b_ephase()                  |  */
/*  |               its mbody part      renov_b_emsys()                   |  */
/*  |                                                                     |  */
/*  |           starting with a mbody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_b_mpos()                    |  */
/*  |               velocities          renov_b_mvel()                    |  */
/*  |               its pbody part      renov_b_mphase()                  |  */
/*  |                                                                     |  */
/*  |           starting with a pbody, we can                             |  */
/*  |                 renovate:           using:                          |  */
/*  |                                                                     |  */
/*  |               positions           renov_b_ppos()                    |  */
/*  |               velocities          renov_b_pvel()                    |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  renov_b_emsys()  can be remembered   |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "renovate e's m part",                      |  */
/*  |                 just as you read  renov_b_mpos()  as                |  */
/*  |                         "renovate m's positions".                   |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  renov_b_pos  --  renovates a system "sys" by replacing its positions by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_pos(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, body_i = sys; body_i - sys < npart; pos_i++, body_i++) 
	SETV(Pos(body_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_vel  --  renovates a system "sys" by replacing its velocities by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_vel(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, body_i = sys; body_i - sys < npart; vel_i++, body_i++) 
	SETV(Vel(body_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_acc  --  renovates a system "sys" by replacing its accelerations by
 *                 those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_acc(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  acc_i;
    
    for (acc_i = csys, body_i = sys; body_i - sys < npart; acc_i++, body_i++) 
	SETV(Acc(body_i), Config(acc_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_phase  --  renovates a system "sys" by replacing its phase space 
 *                     part by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_phase(sys, psys, npart)
bodyptr  sys;
pbodyptr  psys;
int  npart;
    {
    bodyptr  body_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, body_i = sys; body_i-sys < npart; pbody_i++, body_i++)
        {
	SETV(Pos(body_i), Pos(pbody_i));
	SETV(Vel(body_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_msys  --  renovates a system "sys" by replacing its phase space 
 *                    part and its masses by that found in "msys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_msys(sys, msys, npart)
bodyptr  sys;
mbodyptr  msys;
int  npart;
    {
    bodyptr  body_i;
    mbodyptr  mbody_i;
    
    for (mbody_i = msys, body_i = sys; body_i-sys < npart; mbody_i++, body_i++)
        {
	SETV(Pos(body_i), Pos(mbody_i));
	SETV(Vel(body_i), Vel(mbody_i));
        Mass(body_i) = Mass(mbody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_esys  --  renovates a system "sys" by replacing its phase space part,
 *                  masses and potentials by that found in "esys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_esys(sys, esys, npart)
bodyptr  sys;
ebodyptr  esys;
int  npart;
    {
    bodyptr  body_i;
    ebodyptr  ebody_i;
    
    for (ebody_i = esys, body_i = sys; body_i-sys < npart; ebody_i++, body_i++)
        {
	SETV(Pos(body_i), Pos(ebody_i));
	SETV(Vel(body_i), Vel(ebody_i));
        Mass(body_i) = Mass(ebody_i);
        Pot(body_i) = Pot(ebody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_epos  --  renovates a esystem "esys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_epos(esys, csys, npart)
ebodyptr  esys;
cbodyptr  csys;
int  npart;
    {
    ebodyptr  ebody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, ebody_i = esys; ebody_i - esys < npart; 
                                                            pos_i++, ebody_i++)
	SETV(Pos(ebody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_evel  --  renovates a esystem "esys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_evel(esys, csys, npart)
ebodyptr  esys;
cbodyptr  csys;
int  npart;
    {
    ebodyptr  ebody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, ebody_i = esys; ebody_i - esys < npart;
                                                            vel_i++, ebody_i++)
	SETV(Vel(ebody_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_ephase  --  renovates a esystem "esys" by replacing its phase space
 *                    part by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_ephase(esys, psys, npart)
ebodyptr  esys;
pbodyptr  psys;
int  npart;
    {
    ebodyptr  ebody_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, ebody_i = esys; ebody_i-esys < npart;
                                                          pbody_i++, ebody_i++)
        {
	SETV(Pos(ebody_i), Pos(pbody_i));
	SETV(Vel(ebody_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_emsys  --  renovates a esystem "esys" by replacing its phase space
 *                   part and its masses by that found in "msys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_emsys(esys, msys, npart)
ebodyptr  esys;
mbodyptr  msys;
int  npart;
    {
    ebodyptr  ebody_i;
    mbodyptr  mbody_i;
    
    for (mbody_i = msys, ebody_i = esys; ebody_i-esys < npart;
                                                          mbody_i++, ebody_i++)
        {
	SETV(Pos(ebody_i), Pos(mbody_i));
	SETV(Vel(ebody_i), Vel(mbody_i));
        Mass(ebody_i) = Mass(mbody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_mpos  --  renovates a msystem "msys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_mpos(msys, csys, npart)
mbodyptr  msys;
cbodyptr  csys;
int  npart;
    {
    mbodyptr  mbody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, mbody_i = msys; mbody_i - msys < npart; 
                                                            pos_i++, mbody_i++)
	SETV(Pos(mbody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_mvel  --  renovates a msystem "msys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_mvel(msys, csys, npart)
mbodyptr  msys;
cbodyptr  csys;
int  npart;
    {
    mbodyptr  mbody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, mbody_i = msys; mbody_i - msys < npart;
                                                            vel_i++, mbody_i++)
	SETV(Vel(mbody_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_mphase  --  renovates a msystem "msys" by replacing its phase space
 *                    part by that found in "psys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_mphase(msys, psys, npart)
mbodyptr  msys;
pbodyptr  psys;
int  npart;
    {
    mbodyptr  mbody_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, mbody_i = msys; mbody_i-msys < npart;
                                                          pbody_i++, mbody_i++)
        {
	SETV(Pos(mbody_i), Pos(pbody_i));
	SETV(Vel(mbody_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  renov_b_ppos  --  renovates a psystem "psys" by replacing its positions by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_ppos(psys, csys, npart)
pbodyptr  psys;
cbodyptr  csys;
int  npart;
    {
    pbodyptr  pbody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, pbody_i = psys; pbody_i - psys < npart; 
                                                            pos_i++, pbody_i++)
	SETV(Pos(pbody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  renov_b_pvel  --  renovates a psystem "psys" by replacing its velocities by
 *                  those found in "csys".
 *-----------------------------------------------------------------------------
 */
void  renov_b_pvel(psys, csys, npart)
pbodyptr  psys;
cbodyptr  csys;
int  npart;
    {
    pbodyptr  pbody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, pbody_i = psys; pbody_i - psys < npart;
                                                            vel_i++, pbody_i++)
	SETV(Vel(pbody_i), Config(vel_i));
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
/*  |               positions           annex_b_pos()                     |  */
/*  |               velocities          annex_b_vel()                     |  */
/*  |               accelerations       annex_b_acc()                     |  */
/*  |               its pbody part      annex_b_phase()                   |  */
/*  |               its mbody part      annex_b_msys()                    |  */
/*  |               its ebody part      annex_b_esys()                    |  */
/*  |                                                                     |  */
/*  |           starting with a ebody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_b_epos()                    |  */
/*  |               velocities          annex_b_evel()                    |  */
/*  |               its pbody part      annex_b_ephase()                  |  */
/*  |               its mbody part      annex_b_emsys()                   |  */
/*  |                                                                     |  */
/*  |           starting with a mbody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_b_mpos()                    |  */
/*  |               velocities          annex_b_mvel()                    |  */
/*  |               its pbody part      annex_b_mphase()                  |  */
/*  |                                                                     |  */
/*  |           starting with a pbody, we can                             |  */
/*  |                  annex:             using:                          |  */
/*  |                                                                     |  */
/*  |               positions           annex_b_ppos()                    |  */
/*  |               velocities          annex_b_pvel()                    |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |           note: the notations  annex_b_emsys()  can be remembered   |  */
/*  |                 most easily by reading them as                      |  */
/*  |                         "annex e's m part",                         |  */
/*  |                 just as you read  annex_b_mpos()  as                |  */
/*  |                         "annex m's positions".                      |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  annex_b_pos  --  annexes a csystem "csys" by adding its positions by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_pos(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, body_i = sys; body_i - sys < npart; pos_i++, body_i++) 
	INCADDV(Pos(body_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_vel  --  annexes a csystem "csys" by adding its velocities by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_vel(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, body_i = sys; body_i - sys < npart; vel_i++, body_i++) 
	INCADDV(Vel(body_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_acc  --  annexes a csystem "csys" by adding its accelerations by
 *                 those found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_acc(sys, csys, npart)
bodyptr  sys;
cbodyptr  csys;
int  npart;
    {
    bodyptr  body_i;
    cbodyptr  acc_i;
    
    for (acc_i = csys, body_i = sys; body_i - sys < npart; acc_i++, body_i++) 
	INCADDV(Acc(body_i), Config(acc_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_phase  --  annexes a psystem "psys" by adding its phase space part
 *                   by that found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_phase(sys, psys, npart)
bodyptr  sys;
pbodyptr  psys;
int  npart;
    {
    bodyptr  body_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, body_i = sys; body_i-sys < npart; pbody_i++, body_i++)
        {
	INCADDV(Pos(body_i), Pos(pbody_i));
	INCADDV(Vel(body_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_msys  --  annexes a msystem "msys" by adding its phase space part
 *                  and its masses by that found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_msys(sys, msys, npart)
bodyptr  sys;
mbodyptr  msys;
int  npart;
    {
    bodyptr  body_i;
    mbodyptr  mbody_i;
    
    for (mbody_i = msys, body_i = sys; body_i-sys < npart; mbody_i++, body_i++)
        {
	INCADDV(Pos(body_i), Pos(mbody_i));
	INCADDV(Vel(body_i), Vel(mbody_i));
        Mass(body_i) += Mass(mbody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_esys  --  annexes a esystem "esys" by adding its phase space part,
 *                  masses and potentials by that found in the system "sys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_esys(sys, esys, npart)
bodyptr  sys;
ebodyptr  esys;
int  npart;
    {
    bodyptr  body_i;
    ebodyptr  ebody_i;
    
    for (ebody_i = esys, body_i = sys; body_i-sys < npart; ebody_i++, body_i++)
        {
	INCADDV(Pos(body_i), Pos(ebody_i));
	INCADDV(Vel(body_i), Vel(ebody_i));
        Mass(body_i) += Mass(ebody_i);
        Pot(body_i) += Pot(ebody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_epos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_epos(esys, csys, npart)
ebodyptr  esys;
cbodyptr  csys;
int  npart;
    {
    ebodyptr  ebody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, ebody_i = esys; ebody_i - esys < npart; 
                                                            pos_i++, ebody_i++)
	INCADDV(Pos(ebody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_evel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_evel(esys, csys, npart)
ebodyptr  esys;
cbodyptr  csys;
int  npart;
    {
    ebodyptr  ebody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, ebody_i = esys; ebody_i - esys < npart;
                                                            vel_i++, ebody_i++)
	INCADDV(Vel(ebody_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_ephase  --  annexes a psystem "psys" by adding its phase space
 *                    part by that found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_ephase(esys, psys, npart)
ebodyptr  esys;
pbodyptr  psys;
int  npart;
    {
    ebodyptr  ebody_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, ebody_i = esys; ebody_i-esys < npart;
                                                          pbody_i++, ebody_i++)
        {
	INCADDV(Pos(ebody_i), Pos(pbody_i));
	INCADDV(Vel(ebody_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_emsys  --  annexes a msystem "msys" by adding its phase space
 *                   part and its masses by that found in the esystem "esys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_emsys(esys, msys, npart)
ebodyptr  esys;
mbodyptr  msys;
int  npart;
    {
    ebodyptr  ebody_i;
    mbodyptr  mbody_i;
    
    for (mbody_i = msys, ebody_i = esys; ebody_i-esys < npart;
                                                          mbody_i++, ebody_i++)
        {
	INCADDV(Pos(ebody_i), Pos(mbody_i));
	INCADDV(Vel(ebody_i), Vel(mbody_i));
        Mass(ebody_i) += Mass(mbody_i);
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_mpos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_mpos(msys, csys, npart)
mbodyptr  msys;
cbodyptr  csys;
int  npart;
    {
    mbodyptr  mbody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, mbody_i = msys; mbody_i - msys < npart; 
                                                            pos_i++, mbody_i++)
	INCADDV(Pos(mbody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_mvel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_mvel(msys, csys, npart)
mbodyptr  msys;
cbodyptr  csys;
int  npart;
    {
    mbodyptr  mbody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, mbody_i = msys; mbody_i - msys < npart;
                                                            vel_i++, mbody_i++)
	INCADDV(Vel(mbody_i), Config(vel_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_mphase  --  annexes a psystem "psys" by adding its phase space
 *                    part by that found in the msystem "msys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_mphase(msys, psys, npart)
mbodyptr  msys;
pbodyptr  psys;
int  npart;
    {
    mbodyptr  mbody_i;
    pbodyptr  pbody_i;
    
    for (pbody_i = psys, mbody_i = msys; mbody_i-msys < npart;
                                                          pbody_i++, mbody_i++)
        {
	INCADDV(Pos(mbody_i), Pos(pbody_i));
	INCADDV(Vel(mbody_i), Vel(pbody_i));
        }
    }

/*-----------------------------------------------------------------------------
 *  annex_b_ppos  --  annexes a csystem "csys" by adding its positions by
 *                  those found in the psystem "psys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_ppos(psys, csys, npart)
pbodyptr  psys;
cbodyptr  csys;
int  npart;
    {
    pbodyptr  pbody_i;
    cbodyptr  pos_i;
    
    for (pos_i = csys, pbody_i = psys; pbody_i - psys < npart; 
                                                            pos_i++, pbody_i++)
	INCADDV(Pos(pbody_i), Config(pos_i));
    }

/*-----------------------------------------------------------------------------
 *  annex_b_pvel  --  annexes a csystem "csys" by adding its velocities by
 *                  those found in the psystem "psys".
 *-----------------------------------------------------------------------------
 */
void  annex_b_pvel(psys, csys, npart)
pbodyptr  psys;
cbodyptr  csys;
int  npart;
    {
    pbodyptr  pbody_i;
    cbodyptr  vel_i;
    
    for (vel_i = csys, pbody_i = psys; pbody_i - psys < npart;
                                                            vel_i++, pbody_i++)
	INCADDV(Vel(pbody_i), Config(vel_i));
    }

/* endof: bodyconversion.c */
