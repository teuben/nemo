/* systemalgebra.c - */

/*
 *  systemalgebra.c:  for algebraic operations on systems in newton0
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */

#include  "newton0.h"

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |   systemalgebra.c :                                                 |  */
/*  |                                                                     |  */
/*  |           This file contains procedures for operations on xsystems. |  */
/*  |                                                                     |  */
/*  |           xsystems stands for: systems, esystems, msystems,         |  */
/*  |                                psystems, and csystems.              |  */
/*  |                                                                     |  */
/*  |           operations stands for: creating, clearing, copying,       |  */
/*  |                                  multiplying, adding, subtracting,  |  */
/*  |                                  scaling, incrementing and          |  */
/*  |                                  decrementing.                      |  */
/*  |                                                                     |  */
/*  |           these operations preserve the type of xsystem;            |  */
/*  |           other operations which mix different types of x           |  */
/*  |           reside in the file  systemconversion.c .                  |  */
/*  |                                                                     |  */
/*  |           meaning of operations:                                    |  */
/*  |                                                                     |  */
/*  |           multiplying: take a xsystem, multiply all components      |  */
/*  |                        with a constant scale factor, and            |  */
/*  |                        assign the resulting xsystem to a newly      |  */
/*  |                        created xsystem whose pointer is returned;   |  */
/*  |                        the old system remains unchanged.            |  */
/*  |           scaling: as multiplying, but now the results are          |  */
/*  |                    stored in the old system, and nothing is         |  */
/*  |                    returned explicitly by the procedure.            |  */
/*  |           adding: add two xsystems (sub)component-wise, assigning   |  */
/*  |                   the result to a new xsystem a pointer to          |  */
/*  |                   which is returned; the old xsystem remains        |  */
/*  |                   unchanged.                                        |  */
/*  |           incrementing: as adding, but now the results are          |  */
/*  |                         stored in the first of the two systems;     |  */
/*  |                         the other system remains unchanged,         |  */
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
/*  |     TABLE OF CONTENTS:                                              |  */
/*  |                                                                     |  */
/*  |                              PART I                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on systems:            |  */
/*  |                                                                     |  */
/*  |               systptr  mk_empty_system(       )                     |  */
/*  |               systptr  mk_system(              npart )              |  */
/*  |               void     rm_system(        sys  )                     |  */
/*  |               void     clear_system(     sys  )                     |  */
/*  |               systptr  cp_system(        sys  )                     |  */
/*  |               systptr  mul_system(       sys , factor )             |  */
/*  |               systptr  add_system( sys1, sys2 )                     |  */
/*  |               systptr  sub_system( sys1, sys2 )                     |  */
/*  |               void     scale_system(     sys , factor )             |  */
/*  |               void     inc_system( sys ,dsys  )                     |  */
/*  |               void     dec_system( sys ,dsys  )                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART II                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on esystems:           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_esystem(         )                 |  */
/*  |               esystptr  mk_esystem(                npart )          |  */
/*  |               void      rm_eystem(           sys  )                 |  */
/*  |               void      clear_esystem(      esys  )                 |  */
/*  |               esystptr  cp_esystem(         esys  )                 |  */
/*  |               esystptr  mul_esystem(        esys , factor )         |  */
/*  |               esystptr  add_esystem( esys1, esys2 )                 |  */
/*  |               esystptr  sub_esystem( esys1, esys2 )                 |  */
/*  |               void      scale_esystem(      esys , factor )         |  */
/*  |               void      inc_esystem( esys, desys  )                 |  */
/*  |               void      dec_esystem( esys, desys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                             PART III                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on msystems:           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_msystem(         )                 |  */
/*  |               msystptr  mk_msystem(                npart )          |  */
/*  |               void      rm_mystem(           sys  )                 |  */
/*  |               void      clear_msystem(      msys  )                 |  */
/*  |               msystptr  cp_msystem(         msys  )                 |  */
/*  |               msystptr  mul_msystem(        msys , factor )         |  */
/*  |               msystptr  add_msystem( msys1, msys2 )                 |  */
/*  |               msystptr  sub_msystem( msys1, msys2 )                 |  */
/*  |               void      scale_msystem(      msys , factor )         |  */
/*  |               void      inc_msystem( msys, dmsys  )                 |  */
/*  |               void      dec_msystem( msys, dmsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                   (continued)       |  */
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
/*  |           contains procedures for operations on psystems:           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_psystem(         )                 |  */
/*  |               psystptr  mk_psystem(                npart )          |  */
/*  |               void      rm_pystem(           sys  )                 |  */
/*  |               void      clear_psystem(      psys  )                 |  */
/*  |               psystptr  cp_psystem(         psys  )                 |  */
/*  |               psystptr  mul_psystem(        psys , factor )         |  */
/*  |               psystptr  add_psystem( psys1, psys2 )                 |  */
/*  |               psystptr  sub_psystem( psys1, psys2 )                 |  */
/*  |               void      scale_psystem(      psys , factor )         |  */
/*  |               void      inc_psystem( psys, dpsys  )                 |  */
/*  |               void      dec_psystem( psys, dpsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART V                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on csystems:           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_csystem(         )                 |  */
/*  |               csystptr  mk_csystem(                npart )          |  */
/*  |               void      rm_cystem(           sys  )                 |  */
/*  |               void      clear_csystem(      csys  )                 |  */
/*  |               csystptr  cp_csystem(         csys  )                 |  */
/*  |               csystptr  mul_csystem(        csys , factor )         |  */
/*  |               csystptr  add_csystem( csys1, csys2 )                 |  */
/*  |               csystptr  sub_csystem( csys1, csys2 )                 |  */
/*  |               void      scale_csystem(      csys , factor )         |  */
/*  |               void      inc_csystem( csys, dcsys  )                 |  */
/*  |               void      dec_csystem( csys, dcsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_empty_system  --  allocates memory for a many-body system in standard 
 *                       newton0 form, but does not allocate memory for the
 *                       bodies; the number of particles is therefore left
 *                       unspecified.
 *                       returns: new_sys: a pointer to the new system.
 *-----------------------------------------------------------------------------
 */
systptr  mk_empty_system()
    {
    systptr  new_sys;

    new_sys = (systptr) malloc(sizeof(syst));
    if (new_sys == NULL)
	error("mk_empty_system: not enough memory left for a new system\n");

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  mk_system  --  allocates memory for a many-body system in standard newton0
 *                 form, and properly initializes the types of all bodies.
 *                 returns: new_sys: a pointer to the new system.
 *-----------------------------------------------------------------------------
 */
systptr  mk_system(npart)
int  npart;
    {
    systptr  new_sys;

    new_sys = (systptr) malloc(sizeof(syst));
    if (new_sys == NULL)
	error("mk_system: not enough memory left for a new system\n");

    Nbody(new_sys) = npart;
    Bodies(new_sys) = mk_bodies(npart);

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  rm_system  --  deallocates memory for a system, after deallocating
 *                 memory for the bodies.
 *                 accept: old_sys: pointer to system which is to be removed.
 *                 note: the regularized mass matrix contents are not
 *                       automatically removed, since copies of regularized
 *                       systems often point to a shared massmatrix.
 *                 note: the tree root note is not removed -- THIS SHOULD BE
 *                       IMPLEMENTED before using  rm_system()  in a tree
 *                       calculation.
 *-----------------------------------------------------------------------------
 */
void  rm_system(old_sys)
systptr  old_sys;
    {
    free(Bodies(old_sys));
    free(old_sys);
    }

/*-----------------------------------------------------------------------------
 *  clear_system  --  sets all values to zero in a system in standard newton0
 *                    form.
 *                    accepts: sys: pointer to a many-body system
 *-----------------------------------------------------------------------------
 */
void  clear_system(sys)
systptr  sys;
    {
    Tnow(sys) = 0.0;
    
    clear_bodies(Bodies(sys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  cp_system  --  makes a copy of a system in standard newton0 form.
 *                 accepts: old_sys: pointer to a many-body system;
 *                 returns: new_sys: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
systptr  cp_system(old_sys)
systptr  old_sys;
    {
    systptr  new_sys;
    
    new_sys = mk_empty_system();
    Nbody(new_sys) = Nbody(old_sys);
    Tnow(new_sys) = Tnow(old_sys);

    Bodies(new_sys) = cp_bodies(Bodies(old_sys), Nbody(old_sys));    
#ifdef TREE
    Ncell(new_sys) = Ncell(old_sys);
    Root(new_sys) = Root(old_sys);
    Nmaxcell(new_sys) = Nmaxcell(old_sys);
    SETV(Potmincorner(new_sys), Potmincorner(old_sys));
    Potsizesq(new_sys) = Potsizesq(old_sys);
#endif
#ifdef REGULARIZATION
    SETV(Com_Pos(new_sys), Com_Pos(old_sys));
    SETV(Com_Vel(new_sys), Com_Vel(old_sys));
    Massmatrix(new_sys)= Massmatrix(old_sys); /* points to the same old data */
    Regenergy(new_sys) =  Regenergy(old_sys);
    Reglagrangian(new_sys) = Reglagrangian(old_sys);
    Reghamiltonian(new_sys) = Reghamiltonian(old_sys);
#endif

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  mul_system  --  form a new system by multiplying an old system by a scalar.
 *                  accepts: old_sys: pointer to a many-body system in standard
 *                                    newton0 form.
 *                      scale_factor: with which every entry in "old_sys"
 *                                    multiplied.
 *                  returns: new_sys: pointer to the new scaled copy.
 *                  note: EVERYTHING EXCEPT the particle number is scaled,
 *                        even the masses, potentials and accelerations of the
 *                        bodies, as well as the time; if this is not what you 
 *                        want, then use system subset representations such as
 *                        psystem or csystem.
 *-----------------------------------------------------------------------------
 */
systptr  mul_system(old_sys, scale_factor)
systptr  old_sys;
real  scale_factor;
    {
    systptr  new_sys;
    
    new_sys = mk_empty_system();
    Nbody(new_sys) = Nbody(old_sys);
    Tnow(new_sys) = scale_factor * Tnow(old_sys);

    Bodies(new_sys) = mul_bodies(Bodies(old_sys), Nbody(old_sys),scale_factor);

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  add_system  --  adds two systems of bodies together
 *                  accepts: sys1: pointer to a system in standard newton0 form
 *                           sys2: pointer to a second such system;
 *                  returns: new_sys: pointer to a new system which contains
 *                                    the sum of both old ones.
 *                  note: EVERYTHING EXCEPT the particle number is added,
 *                        even the masses, potentials, accelerations and time;
 *                        if this is not what you want, then either put 
 *                        these equal to zero in one of the systems
 *                        so that the sum inherits the values of the other
 *                        system, or use system subset representations such as
 *                        psystem or csystem.
 *-----------------------------------------------------------------------------
 */
systptr  add_system(sys1, sys2)
systptr  sys1;
systptr  sys2;
    {
    systptr  new_sys;
    
    new_sys = mk_empty_system();
    Nbody(new_sys) = Nbody(sys1);
    Tnow(new_sys) = Tnow(sys1) + Tnow(sys2);

    Bodies(new_sys) = add_bodies(Bodies(sys1), Bodies(sys2), Nbody(new_sys));

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  sub_system  --  subtracts two systems of bodies
 *                  accepts: sys1: pointer to a system in standard newton0 form
 *                           sys2: pointer to a second such system;
 *                  returns: new_sys: pointer to a new system which contains
 *                                    the difference of both old ones.
 *                  note: EVERYTHING EXCEPT the particle number is subtracted,
 *                        even the masses, potentials, accelerations and time;
 *                        if this is not what you want, then either put 
 *                        these equal to zero in one of the systems so that 
 *                        the difference inherits the values of the other
 *                        system (modulo minus sign), or use system subset
 *                        representations such as psystem or csystem.
 *            disclaimer: the management is not responsible for such repulsive
 *                        creations as are formed by subtracting large masses
 *                        from small ones ...
 *-----------------------------------------------------------------------------
 */
systptr  sub_system(sys1, sys2)
systptr  sys1;
systptr  sys2;
    {
    systptr  new_sys;
    
    new_sys = mk_empty_system();
    Nbody(new_sys) = Nbody(sys1);
    Tnow(new_sys) = Tnow(sys1) - Tnow(sys2);

    Bodies(new_sys) = sub_bodies(Bodies(sys1), Bodies(sys2), Nbody(new_sys));

    return(new_sys);
    }

/*-----------------------------------------------------------------------------
 *  scale_system  --  multiplies an existing system by a scalar.
 *                    accepts: old_sys: pointer to a many-body system in
 *                                      standard newton0 form.
 *                        scale_factor: with which every entry in "old_sys"
 *                                      multiplied.
 *                    note: EVERYTHING EXCEPT the particle number is scaled,
 *                          even the masses, potentials and accelerations of
 *                          the bodies, as well as the time; if this is not
 *                          what you want, then use system subset
 *                          representations such as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  scale_system(sys, scale_factor)
systptr  sys;
real  scale_factor;
    {
    Tnow(sys) *= scale_factor;

    scale_bodies(Bodies(sys), Nbody(sys), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_system  --  increments a system by adding an other system
 *                  accepts: sys: pointer to a system in standard newton0 form
 *                           d_sys: pointer to the system with which  sys
 *                                 will be incremented.
 *                  note: EVERYTHING EXCEPT the particle number is added,
 *                        even the masses, potentials, accelerations and time;
 *                        if this is not what you want, then either put 
 *                        these equal to zero in one of the systems
 *                        so that the sum inherits the values of the other
 *                        system, or use system subset representations such as
 *                        psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  inc_system(sys, d_sys)
systptr  sys;
systptr  d_sys;
    {
    Tnow(sys) += Tnow(d_sys);

    inc_bodies(Bodies(sys), Bodies(d_sys), Nbody(sys));
    }

/*-----------------------------------------------------------------------------
 *  dec_system  --  decrements a system by subtracting an other system
 *                  accepts: sys: pointer to a system in standard newton0 form
 *                           d_sys: pointer to the system with which  sys
 *                                 will be decremented.
 *                  note: EVERYTHING EXCEPT the particle number is subtracted,
 *                        even the masses, potentials, accelerations and time;
 *                        if this is not what you want, then either put 
 *                        these equal to zero in one of the systems so that 
 *                        the difference inherits the values of the other
 *                        system (modulo minus sign), or use system subset
 *                        representations such as psystem or csystem.
 *            disclaimer: the management is not responsible for such repulsive
 *                        creations as are formed by subtracting large masses
 *                        from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_system(sys, d_sys)
systptr  sys;
systptr  d_sys;
    {
    Tnow(sys) -= Tnow(d_sys);

    dec_bodies(Bodies(sys), Bodies(d_sys), Nbody(sys));
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART II                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on esystems:           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_esystem(         )                 |  */
/*  |               esystptr  mk_esystem(                npart )          |  */
/*  |               void      rm_eystem(           sys  )                 |  */
/*  |               void      clear_esystem(      esys  )                 |  */
/*  |               esystptr  cp_esystem(         esys  )                 |  */
/*  |               esystptr  mul_esystem(        esys , factor )         |  */
/*  |               esystptr  add_esystem( esys1, esys2 )                 |  */
/*  |               esystptr  sub_esystem( esys1, esys2 )                 |  */
/*  |               void      scale_esystem(      esys , factor )         |  */
/*  |               void      inc_esystem( esys, desys  )                 |  */
/*  |               void      dec_esystem( esys, desys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_empty_esystem  --  allocates memory for a many-body esystem in standard 
 *                        newton0 form, but does not allocate memory for the
 *                        ebodies; the number of particles is therefore left
 *                        unspecified.
 *                        returns: new_esys: a pointer to the new esystem.
 *-----------------------------------------------------------------------------
 */
esystptr  mk_empty_esystem()
    {
    esystptr  new_esys;

    new_esys = (esystptr) malloc(sizeof(esyst));
    if (new_esys == NULL)
	error("mk_empty_esystem: not enough memory left for a new esystem\n");

    return(new_esys);
    }

/*-----------------------------------------------------------------------------
 *  mk_esystem  --  allocates memory for a many-body esystem in standard
 *                  newton0 form, and properly initializes the types of
 *                  all ebodies.
 *                  returns: new_esys: a pointer to the new esystem.
 *-----------------------------------------------------------------------------
 */
esystptr  mk_esystem(npart)
int  npart;
    {
    esystptr  new_esys;

    new_esys = (esystptr) malloc(sizeof(esyst));
    if (new_esys == NULL)
	error("mk_esystem: not enough memory left for a new esystem\n");

    Nbody(new_esys) = npart;
    EBodies(new_esys) = mk_ebodies(npart);

    return(new_esys);
    }

/*-----------------------------------------------------------------------------
 *  rm_esystem  --  deallocates memory for a esystem, after deallocating
 *                  memory for the ebodies.
 *                  accept: old_esys: pointer to esystem which is to be removed
 *-----------------------------------------------------------------------------
 */
void  rm_esystem(old_esys)
esystptr  old_esys;
    {
    free(EBodies(old_esys));
    free(old_esys);
    }

/*-----------------------------------------------------------------------------
 *  clear_esystem  --  sets all values to zero in a esystem in standard newton0
 *                     form.
 *                     accepts: esys: pointer to a many-body esystem
 *-----------------------------------------------------------------------------
 */
void  clear_esystem(esys)
esystptr  esys;
    {
    Tnow(esys) = 0.0;
    
    clear_ebodies(EBodies(esys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  cp_esystem  --  makes a copy of a esystem in standard newton0 form.
 *                 accepts: old_esys: pointer to a many-body esystem;
 *                 returns: new_esys: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
esystptr  cp_esystem(old_esys)
esystptr  old_esys;
    {
    esystptr  new_esys;
    
    new_esys = mk_empty_esystem();
    Nbody(new_esys) = Nbody(old_esys);
    Tnow(new_esys) = Tnow(old_esys);

    EBodies(new_esys) = cp_ebodies(EBodies(old_esys), Nbody(old_esys));    

    return(new_esys);
    }
  
/*-----------------------------------------------------------------------------
 *  mul_esystem  --  form a new esystem by multiplying an old esystem by
 *                   a scalar.
 *                   accepts: old_esys: pointer to a many-body esystem in
 *                                      standard newton0 form.
 *                        scale_factor: with which every entry in "old_esys"
 *                                      multiplied.
 *                   returns: new_esys: pointer to the new scaled copy.
 *                   note: EVERYTHING EXCEPT the particle number is scaled,
 *                         even the masses and potentials of the ebodies,
 *                         as well as the time; if this is not what you want
 *                         then use system subset representations such as
 *                         psystem or csystem.
 *-----------------------------------------------------------------------------
 */
esystptr  mul_esystem(old_esys, scale_factor)
esystptr  old_esys;
real  scale_factor;
    {
    esystptr  new_esys;
    
    new_esys = mk_empty_esystem();
    Nbody(new_esys) = Nbody(old_esys);
    Tnow(new_esys) = scale_factor * Tnow(old_esys);

    EBodies(new_esys) = mul_ebodies(EBodies(old_esys), Nbody(old_esys),
                                                                 scale_factor);

    return(new_esys);
    }

/*-----------------------------------------------------------------------------
 *  add_esystem  --  adds two esystems of ebodies together
 *                   accepts: esys1: pointer to a esystem in standard newton0
 *                                   form;
 *                            esys2: pointer to a second such esystem.
 *                   returns: new_esys: pointer to a new esystem which contains
 *                                     the sum of both old ones.
 *                   note: EVERYTHING EXCEPT the particle number is added,
 *                         even the masses, potentials, accelerations and time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the esystems
 *                         so that the sum inherits the values of the other
 *                         esystem, or use system subset representations such
 *                         as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
esystptr  add_esystem(esys1, esys2)
esystptr  esys1;
esystptr  esys2;
    {
    esystptr  new_esys;
    
    new_esys = mk_empty_esystem();
    Nbody(new_esys) = Nbody(esys1);
    Tnow(new_esys) = Tnow(esys1) + Tnow(esys2);

    EBodies(new_esys) = add_ebodies(EBodies(esys1), EBodies(esys2),
                                                              Nbody(new_esys));

    return(new_esys);
    }

/*-----------------------------------------------------------------------------
 *  sub_esystem  --  subtracts two esystems of ebodies
 *                   accepts: esys1: pointer to a esystem in standard newton0
 *                                   form;
 *                            esys2: pointer to a second such esystem.
 *                   returns: new_esys: pointer to a new esystem which contains
 *                                     the difference of both old ones.
 *                   note: EVERYTHING EXCEPT the particle number is subtracted,
 *                         even the masses, potentials, accelerations and time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the esystems so that 
 *                         the difference inherits the values of the other
 *                         esystem (modulo minus sign), or use system subset
 *                         representations such as psystem or csystem.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
esystptr  sub_esystem(esys1, esys2)
esystptr  esys1;
esystptr  esys2;
    {
    esystptr  new_esys;
    
    new_esys = mk_empty_esystem();
    Nbody(new_esys) = Nbody(esys1);
    Tnow(new_esys) = Tnow(esys1) - Tnow(esys2);

    EBodies(new_esys) = sub_ebodies(EBodies(esys1), EBodies(esys2),
                                                              Nbody(new_esys));

    return(new_esys);
    }

/*-----------------------------------------------------------------------------
 *  scale_esystem  --  multiplies an existing esystem by a scalar.
 *                     accepts: old_esys: pointer to a many-body esystem in
 *                                       standard newton0 form.
 *                         scale_factor: with which every entry in "old_esys"
 *                                       multiplied.
 *                     note: EVERYTHING EXCEPT the particle number is scaled,
 *                           even the masses, potentials and accelerations of
 *                           the ebodies, as well as the time; if this is not
 *                           what you want, then use system subset
 *                           representations such as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  scale_esystem(esys, scale_factor)
esystptr  esys;
real  scale_factor;
    {
    Tnow(esys) *= scale_factor;

    scale_ebodies(EBodies(esys), Nbody(esys), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_esystem  --  increments a esystem by adding an other esystem
 *                   accepts: esys: pointer to a esystem in standard newton0
 *                                  form
 *                          d_esys: pointer to the esystem with which  esys
 *                                  will be incremented.
 *                   note: EVERYTHING EXCEPT the particle number is added,
 *                         even the masses, potentials, accelerations and time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the esystems
 *                         so that the sum inherits the values of the other
 *                         esystem, or use system subset representations such
 *                         as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  inc_esystem(esys, d_esys)
esystptr  esys;
esystptr  d_esys;
    {
    Tnow(esys) += Tnow(d_esys);

    inc_ebodies(EBodies(esys), EBodies(d_esys), Nbody(esys));
    }

/*-----------------------------------------------------------------------------
 *  dec_esystem  --  decrements a esystem by subtracting an other esystem
 *                   accepts: esys: pointer to a esystem in standard newton0
 *                                  form
 *                            d_esys: pointer to the esystem with which  esys
 *                                  will be decremented.
 *                   note: EVERYTHING EXCEPT the particle number is subtracted,
 *                         even the masses, potentials, accelerations and time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the esystems so that 
 *                         the difference inherits the values of the other
 *                         esystem (modulo minus sign), or use system subset
 *                         representations such as psystem or csystem.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_esystem(esys, d_esys)
esystptr  esys;
esystptr  d_esys;
    {
    Tnow(esys) -= Tnow(d_esys);

    dec_ebodies(EBodies(esys), EBodies(d_esys), Nbody(esys));
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART III                               |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on msystems.           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_msystem(         )                 |  */
/*  |               msystptr  mk_msystem(                npart )          |  */
/*  |               void      rm_mystem(           sys  )                 |  */
/*  |               void      clear_msystem(      msys  )                 |  */
/*  |               msystptr  cp_msystem(         msys  )                 |  */
/*  |               msystptr  mul_msystem(        msys , factor )         |  */
/*  |               msystptr  add_msystem( msys1, msys2 )                 |  */
/*  |               msystptr  sub_msystem( msys1, msys2 )                 |  */
/*  |               void      scale_msystem(      msys , factor )         |  */
/*  |               void      inc_msystem( msys, dmsys  )                 |  */
/*  |               void      dec_msystem( msys, dmsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_empty_msystem  --  allocates memory for a many-body msystem in standard 
 *                        newton0 form, but does not allocate memory for the
 *                        mbodies; the number of particles is therefore left
 *                        unspecified.
 *                        returns: new_msys: a pointer to the new msystem.
 *-----------------------------------------------------------------------------
 */
msystptr  mk_empty_msystem()
    {
    msystptr  new_msys;

    new_msys = (msystptr) malloc(sizeof(msyst));
    if (new_msys == NULL)
	error("mk_empty_msystem: not enough memory left for a new msystem\n");

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  mk_msystem  --  allocates memory for a many-body msystem in standard
 *                  newton0 form, and properly initializes the types of
 *                  all mbodies.
 *                  returns: new_msys: a pointer to the new msystem.
 *-----------------------------------------------------------------------------
 */
msystptr  mk_msystem(npart)
int  npart;
    {
    msystptr  new_msys;

    new_msys = (msystptr) malloc(sizeof(msyst));
    if (new_msys == NULL)
	error("mk_msystem: not enough memory left for a new msystem\n");

    Nbody(new_msys) = npart;
    MBodies(new_msys) = mk_mbodies(npart);

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  rm_msystem  --  deallocates memory for a msystem, after deallocating
 *                  memory for the mbodies.
 *                  accept: old_msys: pointer to msystem which is to be removed
 *-----------------------------------------------------------------------------
 */
void  rm_msystem(old_msys)
msystptr  old_msys;
    {
    free(MBodies(old_msys));
    free(old_msys);
    }

/*-----------------------------------------------------------------------------
 *  clear_msystem  --  sets all values to zero in a msystem in standard newton0
 *                     form.
 *                     accepts: msys: pointer to a many-body msystem
 *-----------------------------------------------------------------------------
 */
void  clear_msystem(msys)
msystptr  msys;
    {
    Tnow(msys) = 0.0;
    
    clear_mbodies(MBodies(msys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  cp_msystem  --  makes a copy of a msystem in standard newton0 form.
 *                 accepts: old_msys: pointer to a many-body msystem;
 *                 returns: new_msys: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
msystptr  cp_msystem(old_msys)
msystptr  old_msys;
    {
    msystptr  new_msys;
    
    new_msys = mk_empty_msystem();
    Nbody(new_msys) = Nbody(old_msys);
    Tnow(new_msys) = Tnow(old_msys);

    MBodies(new_msys) = cp_mbodies(MBodies(old_msys), Nbody(old_msys));    

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  mul_msystem  --  form a new msystem by multiplying an old msystem by
 *                   a scalar.
 *                   accepts: old_msys: pointer to a many-body msystem in
 *                                      standard newton0 form.
 *                        scale_factor: with which every entry in "old_msys"
 *                                      multiplied.
 *                   returns: new_msys: pointer to the new scaled copy.
 *                   note: EVERYTHING EXCEPT the particle number is scaled,
 *                         masses as well as the time; if this is not what you
 *                         want then use system subset representations such as
 *                         psystem or csystem.
 *-----------------------------------------------------------------------------
 */
msystptr  mul_msystem(old_msys, scale_factor)
msystptr  old_msys;
real  scale_factor;
    {
    msystptr  new_msys;
    
    new_msys = mk_empty_msystem();
    Nbody(new_msys) = Nbody(old_msys);
    Tnow(new_msys) = scale_factor * Tnow(old_msys);

    MBodies(new_msys) = mul_mbodies(MBodies(old_msys), Nbody(old_msys),
                                                                 scale_factor);

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  add_msystem  --  adds two msystems of mbodies together
 *                   accepts: msys1: pointer to a msystem in standard newton0
 *                                   form;
 *                            msys2: pointer to a second such msystem.
 *                   returns: new_msys: pointer to a new msystem which contains
 *                                     the sum of both old ones.
 *                   note: EVERYTHING EXCEPT the particle number is added,
 *                         even the masses and the time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the msystems
 *                         so that the sum inherits the values of the other
 *                         msystem, or use system subset representations such
 *                         as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
msystptr  add_msystem(msys1, msys2)
msystptr  msys1;
msystptr  msys2;
    {
    msystptr  new_msys;
    
    new_msys = mk_empty_msystem();
    Nbody(new_msys) = Nbody(msys1);
    Tnow(new_msys) = Tnow(msys1) + Tnow(msys2);

    MBodies(new_msys) = add_mbodies(MBodies(msys1), MBodies(msys2),
                                                              Nbody(new_msys));

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  sub_msystem  --  subtracts two msystems of mbodies
 *                   accepts: msys1: pointer to a msystem in standard newton0
 *                                   form;
 *                            msys2: pointer to a second such msystem.
 *                   returns: new_msys: pointer to a new msystem which contains
 *                                     the difference of both old ones.
 *                   note: EVERYTHING EXCEPT the particle number is subtracted,
 *                         even the masses and the time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the msystems so that 
 *                         the difference inherits the values of the other
 *                         msystem (modulo minus sign), or use system subset
 *                         representations such as psystem or csystem.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
msystptr  sub_msystem(msys1, msys2)
msystptr  msys1;
msystptr  msys2;
    {
    msystptr  new_msys;
    
    new_msys = mk_empty_msystem();
    Nbody(new_msys) = Nbody(msys1);
    Tnow(new_msys) = Tnow(msys1) - Tnow(msys2);

    MBodies(new_msys) = sub_mbodies(MBodies(msys1), MBodies(msys2),
                                                              Nbody(new_msys));

    return(new_msys);
    }

/*-----------------------------------------------------------------------------
 *  scale_msystem  --  multiplies an existing msystem by a scalar.
 *                     accepts: old_msys: pointer to a many-body msystem in
 *                                       standard newton0 form.
 *                         scale_factor: with which every entry in "old_msys"
 *                                       multiplied.
 *                     note: EVERYTHING EXCEPT the particle number is scaled,
 *                           even the masses as well as the time; if this is
 *                           not what you want, then use system subset
 *                           representations such as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  scale_msystem(msys, scale_factor)
msystptr  msys;
real  scale_factor;
    {
    Tnow(msys) *= scale_factor;

    scale_mbodies(MBodies(msys), Nbody(msys), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_msystem  --  increments a msystem by adding an other msystem
 *                   accepts: msys: pointer to a msystem in standard newton0
 *                                  form
 *                          d_msys: pointer to the msystem with which  msys
 *                                  will be incremented.
 *                   note: EVERYTHING EXCEPT the particle number is added,
 *                         even the masses and the time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the msystems
 *                         so that the sum inherits the values of the other
 *                         msystem, or use system subset representations such
 *                         as psystem or csystem.
 *-----------------------------------------------------------------------------
 */
void  inc_msystem(msys, d_msys)
msystptr  msys;
msystptr  d_msys;
    {
    Tnow(msys) += Tnow(d_msys);

    inc_mbodies(MBodies(msys), MBodies(d_msys), Nbody(msys));
    }

/*-----------------------------------------------------------------------------
 *  dec_msystem  --  decrements a msystem by subtracting an other msystem
 *                   accepts: msys: pointer to a msystem in standard newton0
 *                                  form
 *                            d_msys: pointer to the msystem with which  msys
 *                                  will be decremented.
 *                   note: EVERYTHING EXCEPT the particle number is subtracted,
 *                         even the masses and the time;
 *                         if this is not what you want, then either put 
 *                         these equal to zero in one of the msystems so that 
 *                         the difference inherits the values of the other
 *                         msystem (modulo minus sign), or use system subset
 *                         representations such as psystem or csystem.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_msystem(msys, d_msys)
msystptr  msys;
msystptr  d_msys;
    {
    Tnow(msys) -= Tnow(d_msys);

    dec_mbodies(MBodies(msys), MBodies(d_msys), Nbody(msys));
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART IV                                |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on psystems.           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_psystem(         )                 |  */
/*  |               psystptr  mk_psystem(                npart )          |  */
/*  |               void      rm_pystem(           sys  )                 |  */
/*  |               void      clear_psystem(      psys  )                 |  */
/*  |               psystptr  cp_psystem(         psys  )                 |  */
/*  |               psystptr  mul_psystem(        psys , factor )         |  */
/*  |               psystptr  add_psystem( psys1, psys2 )                 |  */
/*  |               psystptr  sub_psystem( psys1, psys2 )                 |  */
/*  |               void      scale_psystem(      psys , factor )         |  */
/*  |               void      inc_psystem( psys, dpsys  )                 |  */
/*  |               void      dec_psystem( psys, dpsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_empty_psystem  --  allocates memory for a many-body psystem in standard 
 *                        newton0 form, but does not allocate memory for the
 *                        pbodies; the number of particles is therefore left
 *                        unspecified.
 *                        returns: new_psys: a pointer to the new psystem.
 *-----------------------------------------------------------------------------
 */
psystptr  mk_empty_psystem()
    {
    psystptr  new_psys;

    new_psys = (psystptr) malloc(sizeof(psyst));
    if (new_psys == NULL)
	error("mk_empty_psystem: not enough memory left for a new psystem\n");

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  mk_psystem  --  allocates memory for a many-body psystem in standard
 *                  newton0 form, and properly initializes the types of
 *                  all pbodies.
 *                  returns: new_psys: a pointer to the new psystem.
 *-----------------------------------------------------------------------------
 */
psystptr  mk_psystem(npart)
int  npart;
    {
    psystptr  new_psys;

    new_psys = (psystptr) malloc(sizeof(psyst));
    if (new_psys == NULL)
	error("mk_psystem: not enough memory left for a new psystem\n");

    Nbody(new_psys) = npart;
    PBodies(new_psys) = mk_pbodies(npart);

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  rm_psystem  --  deallocates memory for a psystem, after deallocating
 *                  memory for the pbodies.
 *                  accept: old_psys: pointer to psystem which is to be removed
 *-----------------------------------------------------------------------------
 */
void  rm_psystem(old_psys)
psystptr  old_psys;
    {
    free(PBodies(old_psys));
    free(old_psys);
    }

/*-----------------------------------------------------------------------------
 *  clear_psystem  --  sets all values to zero in a psystem in standard newton0
 *                     form.
 *                     accepts: psys: pointer to a many-body psystem
 *-----------------------------------------------------------------------------
 */
void  clear_psystem(psys)
psystptr  psys;
    {
    Tnow(psys) = 0.0;
    
    clear_pbodies(PBodies(psys), Nbody(psys));
    }

/*-----------------------------------------------------------------------------
 *  cp_psystem  --  makes a copy of a psystem in standard newton0 form.
 *                 accepts: old_psys: pointer to a many-body psystem;
 *                 returns: new_psys: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
psystptr  cp_psystem(old_psys)
psystptr  old_psys;
    {
    psystptr  new_psys;
    
    new_psys = mk_empty_psystem();
    Nbody(new_psys) = Nbody(old_psys);
    Tnow(new_psys) = Tnow(old_psys);

    PBodies(new_psys) = cp_pbodies(PBodies(old_psys), Nbody(old_psys));    

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  mul_psystem  --  form a new psystem by multiplying an old psystem by
 *                   a scalar.
 *                   accepts: old_psys: pointer to a many-body psystem in
 *                                      standard newton0 form.
 *                        scale_factor: with which every entry in "old_psys"
 *                                      multiplied.
 *                   returns: new_psys: pointer to the new scaled copy.
 *                   note: the exception for the REGULARIZATION option is not
 *                         very elegant; but it works for the moment.
 *-----------------------------------------------------------------------------
 */
psystptr  mul_psystem(old_psys, scale_factor)
psystptr  old_psys;
real  scale_factor;
    {
    psystptr  new_psys;
    
    new_psys = mk_empty_psystem();
    Nbody(new_psys) = Nbody(old_psys);
    Tnow(new_psys) = scale_factor * Tnow(old_psys);

    PBodies(new_psys) = mul_pbodies(PBodies(old_psys), Nbody(old_psys),
                                                                 scale_factor);

#ifdef REGULARIZATION
    MULVS(Com_Pos(new_psys), Com_Pos(old_psys), scale_factor);
    MULVS(Com_Vel(new_psys), Com_Vel(old_psys), scale_factor);
#endif

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  add_psystem  --  adds two psystems of pbodies together
 *                   accepts: psys1: pointer to a psystem in standard newton0
 *                                   form;
 *                            psys2: pointer to a second such psystem.
 *                   returns: new_psys: pointer to a new psystem which contains
 *                                     the sum of both old ones.
 *                     note: the exception for the REGULARIZATION option is not
 *                           very elegant; but it works for the moment.
 *-----------------------------------------------------------------------------
 */
psystptr  add_psystem(psys1, psys2)
psystptr  psys1;
psystptr  psys2;
    {
    psystptr  new_psys;
    
    new_psys = mk_empty_psystem();
    Nbody(new_psys) = Nbody(psys1);
    Tnow(new_psys) = Tnow(psys1) + Tnow(psys2);

    PBodies(new_psys) = add_pbodies(PBodies(psys1), PBodies(psys2),
                                                              Nbody(new_psys));

#ifdef REGULARIZATION
    ADDV(Com_Pos(new_psys), Com_Pos(psys1), Com_Pos(psys2));
    ADDV(Com_Vel(new_psys), Com_Vel(psys1), Com_Vel(psys2));
#endif

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  sub_psystem  --  subtracts two psystems of pbodies
 *                   accepts: psys1: pointer to a psystem in standard newton0
 *                                   form;
 *                            psys2: pointer to a second such psystem.
 *                   returns: new_psys: pointer to a new psystem which contains
 *                                     the difference of both old ones.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
psystptr  sub_psystem(psys1, psys2)
psystptr  psys1;
psystptr  psys2;
    {
    psystptr  new_psys;
    
    new_psys = mk_empty_psystem();
    Nbody(new_psys) = Nbody(psys1);
    Tnow(new_psys) = Tnow(psys1) - Tnow(psys2);

    PBodies(new_psys) = sub_pbodies(PBodies(psys1), PBodies(psys2),
                                                              Nbody(new_psys));

    return(new_psys);
    }

/*-----------------------------------------------------------------------------
 *  scale_psystem  --  multiplies an existing psystem by a scalar.
 *                     accepts: old_psys: pointer to a many-body psystem in
 *                                       standard newton0 form.
 *                         scale_factor: with which every entry in "old_psys"
 *                                       multiplied.
 *                     note: the exception for the REGULARIZATION option is not
 *                           very elegant; but it works for the moment.
 *-----------------------------------------------------------------------------
 */
void  scale_psystem(psys, scale_factor)
psystptr  psys;
real  scale_factor;
    {
    Tnow(psys) *= scale_factor;

#ifdef REGULARIZATION
    INCMULVS(Com_Pos(psys), scale_factor);
    INCMULVS(Com_Vel(psys), scale_factor);
#endif

    scale_pbodies(PBodies(psys), Nbody(psys), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_psystem  --  increments a psystem by adding an other psystem
 *                   accepts: psys: pointer to a psystem in standard newton0
 *                                  form
 *                          d_psys: pointer to the psystem with which  psys
 *                                  will be incremented.
 *-----------------------------------------------------------------------------
 */
void  inc_psystem(psys, d_psys)
psystptr  psys;
psystptr  d_psys;
    {
    Tnow(psys) += Tnow(d_psys);

    inc_pbodies(PBodies(psys), PBodies(d_psys), Nbody(psys));
    }

/*-----------------------------------------------------------------------------
 *  dec_psystem  --  decrements a psystem by subtracting an other psystem
 *                   accepts: psys: pointer to a psystem in standard newton0
 *                                  form
 *                            d_psys: pointer to the psystem with which  psys
 *                                  will be decremented.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_psystem(psys, d_psys)
psystptr  psys;
psystptr  d_psys;
    {
    Tnow(psys) -= Tnow(d_psys);

    dec_pbodies(PBodies(psys), PBodies(d_psys), Nbody(psys));
    }

/*===========================================================================*/

/*****************************************************************************/
/*   _____________________________________________________________________   */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |                              PART V                                 |  */
/*  |                                                                     |  */
/*  |           contains procedures for operations on csystems.           |  */
/*  |                                                                     |  */
/*  |               csystptr  mk_empty_csystem(         )                 |  */
/*  |               csystptr  mk_csystem(                npart )          |  */
/*  |               void      rm_cystem(           sys  )                 |  */
/*  |               void      clear_csystem(      csys  )                 |  */
/*  |               csystptr  cp_csystem(         csys  )                 |  */
/*  |               csystptr  mul_csystem(        csys , factor )         |  */
/*  |               csystptr  add_csystem( csys1, csys2 )                 |  */
/*  |               csystptr  sub_csystem( csys1, csys2 )                 |  */
/*  |               void      scale_csystem(      csys , factor )         |  */
/*  |               void      inc_csystem( csys, dcsys  )                 |  */
/*  |               void      dec_csystem( csys, dcsys  )                 |  */
/*  |                                                                     |  */
/*  |                                                                     |  */
/*  |_____________________________________________________________________|  */
/*                                                                           */
/*****************************************************************************/

/*-----------------------------------------------------------------------------
 *  mk_empty_csystem  --  allocates memory for a many-body csystem in standard 
 *                        newton0 form, but does not allocate memory for the
 *                        cbodies; the number of particles is therefore left
 *                        unspecified.
 *                        returns: new_csys: a pointer to the new csystem.
 *-----------------------------------------------------------------------------
 */
csystptr  mk_empty_csystem()
    {
    csystptr  new_csys;

    new_csys = (csystptr) malloc(sizeof(csyst));
    if (new_csys == NULL)
	error("mk_empty_csystem: not enough memory left for a new csystem\n");

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  mk_csystem  --  allocates memory for a many-body csystem in standard
 *                  newton0 form, and properly initializes the types of
 *                  all cbodies.
 *                  returns: new_csys: a pointer to the new csystem.
 *-----------------------------------------------------------------------------
 */
csystptr  mk_csystem(npart)
int  npart;
    {
    csystptr  new_csys;

    new_csys = (csystptr) malloc(sizeof(csyst));
    if (new_csys == NULL)
	error("mk_csystem: not enough memory left for a new csystem\n");

    Nbody(new_csys) = npart;
    CBodies(new_csys) = mk_cbodies(npart);

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  rm_csystem  --  deallocates memory for a csystem, after deallocating
 *                  memory for the cbodies.
 *                  accept: old_csys: pointer to csystem which is to be removed
 *-----------------------------------------------------------------------------
 */
void  rm_csystem(old_csys)
csystptr  old_csys;
    {
    free(CBodies(old_csys));
    free(old_csys);
    }

/*-----------------------------------------------------------------------------
 *  clear_csystem  --  sets all values to zero in a csystem in standard newton0
 *                     form.
 *                     accepts: csys: pointer to a many-body csystem
 *-----------------------------------------------------------------------------
 */
void  clear_csystem(csys)
csystptr  csys;
    {
    Tnow(csys) = 0.0;
    
    clear_cbodies(CBodies(csys), Nbody(csys));
    }

/*-----------------------------------------------------------------------------
 *  cp_csystem  --  makes a copy of a csystem in standard newton0 form.
 *                 accepts: old_csys: pointer to a many-body csystem;
 *                 returns: new_csys: pointer to the new copy
 *-----------------------------------------------------------------------------
 */
csystptr  cp_csystem(old_csys)
csystptr  old_csys;
    {
    csystptr  new_csys;
    
    new_csys = mk_empty_csystem();
    Nbody(new_csys) = Nbody(old_csys);
    Tnow(new_csys) = Tnow(old_csys);

    CBodies(new_csys) = cp_cbodies(CBodies(old_csys), Nbody(old_csys));    

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  mul_csystem  --  form a new csystem by multiplying an old csystem by
 *                   a scalar.
 *                   accepts: old_csys: pointer to a many-body csystem in
 *                                      standard newton0 form.
 *                        scale_factor: with which every entry in "old_csys"
 *                                      multiplied.
 *                   returns: new_csys: pointer to the new scaled copy.
 *-----------------------------------------------------------------------------
 */
csystptr  mul_csystem(old_csys, scale_factor)
csystptr  old_csys;
real  scale_factor;
    {
    csystptr  new_csys;
    
    new_csys = mk_empty_csystem();
    Nbody(new_csys) = Nbody(old_csys);
    Tnow(new_csys) = scale_factor * Tnow(old_csys);

    CBodies(new_csys) = mul_cbodies(CBodies(old_csys), Nbody(old_csys),
                                                                 scale_factor);

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  add_csystem  --  adds two csystems of cbodies together
 *                   accepts: csys1: pointer to a csystem in standard newton0
 *                                   form;
 *                            csys2: pointer to a second such csystem.
 *                   returns: new_csys: pointer to a new csystem which contains
 *                                     the sum of both old ones.
 *-----------------------------------------------------------------------------
 */
csystptr  add_csystem(csys1, csys2)
csystptr  csys1;
csystptr  csys2;
    {
    csystptr  new_csys;
    
    new_csys = mk_empty_csystem();
    Nbody(new_csys) = Nbody(csys1);
    Tnow(new_csys) = Tnow(csys1) + Tnow(csys2);

    CBodies(new_csys) = add_cbodies(CBodies(csys1), CBodies(csys2),
                                                              Nbody(new_csys));

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  sub_csystem  --  subtracts two csystems of cbodies
 *                   accepts: csys1: pointer to a csystem in standard newton0
 *                                   form;
 *                            csys2: pointer to a second such csystem.
 *                   returns: new_csys: pointer to a new csystem which contains
 *                                     the difference of both old ones.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
csystptr  sub_csystem(csys1, csys2)
csystptr  csys1;
csystptr  csys2;
    {
    csystptr  new_csys;
    
    new_csys = mk_empty_csystem();
    Nbody(new_csys) = Nbody(csys1);
    Tnow(new_csys) = Tnow(csys1) - Tnow(csys2);

    CBodies(new_csys) = sub_cbodies(CBodies(csys1), CBodies(csys2),
                                                              Nbody(new_csys));

    return(new_csys);
    }

/*-----------------------------------------------------------------------------
 *  scale_csystem  --  multiplies an existing csystem by a scalar.
 *                     accepts: old_csys: pointer to a many-body csystem in
 *                                       standard newton0 form.
 *                         scale_factor: with which every entry in "old_csys"
 *                                       multiplied.
 *-----------------------------------------------------------------------------
 */
void  scale_csystem(csys, scale_factor)
csystptr  csys;
real  scale_factor;
    {
    Tnow(csys) *= scale_factor;

    scale_cbodies(CBodies(csys), Nbody(csys), scale_factor);
    }

/*-----------------------------------------------------------------------------
 *  inc_csystem  --  increments a csystem by adding an other csystem
 *                   accepts: csys: pointer to a csystem in standard newton0
 *                                  form
 *                          d_csys: pointer to the csystem with which  csys
 *                                  will be incremented.
 *-----------------------------------------------------------------------------
 */
void  inc_csystem(csys, d_csys)
csystptr  csys;
csystptr  d_csys;
    {
    Tnow(csys) += Tnow(d_csys);

    inc_cbodies(CBodies(csys), CBodies(d_csys), Nbody(csys));
    }

/*-----------------------------------------------------------------------------
 *  dec_csystem  --  decrements a csystem by subtracting an other csystem
 *                   accepts: csys: pointer to a csystem in standard newton0
 *                                  form
 *                            d_csys: pointer to the csystem with which  csys
 *                                  will be decremented.
 *             disclaimer: the management is not responsible for such repulsive
 *                         creations as are formed by subtracting large masses
 *                         from small ones ...
 *-----------------------------------------------------------------------------
 */
void  dec_csystem(csys, d_csys)
csystptr  csys;
csystptr  d_csys;
    {
    Tnow(csys) -= Tnow(d_csys);

    dec_cbodies(CBodies(csys), CBodies(d_csys), Nbody(csys));
    }

/*===========================================================================*/

/* endof: systemalgebra.c */
