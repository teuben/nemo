/* save.c - init_save_state, restore_ctrl, restore_diag, restore_spec, 
            restore_state, restore_system, save_ctrl, save_diag, save_spec, 
            save_state, save_system, write_restart_command */

/*
 *  save.c: saving/restoring module for newton0.c : equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *      Feb  2001  -  PJT  fixed specsptr->diagptr in definition of restore_diag
 */
   
#include  "newton0.h"
#include  <filestruct.h>

static void write_restart_command(string restartfile);
static void save_system(systptr sys, stream savestr);
static void save_spec(specptr specs, stream savestr);
static void save_ctrl(ctrlptr ctr, stream savestr);
static void save_diag(diagptr diags, stream savestr);
static void save_regsystem(systptr regsys, stream savestr);  
static systptr restore_system(stream restorestr);
static void restore_spec(specptr specs, stream restorestr);
static void restore_ctrl(ctrlptr ctr, stream restorestr);
static void restore_diag(diagptr diags, stream restorestr);  
static systptr restore_regsystem(stream restorestr);

/*-----------------------------------------------------------------------------
 *  init_save_state  --  checks whether a safefile is specified; if so, starts
 *                       a first system saving, and checks whether a
 *                       restartfile has been specified; if so, outputs the
 *                       appropriate restart command in the restartfile, to
 *                       enable automatic restart after a system crash.
 *-----------------------------------------------------------------------------
 */
void  init_save_state(the_state)
stateptr  the_state;
    {
    if (*Savefile(Ctrls(the_state)) == 0)                  /* no saving?     */
	{                                                  /* then avoid     */
        Tsave(Ctrls(the_state)) = Tend(Ctrls(the_state));  /* future         */
        DTsave(Ctrls(the_state)) = 0.0;                    /* attempts and   */
	return;                                            /* return quietly.*/
	}

    if (*Restartfile(Ctrls(the_state)) != 0)               /* restart desired?   */
        write_restart_command(Restartfile(Ctrls(the_state)));
    }

/*-----------------------------------------------------------------------------
 *  write_restart_command  --  write in the restartfile the same command with
 *                             which the present run was started, but with the
 *                             modification of using the restartfilename in
 *                             the restore option.
 *-----------------------------------------------------------------------------
 */
local void  write_restart_command(restartfile)
string  restartfile;
    {
    stream restartstr;

    restartstr = stropen(restartfile, "w");

    printf("\n write_restart_command: not yet implemented (ask Josh)\n");
    }

/*-----------------------------------------------------------------------------
 *  save_state  --  saves the full state of the system, which contains all the
 *                  information needed to restart and resume the integration.
 *                  save_state()  is invoked, when desired, at times in between
 *                  calls to  propagate() , as orchestrated by  evolve()  in
 *                  the file  orbit.c .
 *                  the order of output is irrelevant, as long as it is
 *                  duplicated exactly in  restore_state()  below.
 *-----------------------------------------------------------------------------
 */
void  save_state(the_state)
stateptr  the_state;
    {
    stream savestr;

    savestr = stropen(Savefile(Ctrls(the_state)), "a");
    fseek(savestr, 0L, 0);			/* rewind to origin */

    save_system(System(the_state), savestr);
    save_spec(Specs(the_state), savestr);
    save_ctrl(Ctrls(the_state), savestr);
    save_diag(Diags(the_state), savestr);
#ifdef REGULARIZATION
    save_regsystem(Regsystem(the_state), savestr);
#endif

    strclose(savestr);
/*
 * the following announcement is channeled through a procedure in the file
 * out.c  for modularity's sake, leaving final output to  out.c .
 */
    announce("\n\tSystem saving done at time t_now = %lf\n",
                                                      Tnow(System(the_state)));
/*
 * schedule the next time for saving the system state:
 */
    Tsave(Ctrls(the_state)) = Tnow(System(the_state))+DTsave(Ctrls(the_state));
    }

/*-----------------------------------------------------------------------------
 *  restore_state  --  restores the full state of the system, from a dump made
 *                     by save_state() in a previous run.
 *-----------------------------------------------------------------------------
 */
stateptr  restore_state(old_file)
string  old_file;
    {
    stateptr  old_state;
    stream restorestr;
    string local_announcement = "\tRestoring an old system state";

    if (*old_file == 0)
	error("restore specified without old_file name");
    restorestr = stropen(old_file, "r");

    old_state = mk_state();                           /* in  statealgebra.c  */

    System(old_state) = restore_system(restorestr);
    restore_spec(Specs(old_state), restorestr);
    restore_ctrl(Ctrls(old_state), restorestr);
    restore_diag(Diags(old_state), restorestr);
#ifdef REGULARIZATION
    Regsystem(old_state) = restore_regsystem(restorestr);
#endif

    strclose(restorestr);

    Announcement(Ctrls(old_state)) = local_announcement;

    return(old_state);
    }

/*-----------------------------------------------------------------------------
 *  save_system  --  saves the system part of a state
 *-----------------------------------------------------------------------------
 */
local void  save_system(sys, savestr)
systptr  sys;
stream  savestr;
    {
    put_data(savestr, "nbody", IntType, &Nbody(sys), 0);
    put_data(savestr, "t_now", RealType, &Tnow(sys), 0);
    put_data(savestr, "bodies", AnyType, Bodies(sys), Nbody(sys),
                                                              sizeof(body), 0);
#ifdef TREE
    put_data(savestr, "ncell", IntType, &Ncell(sys), 0);
    error("TREE saving not yet implemented, sorry\n");
#endif
    }

/*-----------------------------------------------------------------------------
 *  save_spec  --  saves the specifications part of a state
 *-----------------------------------------------------------------------------
 */
local void  save_spec(specs, savestr)
specptr  specs;
stream  savestr;
    {
    int  k;
    char  label[32];

    put_data(savestr, "specifications", AnyType, specs, sizeof(spec), 0);

    for (k = 0; k < Number_of_methods; k++)
        {
	sprintf(label, "method%d", k);
	put_string(savestr, label, Methods(specs)[k]);
	}
    }

/*-----------------------------------------------------------------------------
 *  save_ctrl  --  saves the control part of a state
 *-----------------------------------------------------------------------------
 */
local void  save_ctrl(ctr, savestr)
ctrlptr  ctr;
stream  savestr;
    {
    int  k;
    char  label[32];

    put_data(savestr, "controls", AnyType, ctr, sizeof(ctrl), 0);

    for (k = 0; k < Number_of_filenames; k++)
        {
	sprintf(label, "filename%d", k);
	put_string(savestr, label, Filenames(ctr)[k]);
	}

    for (k = 0; k < Number_of_messages; k++)
        {
	sprintf(label, "message%d", k);
	put_string(savestr, label, Messages(ctr)[k]);
	}
    }

/*-----------------------------------------------------------------------------
 *  save_diag  --  saves the diagnostics part of a state
 *-----------------------------------------------------------------------------
 */
local void  save_diag(diags, savestr)
diagptr  diags;
stream  savestr;
    {
    put_data(savestr, "diagnostics", AnyType, diags, sizeof(diag), 0);
    }

#ifdef REGULARIZATION

/*-----------------------------------------------------------------------------
 *  save_regsystem  --  saves the regularized system part of a state
 *-----------------------------------------------------------------------------
 */
local void  save_regsystem(regsys, savestr)
systptr  regsys;
stream  savestr;
    {
    put_data(savestr, "nbody", IntType, &Nbody(regsys), 0);
    put_data(savestr, "t_now", RealType, &Tnow(regsys), 0);
    put_data(savestr, "bodies", AnyType, Bodies(regsys), Nbody(regsys),
                                                              sizeof(body), 0);
    put_data(savestr, "massmatrix", RealType, Massmatrix(regsys),
             Nbody(regsys) * Nbody(regsys), 0);
    put_data(savestr, "regenergy", RealType, &Regenergy(regsys), 0);
    put_data(savestr, "reglagrangian", RealType, &Reglagrangian(regsys), 0);
    put_data(savestr, "reghamiltonian", RealType, &Reghamiltonian(regsys), 0);
    }

#endif

/*-----------------------------------------------------------------------------
 *  restore_system  --  restores the system part of a state
 *-----------------------------------------------------------------------------
 */
local systptr  restore_system(restorestr)
stream  restorestr;
    {
    int  nbody;
    systptr  old_sys;

    get_data(restorestr, "nbody", IntType, &nbody, 0);

    old_sys = mk_system(nbody);

    get_data(restorestr, "t_now", RealType, &Tnow(old_sys), 0);
    get_data(restorestr, "bodies", AnyType, Bodies(old_sys), nbody,
             sizeof(body), 0);

    return(old_sys);
    }

/*-----------------------------------------------------------------------------
 *  restore_spec  --  restores the specifications part of a state
 *-----------------------------------------------------------------------------
 */
local void  restore_spec(specs, restorestr)
specptr  specs;
stream  restorestr;
    {
    int  k;
    char  label[32];

    get_data(restorestr, "specifications", AnyType, specs, sizeof(spec), 0);

    for (k = 0; k < Number_of_methods; k++)
        {
	sprintf(label, "method%d", k);
	Methods(specs)[k] = get_string(restorestr, label);
	}
    }

/*-----------------------------------------------------------------------------
 *  restore_ctrl  --  restores the controls part of a state
 *-----------------------------------------------------------------------------
 */
local void  restore_ctrl(ctr, restorestr)
ctrlptr  ctr;
stream  restorestr;
    {
    int  k;
    char  label[32];

    get_data(restorestr, "controls", AnyType, ctr, sizeof(ctrl), 0);

    for (k = 0; k < Number_of_filenames; k++)
        {
	sprintf(label, "filename%d", k);
	Filenames(ctr)[k] = get_string(restorestr, label);
	}

    for (k = 0; k < Number_of_messages; k++)
        {
	sprintf(label, "message%d", k);
	Messages(ctr)[k] = get_string(restorestr, label);
	}
    }

/*-----------------------------------------------------------------------------
 *  restore_diag  --  restores the specifications part of a state
 *-----------------------------------------------------------------------------
 */
local void  restore_diag(diags, restorestr)
diagptr  diags;
stream  restorestr;
    {
    get_data(restorestr, "diagnostics", AnyType, diags, sizeof(diag), 0);
    }

#ifdef REGULARIZATION

/*-----------------------------------------------------------------------------
 *  restore_regsystem  --  restores the regularized system part of a state
 *-----------------------------------------------------------------------------
 */
local systptr  restore_regsystem(restorestr)
stream  restorestr;
    {
    int  nbody;
    systptr  old_regsys;

    get_data(restorestr, "nbody", IntType, &nbody, 0);

    old_regsys = mk_system(nbody);

    get_data(restorestr, "t_now", RealType, &Tnow(old_regsys), 0);
    get_data(restorestr, "bodies", AnyType, Bodies(old_regsys), nbody,
             sizeof(body), 0);

    Massmatrix(old_regsys) = mk_massmatrix(nbody);

    get_data(restorestr, "massmatrix", RealType, Massmatrix(old_regsys),
             nbody * nbody, 0);
    get_data(restorestr, "regenergy", RealType, &Regenergy(old_regsys), 0);
    get_data(restorestr, "reglagrangian", RealType,
             &Reglagrangian(old_regsys), 0);
    get_data(restorestr, "reghamiltonian", RealType,
             &Reghamiltonian(old_regsys), 0);

    return(old_regsys);
    }

#endif

/* endof: save.c */
