/* out.c - announce, ascii_out, init_ascii_out, maj_out, min_out,
           talkin_talkin,  write_initial_output, write_final_output */

/*
 *  out.c: output module for newton0.c : equal time steps
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
#include  "newton0.h"

static real proj_free_time(systptr sys, real r0_soft);

extern real  cputime();    /* a procedure using a UNIX times(3C) call to return */
                        /* the CPU time used at the time of invoking         */
                        /* cputime() while executing instructions in the     */
                        /* in the user space of the current process.         */


/*-----------------------------------------------------------------------------
 *  write_initial_output  --  start with major/minor output and save the state
 *-----------------------------------------------------------------------------
 */
void  write_initial_output(new_state)
stateptr  new_state;
    {
    ctrlptr  ctr;
    
    ctr = Ctrls(new_state);
/*
 * open the output file, if specified:
 */
    if (*Outfile(ctr) != 0)                            /* output file given? */
        open_binary_output(Outfile(ctr), Headline(ctr));  /* in  binaryout.c */
/*
 * announce the present run, and print the most interesting parameter values:
 */
    init_ascii_out(new_state);
/*
 * decide which, if any, type of output to perform (major, minor, or none):
 */
    if (too_short(DTmajor(ctr)))                /* no major output desired?  */
        {
	Tmajor(ctr) = Tend(ctr);                /* then so be it.            */
        if (too_short(DTminor(ctr)))            /* no minor output desired?  */
	    Tminor(ctr) = Tend(ctr);            /* then so be it,            */
        else                                    /* otherwise schedule it,    */
            min_out(new_state, TRUE);           /* replacing a major output. */
        }
    else                                        /* schedule major output.    */
        {
        maj_out(new_state, TRUE);
        if (too_short(DTminor(ctr)))            /* no minor output desired?  */
	    Tminor(ctr) = Tend(ctr);            /* then we are done, since   */
	}                                       /* maj_out() schedules minor */
/*
 * decide whether to start saving the system state:
 */	
    if (too_short(DTsave(ctr)))                 /* no system saving desired? */
	Tsave(ctr) = Tend(ctr);                 /* then so be it,            */
    else                                        /* otherwise schedule it.    */
        init_save_state(new_state);
    }

/*-----------------------------------------------------------------------------
 *  write_final_output  --  perform the last output and close all open files
 *-----------------------------------------------------------------------------
 */
void  write_final_output(old_state)
stateptr  old_state;
    {
    ctrlptr  ctr;

    ctr = Ctrls(old_state);
/*
 * first do a last major or minor output, if desired:
 */
    if (DIAGnsteps(Diags(old_state)) > 0)  /* if not, no need for 2nd output */
        {
        if (! too_short(DTmajor(ctr)))          /* major output desired?     */
            maj_out(old_state, FALSE);
        else if (! too_short(DTminor(ctr)))     /* minor output desired?     */
            min_out(old_state, FALSE);
        }

/*
 * anything more required besides the normal minor or major output?
 */
    final_diagnostics(old_state);          /* resides in file  diagnose.c   */

    /* ...... */                  /* write final diagnostics */
/*
 * when done with the final output, close the output file:
 */
    if (*Outfile(ctr) != 0)                           /* output file given?  */
        close_binary_output();                  /*   in  binaryout.c   */

    if (*Lastoutfile(ctr) != 0)                   /* last-output file given? */
        {
        open_binary_output(Lastoutfile(ctr), Headline(ctr));  /* binaryout.c */

        diagnostics(old_state, FALSE);   /* not clear how this connects with */
                                         /*  final_diagnostics()  above.     */
#ifdef REGULARIZATION
        reg_to_nonreg(Regsystem(old_state), System(old_state));
#endif
        put_full_snapshot(old_state);

        close_binary_output();                  /*   in  binaryout.c   */
	}
/*
 * the last save operation should be performed at the very last,
 * since we should be able to continue another integration,
 * restarting from the final saved system state
 * without any loss of either binary or standard output:
 */
    if (! too_short(DTsave(ctr)))               /* no system saving desired? */
        save_state(old_state);                  /* save the state.           */
    }

/*-----------------------------------------------------------------------------
 *  maj_out  --  takes care of major output operation
 *-----------------------------------------------------------------------------
 */
void  maj_out(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    systptr  sys;
    ctrlptr  ctr;
    diagptr  diags;
    
    sys = System(the_state);
    ctr = Ctrls(the_state);
    diags = Diags(the_state);

    diagnostics(the_state, init_flag);

    ascii_out(the_state);
/*
 * binary output, if desired:
 */
    if (*Outfile(ctr) != 0)                           /* output file given?  */
        {
#ifdef REGULARIZATION
        reg_to_nonreg(Regsystem(the_state), sys);
#endif
        put_full_snapshot(the_state);
        }
/*
 * while a major output makes a simultaneous minor output unnecessary, 
 * it has to reschedule the next minor output, as well as the major output:
 */
    Tmajor(ctr) = Tnow(sys) + DTmajor(ctr);
    if (! too_short(DTminor(ctr)))                 /* minor output desired?  */
        Tminor(ctr) = Tnow(sys) + DTminor(ctr);    /* then schedule it.      */
    Nmajor(ctr) = DIAGnsteps(diags) + DNmajor(ctr);
    }


/*-----------------------------------------------------------------------------
 *  min_out  --  takes care of minor output operation
 *-----------------------------------------------------------------------------
 */
void  min_out(the_state, init_flag)
stateptr  the_state;
bool  init_flag;
    {
    systptr  sys;
    ctrlptr  ctr;
    diagptr  diags;
    
    sys = System(the_state);
    ctr = Ctrls(the_state);
    diags = Diags(the_state);

    diagnostics(the_state, init_flag);

    ascii_out(the_state);
/*
 * binary output, if desired:
 */
    if (*Outfile(ctr) != 0)                          /* output file given?  */
        put_diagnostics_snapshot(the_state);
/*
 * schedule the next minor output:
 *   note the following minor point:  t_min_out += dt_min_out would introduce
 *   the possibility that successive ouput times get more and more out of sync,
 *   if t_maj_out and t_save slightly differ from t_min_out. Therefore here
 *   and in other output procedures the use of t_now is preferred.
 */
    Tminor(ctr) = Tnow(sys) + DTminor(ctr);
    }

/*-----------------------------------------------------------------------------
 *  init_ascii_out  --  prints initial information on the standard output
 *-----------------------------------------------------------------------------
 */
void  init_ascii_out(new_state)
stateptr  new_state;
    {
    systptr  sys;
    specptr  specs;
    ctrlptr  ctr;
    
    sys = System(new_state);
    specs = Specs(new_state);
    ctr = Ctrls(new_state);

/*
 * announce the present run, and print the most interesting parameter values:
 */
    printf("\n\t%s,\n", Headline(ctr));  /* print headline      */
#ifndef EXTRAPOLATION
    printf("\t             integration scheme: %s\n",Integrationscheme(specs));
    printf("\t            time step criterion: %s\n", Timestepmethod(specs));
#else                                         
    printf("\t             integration scheme: predictor-semicorrector\n");
    printf("\t            time step criterion: Aarseth's force criterion\n");
#endif
    if(Softparam(specs) > 0.0)
        printf("\t               softening method: %s\n", Softfocus(specs));
    if (*Outfile(ctr) == 0 && *Lastoutfile(ctr) == 0)  /* output file? */
	printf("\n\twarning: no binary output file name provided\n");

    printf("\n\tnbody = %d    eta_acc = %.2g    r0_soft = %.2g",
                               Nbody(sys), Stepparam(specs), Softparam(specs));
#ifdef TREE
    printf("    tol = %.2g", sqrt(Tolsqparam(specs)));
#endif
#ifdef REGULARIZATION
    printf("    niter = %d", Ntimingiter(ctr));
#endif
    printf("\n");

/*
 * announce the initial set-up activities which have been performed:
 */
    printf("\n%s\n", Announcement(ctr));
    }

/*-----------------------------------------------------------------------------
 *  ascii_out  --  prints information on the standard output
 *-----------------------------------------------------------------------------
 */
void  ascii_out(the_state)
stateptr  the_state;
    {
    systptr  sys;
    specptr  specs;
    diagptr  diags;
#ifdef TREE
    real  nbavg;  /* average number of body-body force calculations per body */
    real  ncavg;  /* average number of body-cell force calculations per body */
    real  fcell;      /* ratio of numbers of cells and bodies                */
#endif
    
    sys = System(the_state);
    specs = Specs(the_state);
    diags = Diags(the_state);

#ifdef TREE
    nbavg = ((real) N2bcalc(specs) / (real) Nfcalc(specs));
    ncavg = ((real) Nbccalc(specs) / (real) Nfcalc(specs));
    fcell = ((real) Ncell(sys) / (real) Nbody(sys));
/*
 * reset to start measuring till next major or minor output:
 */
    Nfcalc(specs) = N2bcalc(specs) = Nbccalc(specs) = 0;
#endif

/*
 * diagnostics on the standard output:
 */
#ifndef TREE
    printf("\n  %8s%7s%7s%12s%7s%8s%10s%9s\n", "t_now",
	   "T+U", "T/U", "(E-E0)/E0", "dE/E0", "nsteps", "pair_min","cputime");
    printf("  %8.3f%8.4f%8.4f%9.2g%9.2g%7d%9.2g%10.2g\n",
           DIAGtime(diags)[CURRENT], DIAGetot(diags)[CURRENT],
           DIAGekin(diags)[CURRENT] / DIAGepot(diags)[CURRENT], 
           RELDRIFT(DIAGetot(diags)), RELINCREMENT(DIAGetot(diags)), 
           DIAGnsteps(diags), SPECmin_pair(specs), cputime());
#endif
#ifdef TREE
    printf("\n  %8s%7s%7s%12s%7s%8s%7s%7s%6s%8s\n", "t_now",
	      "T+U", "T/U", "(E-E0)/E0", "dE/E0", "nsteps",
                                           "<bb>", "<bc>", "fcell", "cputime");
    printf("  %8.3f%8.4f%8.4f%9.2g%9.2g%7d%7.3g%7.3g%6.2g%8.2g\n",
           DIAGtime(diags)[CURRENT], DIAGetot(diags)[CURRENT],
           DIAGekin(diags)[CURRENT] / DIAGepot(diags)[CURRENT], 
           RELDRIFT(DIAGetot(diags)), RELINCREMENT(DIAGetot(diags)), 
           DIAGnsteps(diags), nbavg, ncavg, fcell, cputime());
#endif
    }


/*-----------------------------------------------------------------------------
 *  talkin_talkin  --  provides lots of information, in case you're interested
 *-----------------------------------------------------------------------------
 */
void  talkin_talkin()
    {
    announce("\ntalkin_talkin: not yet implemented\n", 0.0);
    }

/*-----------------------------------------------------------------------------
 *  announce  --  prints a string and a real number (generally the time)
 *                describing a binary output action on the standard
 *                output file.
 *                This procedure collects such announcements from procedures
 *                in several different files, and acts as a barrier between
 *                levels of abstraction: procedures in other files are not
 *                allowed to print directly on the standard output, but have to
 *                use the present procedure as an intermediary.
 *-----------------------------------------------------------------------------
 */
void  announce(what_happened, when_did_it_happen)
string  what_happened;
real  when_did_it_happen;
    {
    printf(what_happened, when_did_it_happen);
    }    

/* endof: out.c */
