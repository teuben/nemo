/* ctrl.h - Announcement, CPUlastcall, CPUmax, DEmax, DNmajor, DTimes, DTmajor,
            DTmax, DTminor, DTsave, DTstep, Filenames, Forwards, Headline, 
            Ilimits, Infile, Lastoutfile, Messages, Min_pairdist, Nmajor,
            Nmaxstep, Ntimingiter, Number_of_dtimes, Number_of_filenames,
            Number_of_ilimits, Number_of_messages, Number_of_ntimes, 
            Number_of_rlimits, Number_of_times, Outfile, Restorefile, 
            Resumefile, Rlimits, Savefile, Times, Tend, Tmajor, Tminor,
            Tsave, ctrl, ctrlptr */

/*
 *  ctrl.h: for "ctrl", a substructure of a "state".
 *
 *      June 1987  -  Piet Hut  @ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 */
   
/*-----------------------------------------------------------------------------
 *  ctrl, ctrlptr  --  contains parameters which control the overall
 *                     activities: duration of orbit integration,
 *                     scheduling and directing input/output, etc.
 *                     note: recipe for adding an extra element:
 *                           1) increase the appropriate array length by
 *                              adding 1 to Number_of_... ;
 *                           2) introduce a macro below to provide a handle
 *                              to access the new element;
 *                           3) provide a short description at the end of this
 *                              file to document the meaning of the new element
 *                           in addition, if the new element should be read
 *                           from the command line (which is generally the
 *                           case) then:
 *                           4) add an entry to the  defv[]  string at the top
 *                              of the file  newton0.c;
 *                           5) add a line to the procedure  set_controls()
 *                              in the file  newton0.c , using one of the 
 *                              getparam tools.
 *                           no change has to be made to the procedures in
 *                           the file  save.c , unless a new type of array is
 *                           added to  ctrl .
 *-----------------------------------------------------------------------------
 */
#define    Number_of_filenames    7
#define    Number_of_messages	  2
#define    Number_of_times	  4
#define    Number_of_dtimes	  4
#define    Number_of_ntimes	  1
#define    Number_of_dntimes	  1
#define    Number_of_rlimits	  5
#ifndef REGULARIZATION
#define    Number_of_ilimits	  1
#endif
#ifdef REGULARIZATION
#define    Number_of_ilimits	  2
#endif

typedef struct
    {
    string  filenames[ Number_of_filenames ];
    string  messages[  Number_of_messages  ];
    real  times[       Number_of_times     ];
    real  dtimes[      Number_of_dtimes    ];
    int   ntimes[      Number_of_ntimes    ];
    int   dntimes[     Number_of_dntimes   ];
    real  rlimits[     Number_of_rlimits   ];
    int   ilimits[     Number_of_ilimits   ];
    bool  forwards;
    } ctrl, *ctrlptr;

/*-----------------------------------------------------------------------------
 *  macros to extract arrays of components from a ctrl structure:
 *  NOTE:  ptr should be a pointer of type "ctrlptr":
 *-----------------------------------------------------------------------------
 */
#define  Filenames(ptr)      ((ptr)->filenames)          /* type: stringptr  */
#define  Messages(ptr)       ((ptr)->messages)           /* type: stringptr  */
#define  Times(ptr)          ((ptr)->times)              /* type: realptr    */
#define  DTimes(ptr)         ((ptr)->dtimes)             /* type: realptr    */
#define  NTimes(ptr)         ((ptr)->ntimes)             /* type: intptr     */
#define  DNTimes(ptr)        ((ptr)->dntimes)            /* type: intptr     */
#define  Rlimits(ptr)        ((ptr)->rlimits)            /* type: realptr    */
#define  Ilimits(ptr)        ((ptr)->ilimits)            /* type: intptr     */

/*-----------------------------------------------------------------------------
 *  macros to extract individual components from a ctrl structure (see next
 *  page for a description of their meaning):
 *  NOTE:  ptr should be a pointer of type "ctrlptr"
 *-----------------------------------------------------------------------------
 */
#define  Forwards(ptr)           ((ptr)->forwards)         /* type: bool     */

#define  Infile(ptr)             ((ptr)->filenames[0])     /* type: string   */
#define  Outfile(ptr)            ((ptr)->filenames[1])     /* type: string   */
#define  Lastoutfile(ptr)        ((ptr)->filenames[2])     /* type: string   */
#define  Savefile(ptr)           ((ptr)->filenames[3])     /* type: string   */
#define  Resumefile(ptr)         ((ptr)->filenames[4])     /* type: string   */
#define  Restorefile(ptr)        ((ptr)->filenames[5])     /* type: string   */
#define  Restartfile(ptr)        ((ptr)->filenames[6])     /* type: string   */

#define  Headline(ptr)           ((ptr)->messages[0])      /* type: string   */
#define  Announcement(ptr)       ((ptr)->messages[1])      /* type: string   */

#define  Tend(ptr)               ((ptr)->times[0])         /* type: real     */
#define  Tmajor(ptr)             ((ptr)->times[1])         /* type: real     */
#define  Tminor(ptr)             ((ptr)->times[2])         /* type: real     */
#define  Tsave(ptr)              ((ptr)->times[3])         /* type: real     */

#define  DTmajor(ptr)            ((ptr)->dtimes[0])        /* type: real     */
#define  DTminor(ptr)            ((ptr)->dtimes[1])        /* type: real     */
#define  DTsave(ptr)             ((ptr)->dtimes[2])        /* type: real     */
#define  DTstep(ptr)             ((ptr)->dtimes[3])        /* type: real     */

#define  Nmajor(ptr)             ((ptr)->ntimes[0])        /* type: int      */

#define  DNmajor(ptr)            ((ptr)->dntimes[0])       /* type: int      */

#define  DTmax(ptr)              ((ptr)->rlimits[0])       /* type: real     */
#define  CPUmax(ptr)             ((ptr)->rlimits[1])       /* type: real     */
#define  CPUlastcall(ptr)        ((ptr)->rlimits[2])       /* type: real     */
#define  Min_pairdist(ptr)       ((ptr)->rlimits[3])       /* type: real     */
#define  DEmax(ptr)              ((ptr)->rlimits[4])       /* type: real     */

#define  Nmaxstep(ptr)           ((ptr)->ilimits[0])       /* type: int      */
#ifdef REGULARIZATION
#  define  Ntimingiter(ptr)      ((ptr)->ilimits[1])       /* type: int      */
#endif

/*-----------------------------------------------------------------------------
 *  short descriptions of the various components of a "ctrl" structure:
 *
 *
 *      Forwards()           indicates the time direction of the integration;
 *                           initially set to  TRUE  if   Tnow()  <  Tend() ;
 *                           initially set to  FALSE if   Tnow()  >  Tend() .
 *
 *      Infile()             binary input file.
 *      Outfile()            binary output file.
 *      Lastoutfile()        file to which the last binary output is delivered.
 *      Savefile()           file in which to save data from the present run.  
 *      Resumefile()         file containing data saved from a previous run,   
 *                           which will be resumed with new control parameters.
 *      Restorefile()        file containing data saved from a previous run,   
 *                           which will be continued without any change.       
 *      Restartfile()        file in which to write the restart command for the
 *                           current run, in case the current run is aborted.
 *
 *      Headline()           output message describing the type of calculation,
 *                           initialized in the file  newton0.c
 *      Announcement()       message describing initial set-up activities,
 *                           initialized in one of the files
 *                           create.c , in.c , or  save.c . 
 *
 *      Tend()               time at which the integration should halt
 *      Tmajor()             time of next major output
 *      Tminor()             time of next minor output
 *      Tsave()              time at which system state will be saved
 *
 *      DTmajor()            time interval between major outputs
 *      DTminor()            time interval between minor outputs
 *      DTsave()             time interval between saving system state
 *      DTstep()             common time step size for all particles
 *
 *      Nmajor()             number of steps after which integration should
 *                           halt.
 *
 *      DNmajor()            number of integration steps between major outputs.
 *
 *      DTmax()              maximum integration time step length.
 *      CPUmax()             maximum allowed amount of CPU time (in minutes).
 *      CPUlastcall()        amount of CPU time after which integration will
 *                           be halted at next major-output time (in minutes).
 *      Min_pairdist()       minimum distance allowed between any pair of
 *                           particles in the system; passing this limit causes
 *                           integration to be halted.
 *      DEmax()              maximum drift allowed in the total energy, before
 *                           integration is interrupted.
 *
 *      Nmaxstep()           maximum number of integration steps allowed
 *      Ntimingiter()        number of iterations for arriving at the desired
 *                           output times using regularization, by choosing
 *                           increasingly smaller steps in pseudo-time.
 *-----------------------------------------------------------------------------
 */

/* endof: ctrl.h */
