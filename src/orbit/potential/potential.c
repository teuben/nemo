/* potential.c - proc get_potential, local proc load_potential */
/*               double get_pattern                            */

/*------------------------------------------------------------------------------
 *  POTENTIAL:  procedures for finding fixed potential, initially written 
 *		for newton0
 *              It has some nasty assumptions if fortran code is loaded,
 *              in particular names of symbols, and f2c interfaces
 *
 *      July 1987  -  Peter Teuben  @ Inst. f. Adv. Study, 
 *						  Princeton, NJ 08540, USA
 *	April 1988 -  V2.0 added dataname parameter to get_potential
 *	March 1989 -  V2.1 added dataname length in case fortran was called 
 *			   (only BSD interface fortran-C)
 *	Febr. 1990 -  V3.0 added time as extra parameter in potential() call 
 *			   although no extra code needed here...
 *			   drange() replaced by nemoinpd()		 - PJT
 *	14-mar-90     V3.1 made GCC happy - increased MAXPAR for Stefano - PJT
 *      25-oct-90     V3.2 allow multiple get_potential() calls		   PJT
 *	 6-nov-90     V3.2a force link of spline(3NEMO) routines	   PJT
 *	20-nov-90         b NEMOPATH->NEMO and few other cosmetics         PJT
 *      14-oct-91     V3.3  add image I/O to make potfft a potential(5)    PJT
 *	 7-mar-92     V3.3a happy gcc2.0                                   PJT
 *	24-may-92     V3.4 changed names of local variables to accomodate
 *			   the new 'potential' data type.		   PJT
 *      11-oct-93     V5.0 get_pattern                                     PJT
 *	20-jan-94     V5.1 using mapsys(); needed for solaris		   pjt
 *	22-feb-94     V5.1a ansi
 *	23-mar-95         b prototypes	- fixed 0-strlen bug in f77 I/O    pjt
 *	17-apr-95           use -DNO_IMAGE				   pjt
 *      13-dec-95     V5.2  Now using ldso to access shared objects        pjt 
 *                          fixed small CFLAGS bug
 *	 1-apr-01     V5.3  converted for NEMO V3 with .so files           pjt
 *      13-sep-01     V5.4  support for potential_double and potential_float   pjt
 *      18-sep-01           auto-detecting which type is present
 *      17-may-02     V5.4a fix ambiguity about float/double if generic present WD
 *      25-apr-04         b fool optimizing compilers                           PJT
 *      20-may-04         c add sqr() for dummy linker
 *      14-jul-05         d made dummy functions global, for new (FC4) linker 
 *      18-sep-08         e make 'r' == SINGLEPREC? 'f' : 'd'              WD
 *       
 *------------------------------------------------------------------------------
 */

#include  <stdinc.h>
#include  <getparam.h>
#include  <loadobj.h>
#include  <filefn.h>
#include  <potential.h>

#define MAXPAR 64

local double local_par[MAXPAR]; /* NOTE: first par reserved for pattern speed*/
local int  local_npar=0;        /* actual used number of par's               */
local real local_omega=0;	/* pattern speed                             */
local proc l_potential=NULL;    /* actual storage of pointer to exter worker */
local proc l_inipotential=NULL; /* actual storage of pointer to exter inits  */
local bool Qfortran = FALSE;    /* was a fortran routine used ? -- a hack -- */
local bool first = TRUE;        /* see if first time called for mysymbols()  */

void potential_dummy_for_c(void);

/* forward declarations */

local proc load_potential(string, string, string, char); /* load by name    */

/*-----------------------------------------------------------------------------
 *  get_potential --  returns the pointer ptr to the function which carries out
 *                the calculation of potential and accelerations.
 *-----------------------------------------------------------------------------
 */
proc  get_potential(string potname, string potpars, string potfile)
{
    if (potname == NULL || *potname == 0)	/* if no name provided */
        return NULL;				/* return no potential */
    l_potential = load_potential(potname, potpars, potfile,'r');
    return l_potential;
}

potproc_double get_potential_double(string potname, string potpars, string potfile)
{
    if (potname == NULL || *potname == 0)	/* if no name provided */
        return NULL;				/* return no potential */
    l_potential = load_potential(potname, potpars, potfile,'d');
    return (potproc_double) l_potential;
}

potproc_float  get_potential_float(string potname, string potpars, string potfile)
{
    if (potname == NULL || *potname == 0)	/* if no name provided */
        return NULL;				/* return no potential */
    l_potential = load_potential(potname, potpars, potfile,'f');
    return (potproc_float) l_potential;
}

/*-----------------------------------------------------------------------------
 *  get_inipotential --  returns the pointer ptr to the last inipotential
 *          function which initializes the potential
 *-----------------------------------------------------------------------------
 */
proc get_inipotential()
{
    if (first) error("get_inipotential: get_potential not called yet");
    return l_inipotential;
}

/*-----------------------------------------------------------------------------
 *  get_pattern --  return the last set pattern speed
 *		    note that a 0.0 parameter does *not* reset it !!!
 *-----------------------------------------------------------------------------
 */

real get_pattern()
{
    if (first) error("get_pattern: get_potential not called yet");
    return local_omega;
}

/*-----------------------------------------------------------------------------
 *  load_potential -- load the potential from an object file
 *       This routine depends heavily on the object-loader (loadobj(3NEMO))
 *	BUG: if dataname is NULL, the routines dies when trying f_c interface 
 *		since it can't take strlen(NULL)
 *-----------------------------------------------------------------------------
 */
local proc load_potential(string fname, string parameters, string dataname, char type)
{
    char  name[256], cmd[256], path[256], pname[32];
    char  *fullname, *nemopath, *potpath;
    proc  pot, ini_pot;
    int never=0;

    if (parameters!=NULL && *parameters!=0) {              /* get parameters */
	local_npar = nemoinpd(parameters,local_par,MAXPAR);
        if (local_npar>MAXPAR)
            error ("get_potential: potential has too many parameters (%d)",
                            local_npar);
        if (local_npar<0)
            warning("get_potential: parsing error in: %s",parameters);
    } else
        local_npar=0;
    if (local_npar > 0 && local_par[0] != 0.0) { /* aid multiple potentials */
        local_omega = local_par[0];
        dprintf(1,"get_potential: setting local_omega = %g\n",local_omega);
    }

    if (first) {
        mysymbols(getparam("argv0"));      /* get symbols for this program */
        first = FALSE;			   /* and tell it we've initialized */
    }
    potpath = getenv("POTPATH");	     /* is there customized path ? */
    if (potpath==NULL) {			/* use default path */
       potpath = path;
       strcpy (path,".");
       nemopath = getenv("NEMO");
       if (nemopath==NULL) {
	 strcat(path,":");
	 strcat (path,nemopath);
	 strcat (path,"/obj/potential");		/* ".:$NEMO/obj/potential" */
       }
    }
    strcpy (name,fname);
    strcat (name,".so");
    fullname = pathfind (potpath, name);
    if (fullname!=NULL) {			/* .o found !! */
        dprintf (2,"Attempt to load potential from %s\n",name);
    	loadobj(fullname);
    } else {					/* no .o found */
	strcpy (path,".");	/* just look in current directory */
	strcpy (name,fname);
	strcat (name,".c");
	if (pathfind(".",name)==NULL)
	    error("get_potential: no potential %s found",name);
	dprintf (0,"[Compiling potential %s]\n",name);	
	sprintf (cmd,"make -f $NEMOLIB/Makefile.lib %s.so",name);
	dprintf (1,"%s\n",cmd);
	if (system(cmd)!=0)
	    error ("Error in compiling potential file");
	strcpy (name,fname);
	strcat (name,".so");
	fullname = name;	/* or use: findpath ? */
	loadobj (name);
    }

    /*
     * changed code 17/05/02 WD
     * debugged     18/09/08 WD
     *
     * after finding the .so file, we still need to find the proper 
     * routine to link with. This is controlled by the type.
     * if type = 'r',  behaviour depends on the macro SINGLEPREC:
     *                 If defined, we behave as for type='f', otherwise
     *                 as for type='d'.
     * if type = 'r',  we look for "potential" and restore to 
     *                 "potential_double" if no "potential" found
     * if type = 'd',  we first look for "potential_double" and
     *                 restore to "potential" if no "potential_double" found
     * if type = 'f',  we only look for "potential_float"
     *
     * macro FIND(POTENTIAL) tries to find routine POTENTIAL          WD
     */
#define FIND(POTENTIAL) {							\
  strcpy(pname,POTENTIAL);                 /*   try potential           */	\
  mapsys(pname);								\
  pot = (proc) findfn (pname);             /*   try C-routine           */ 	\
  if (pot==NULL) {                         /*   IF not found          > */ 	\
    strcat(pname,"_");								\
    pot = (proc) findfn (pname);           /*     try F77-routine       */	\
    if (pot)                               /*     found!                */ 	\
      Qfortran = TRUE;	 	     /*       must be F77 then    */ 		\
  }                                        /*   <                       */ 	\
  if(pot) dprintf(1,"\"%s\" loaded from file \"%s\"\n",POTENTIAL,fullname);	\
}

    char search_type = type=='r'?
#ifdef SINGLEPREC
      'f'
#else
      'd'
#endif
      : type;

    if(search_type=='f') {                   /* > ELSE, type=r          > */
      FIND("potential_float");               /*   try "potential_float"   */
    } else if(search_type=='d') {            /* > ELSE, type=d          > */
      FIND("potential_double");              /*   try "potential_double"  */
      if( pot==NULL) {                       /*   IF not found          > */
	FIND("potential");                   /*   try "potential"         */
      }
    } else
      error("unknown data type '%c'\n",type);
#undef FIND
    /* it is perhaps possible that some fortran compilers will add __ for */
    /* routines that have embedded _ in their name.... we're not catching */
    /* those here yet !!!                                                 */
    /* e.g.  g77 options:  -fno-underscoring and -fno-second-underscore   */
    /* will fix this problem                                              */
    if (pot==NULL) {
      error("Could not find a suitable potential type %c in %s",type,fname);
      return NULL;
    }

    strcpy(pname,"inipotential");
    mapsys(pname);
    ini_pot = (proc) findfn (pname);              		/* C */
    if (ini_pot==NULL) {
        strcpy(pname,"ini_potential");
        mapsys(pname);
        ini_pot = (proc) findfn (pname);	        	/* C */
        if (ini_pot==NULL) {
            strcpy(pname,"inipotential_");
            mapsys(pname);
            ini_pot = (proc) findfn (pname);		        /* F77 */
            if (ini_pot==NULL) {
                strcpy(pname,"ini_potential_");
                mapsys(pname);
                ini_pot = (proc) findfn (pname);        	/* F77 */
            }
        }
    }
    if (ini_pot)
        if (!Qfortran)
            (*ini_pot)(&local_npar,local_par,dataname); 	/* C */
        else {
            if (dataname==NULL)
                (*ini_pot)(&local_npar,local_par,dataname,0);   /* F77 */
            else
                (*ini_pot)(&local_npar,local_par,dataname,strlen(dataname)); /* F77 */

        }
    else {
        printf ("Warning: inipotential(_) not present in %s", fname);
        printf (",default taken\n");
    }
    l_potential = pot;            /* save these two for later references */
    l_inipotential = ini_pot;
    if (local_npar > 0 && local_par[0] != local_omega) {
    	local_omega = local_par[0];
    	dprintf(1,"get_potential: modified omega=%g\n",local_omega);
    }
    if (pot==NULL) potential_dummy_for_c();    /* should never be called */
    return pot;
}
/* endof: potential.c */
 
/******   now some junk needed to force linker to load some extra ******/

#include <mathlinker.h>		/* extra math */
#include <fslinker.h>		/* extra routines from filestruct() */

void potential_dummy_for_c(void)
{
    double a,b,c;
    int spline();
    double bessi0(), bessk0(), bessi1(), bessk1();
    void get_atable();
    void read_image();

    error("potential.c: Cannot call dummy_for_c - included to fool linkers");
    (void) spline();
    (void) bessi0();
    (void) bessk0();
    (void) bessi1();
    (void) sqr(1.0);
    stropen("/dev/null","w");
#ifndef NO_IMAGE
    read_image();
#endif    
    get_atable();
}

#if defined(FORTRAN)
void potential_dummy_for_fortran()
{
  extern void zzzzzz_(void);
  zzzzzz_();			/* force loading of Fortran I/O */
}
#endif

void potential_dummy_for_fortran_math()
{
  extern void fmath_(void);
  fmath_();
}
