/* -------------------------------------------------------------- *\
|* io_nemo.c :	      JC LAMBERT	19-Nov-96	V1.3
|*					
|*
|* V1.0  : 19-Nov-96	created	
|* V1.1  : 21-Feb-97	memory optimisation,bugs fixed
|* V1.2  : 20-Sep-97    save history, added
|* V1.3  : 17-Mar-98    selected particles, added    
|*
|* Perform I/O NEMO for C programs.
|*
|* Materiel     : Sparc STATION
|* OS           : Solaris 2.x
|* Langage      : C
|* Version Nemo : 2.0 Juillet 1994 
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>

#include <string.h>
#include <stdio.h>
#include <stdarg.h>

#include "maxio.h"
#include "io_init.h" 
#include "parameters.h"
#include "io_get_put.h"
#include "check_file.h"
#include "io_nemo.h"

/* -------------------------------------------------------------- *\
|* Shared variables
\* -------------------------------------------------------------- */

/* >>> "file_data.h" */
/* Variables used for reading operations */
char * io_in[MAXIO]; 
FILE * instr[MAXIO];
bool   read_one[MAXIO];

/* Variables used for writing operations */
char * io_out[MAXIO]; 
FILE * outstr[MAXIO];
bool   save_one[MAXIO];

/* variables used to keep track of history */
bool set_history[MAXIO];
char * history_prog=NULL;
char * hist_file;        /* history file name */

/* <<< */

/* >>> "snapshot_data.h" */
/* Snapshot datas */
char         * pos, ** pos_p,    /* position */
             * vel, ** vel_p,    /* velocity */
             * pot, ** pot_p,    /* potential */
             * acc, ** acc_p,    /* acceleration */
             * mass,** mass_p,   /* mass */
             * timu,** time_p,   /* time steps */
             * selt,** selt_p,   /* selected time */
			 *SelectionString123;
             
int   * nbody,** nbody_p;   /* nbody */
/* <<< */

/* >>> "flags_data.h" */
/* flags parameters */
int N_io, T_io, M_io, X_io, V_io, P_io, A_io,
           C_io, S_io, R_io, I_io, H_io, ST_io, SP_io;
/* <<< */

/* control */
bool first=TRUE;

/* -------------------------------------------------------------- *\ 
|* allocate_pointer :
|* check if the pointer is already allocated. if not it is
|* allocated.
\* -------------------------------------------------------------- */ 
char * allocate_pointer(char * p, int lg)
{ char * ptr;

  if (p == NULL)
    { 
      ptr = (char *) malloc(sizeof(char)*lg);
      if (!ptr)
	{ 
	  fprintf(stderr,
		  "[allocate_pointer], allocation memory error, aborted\n");
          exit(1);
        }
      return ptr;
    }
  else /* assume p is already  allocated */
      return p;
} 
/* -------------------------------------------------------------- *\ 
|* reajust_ptr :
|* match for each pointers passing to "io_nemo" the pointers
|* allocated.
\* -------------------------------------------------------------- */ 
void reajust_ptr()
{ 
  if (N_io)    *nbody_p    = nbody;
  if (T_io)    *time_p     = timu;
  if (M_io)    *mass_p     = mass;
  if (X_io)    *pos_p      = pos;
  if (V_io)    *vel_p      = vel;
  if (P_io)    *pot_p      = pot;
  if (A_io)    *acc_p      = acc; 
  /*  if (ST_io)   *selt_p     = selt; */
}

/* -------------------------------------------------------------- *\ 
|* io_nemo :
|* function called from a C program to perform I/O on a NEMO file.
\* -------------------------------------------------------------- */
int io_nemo(char * iofile,
            char * param, ...)
{
  va_list pa;

  char * p, * champs;
  int code; /* return value, depend the operation realized */
  bool io_op=TRUE; /* TRUE -> READ,   FALSE -> SAVE */
  int rtype; /* rtype = 1 (float), or = 2 (double) */


  /* init flag */
  init_flag_io();

  /* init I/O variable */
  if (first)
    { 
      init_io_one();
      first = FALSE;
    }

  /* information flag */
  I_io = 0;

  /* get the first parameter from the generic list */
  va_start(pa, param);
  p = param;

#ifdef DEBUG
 fprintf(stderr,"Param >%s<\n",param);
#endif

  while (*p)
     { 
       champs = (char *) get_champs(&p);
#ifdef DEBUG
       fprintf(stderr,"Champs >%s<\n",champs);
       fprintf(stderr,"get_case(%s)=%d\n",champs,get_case(champs));
#endif
       switch (get_case(champs))
        { case '\001' : N_io = 1;
                     nbody_p = va_arg(pa, int **);
		     nbody   = *nbody_p; /* match local and parameter pointer */
                     break;
 
          case '\002' : T_io = 1;
                     time_p  = va_arg(pa, char **);
		     timu    = *time_p; /* match local and parameter pointer */
                     break;

          case '\003' : M_io = 1;
                     mass_p  = va_arg(pa, char **);
		     mass    = *mass_p; /* match local and parameter pointer */
                     break;

          case '\004' : X_io = 1;
                     pos_p   = va_arg(pa, char **);
		     pos     = *pos_p; /* match local and parameter pointer */
		     break;

          case '\005' : V_io = 1;
                     vel_p   = va_arg(pa, char **);
		     vel     = *vel_p; /* match local and parameter pointer */
		     break; 

          case '\006' : P_io = 1;
		     pot_p   = va_arg(pa, char **);
		     pot     = *pot_p; /* match local and parameter pointer */
                     break;

          case '\012' : A_io = 1;
		     acc_p   = va_arg(pa, char **);
		     acc     = *acc_p; /* match local and parameter pointer */
                     break;

          case '\017' : ST_io = 1;
		     selt     = va_arg(pa, char *);
                     break;

          case '\007' : C_io = 1;
                     break;

          case '\010' : io_op = FALSE;  /* save the snapshot */
                     break;

          case '\011' : io_op = TRUE;   /* read the snapshot */
                     break;

          case '\013' : rtype = 1;     /* float dimension */
                     break;

          case '\014' : rtype = 2;     /* double dimension */
                     break;

          case '\015' : I_io = 1;      /* print information */
                     break;

          case '\016' : H_io = 1;      /* History file name */
	             hist_file = va_arg(pa, char *);
                     break;
	  
	  case '\020' : SP_io=1;
		     SelectionString123=va_arg(pa, char*);
		     break;

          default  : fprintf(stderr,
                    "Parameter error ## [io_snap_nemo] \"%s\" unknown\n",
                     champs);
                     exit(1);
                     break;
        }

       free(champs);
     }

  va_end(pa);

#ifdef DEBUG
    fprintf(stderr,"Debut de get_data_select\n");
#endif
  if (C_io)  /* close the snapshot */
    code=close_io_nemo(iofile);
  else
    if (io_op) /* operation on snaphot is reading */
      { 
	code = get_data_select(iofile, rtype);
	/* reajust pointers */
	reajust_ptr();
      }
    else
      { 
	if (!N_io)
	  { 
	    fprintf(stderr,"Parameter error ## [io_nemo] param : \"%s\"\n",
		    param);
	    fprintf(stderr,"You must specify \"nbody\" in the field parameter for SAVE operation \n");
	    exit(1);
	  }
	code = put_data_select(iofile, rtype);
      }
  
  return code;
}

/* -------------------------------------------------------------- *\ 
|* close_io_nemo :
|* Close the opening snapshot.
|* Return 0 if the file was not open, otherwise 1
\* -------------------------------------------------------------- */
int close_io_nemo(char * iofile)
{
  int no_io,code;
  
  /* check if the file is already open */

  if ((no_io = get_old_file(iofile,io_in,read_one,instr)) < 0) 
    { 
      if ((no_io = get_old_file(iofile,io_out,save_one,outstr)) < 0)
	{ 
	  fprintf(stderr,"WARNING!! snapshot [%s] not OPEN, unable to close it\n",iofile);
	  code=0;
	} 
      else
	{ 
	  /* close the file open for writing */
	  strclose(outstr[no_io]);

	  /* RAZ variables */
	  save_one[no_io] = FALSE;
	  set_history[no_io] = FALSE;
	  
	  free((char *) io_out[no_io]);
	  code=1;
	}
    }
  else
    { 
      /* close the file open for reading */
      strclose(instr[no_io]);

      /* RAZ variables */
      read_one[no_io] = FALSE;
      free((char *) io_in[no_io]);
      code=1;
    }

  return code;
}
/* -------------------------------------------------------------- *\ 
|* End of io_nemo.c 
\* -------------------------------------------------------------- */
