/* -------------------------------------------------------------- *\
|* $Id$
|*
|* Perform I/O operations on NEMO structure data from a C program
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

#include <stdarg.h>

#include "io_init.h" 
#include "parameters.h"
#include "io_get_put.h"
#include "check_file.h"
#include "io_nemo.h"
#include "maxio.h"

/* extern variables */
#include "flags_data.h"

/* -------------------------------------------------------------- *\
|* Shared variables
\* -------------------------------------------------------------- */

/* Variables used for reading operations */
char * io_in[MAXIO]; 
FILE * instr[MAXIO];
bool   read_one[MAXIO];

/* Variables used for writing operations */
char * io_out[MAXIO]; 
FILE * outstr[MAXIO];
bool   save_one[MAXIO];

/* variables used to keep track of history (EXPORTED) */
bool   set_history[MAXIO];
char * history_prog=NULL;
char * hist_file;         /* history file name  */

/* io_nemo_data (EXPORTED) */
char         
  * pos,   ** pos_p,      /* position           */
  * vel,   ** vel_p,      /* velocity           */
  * phase, ** phase_p,    /* phase              */
  * pot,   ** pot_p,      /* potential          */
  * acc,   ** acc_p,      /* acceleration       */
  * mass,  ** mass_p,     /* mass               */
  * keys,  ** keys_p,     /* keys               */
  * timu,  ** time_p,     /* time steps         */
  * selt,  ** selt_p,     /* selected time      */
  *SelectionString123;    /* selected particles */
int   
  * nbody, ** nbody_p;    /* nbody              */

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
  if (XV_io)   *phase_p    = phase;
  if (P_io)    *pot_p      = pot;
  if (A_io)    *acc_p      = acc; 
  if (K_io)    *keys_p     = keys; 
  /*  if (ST_io)   *selt_p     = selt; */
}

/* -------------------------------------------------------------- *\ 
|* io_nemo :
|* function called from a C program to perform I/O on a NEMO file.
\* -------------------------------------------------------------- */
int io_nemo(char * iofile,
            char * param, ...)
{
  va_list pa;      /* variable parameter list                 */

  char * p, * field;
  int code;        /* return value according to I/O performed */
  bool io_op=TRUE; /* TRUE -> READ,   FALSE -> SAVE           */
  int rtype;       /* rtype = 1 (float), or = 2 (double)      */
  /* control */
  static bool first=TRUE;


  /* init flag */
  init_flag_io();

  /* init I/O variable */
  if (first) {
    init_io_one(read_one, save_one, set_history, &history_prog, MAXIO);
    first = FALSE;
  }

  /* get the first parameter from the variable list */
  va_start(pa, param);
  p = param;    /* point to the first element of the variable list */

  while (*p) {

    field = (char *) get_field(&p); /* get the next element */
    switch (get_case(field)) {

    case 1 : N_io = 1;
      nbody_p = va_arg(pa, int **);
      nbody   = *nbody_p; /* match local and parameter pointer */
      break;

    case 2 : T_io = 1;
      time_p  = va_arg(pa, char **);
      timu    = *time_p; /* match local and parameter pointer */
      break;

    case 3 : M_io = 1;
      mass_p  = va_arg(pa, char **);
      mass    = *mass_p; /* match local and parameter pointer */
      break;

    case 4 : X_io = 1;
      pos_p   = va_arg(pa, char **);
      pos     = *pos_p; /* match local and parameter pointer */
      break;

    case 5 : V_io = 1;
      vel_p   = va_arg(pa, char **);
      vel     = *vel_p; /* match local and parameter pointer */
      break; 

    case 6 : P_io = 1;
      pot_p   = va_arg(pa, char **);
      pot     = *pot_p; /* match local and parameter pointer */
      break;

    case 7 : A_io = 1;
      acc_p   = va_arg(pa, char **);
      acc     = *acc_p; /* match local and parameter pointer */
      break;

    case 8 : K_io = 1;
      keys_p   = va_arg(pa, char **);
      keys     = *keys_p; /* match local and parameter pointer */
      break;

    case 57 : ST_io = 1;
      selt     = va_arg(pa, char *);
      break;

    case 52 : io_op = FALSE;  /* save the snapshot */
      break;

    case 53 : io_op = TRUE;   /* read the snapshot */
      break;

    case 54 : rtype = 1;     /* float dimension */
      break;

    case 55 : rtype = 2;     /* double dimension */
      break;

    case 56 : I_io = 1;      /* print information */
      break;

    case 59 : H_io = 1;      /* History file name */
      hist_file = va_arg(pa, char *);
      break;
	  
    case 58 : SP_io=1;
      SelectionString123=va_arg(pa, char*);
      break;

    case 60 : C_io = 1;
      break;

    default  : fprintf(stderr,
		       "Parameter error ## [io_nemo] \"%s\" unknown\n",
		       field);
    exit(1);
    break;
    }

    free(field);
  }

  va_end(pa);

  if (C_io)  /* close the snapshot */
    code=close_io_nemo(iofile);
  else
    if (io_op) { /* operation on snaphot is reading */
      code = get_data_select(iofile, rtype, io_in, read_one, instr, MAXIO);
      /* reajust pointers */
      reajust_ptr();
    }
    else
      { 
	if (!N_io) {
	  fprintf(stderr,"Parameter error ## [io_nemo] param : \"%s\"\n",
		  param);
	  fprintf(stderr,
		  "You must specify \"nbody\" "
		  "in the field parameter for SAVE operation \n");
	  exit(1);
	}
	code = put_data_select(iofile, rtype, io_out, save_one, outstr, MAXIO);
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
  if ((no_io = get_old_file(iofile,io_in,read_one,instr,MAXIO)) < 0) {
    if ((no_io = get_old_file(iofile,io_out,save_one,outstr,MAXIO)) < 0) {
      fprintf(stderr,
	      "WARNING!! snapshot [%s] not OPEN, unable to close it\n",
	      iofile);
      code=0;
    } 
    else {
      /* close the file open for writing */
      strclose(outstr[no_io]);

      /* RAZ variables */
      save_one[no_io] = FALSE;
      set_history[no_io] = FALSE;
	  
      free((char *) io_out[no_io]);
      code=1;
    }
  }
  else {
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
