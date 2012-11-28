/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Perform I/O operations on NEMO structure data from a C program   
| ==================================================================
| history                                                           
|                *            *             *                       
| 15-Nov-96	 V1.0 : created                                  JCL
| 21-Feb-97	 V1.1 : memory optimisation                      JCL
| 16-Apr-97	 V1.11: manual created                           JCL
| 19-Jul-02	 V1.2 : io_nemo and io_nemo_f unified            JCL
| 18-Mar-04	 V1.21: bugs fixed, softening added              JCL
| 03-Mar-05	 V1.30: memory bugs fixed, nemo control bits     JCL
|                       added, valgrind mem/leak safe               
| 24-Apr-06      V1.31: memory leak fixed                        JCL
| 19-Jun-06      V1.32: happy gfortran                           JCL
| 29-May-07      V1.42: handle snapshot with different #bodies   JCL
| 29-Feb-07      V1.50: Aux and Dens field added                 JCL
+----------------------------------------------------------------  */

/* -----------------------------------------------------------------
|  Include files                                                    
+----------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>

#include <stdio.h>
#include <stdarg.h>

#include "io_init.h" 
#include "parameters.h"
#include "io_get_put.h"
#include "check_file.h"
#include "io_nemo.h"
#include "maxio.h"

/* extern variables */
#include "flags_data.h"

/* -----------------------------------------------------------------
|  Shared variables                                                 
+----------------------------------------------------------------- */

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

/* variables to store max #bodies per file */
int   maxbodies[MAXIO];

/* ------------------------------------------------------------------
|  reajust_ptr :                                                     
|  match for each pointers passing to "io_nemo" the pointers         
|  allocated.                                                        
+------------------------------------------------------------------ */ 
void reajust_ptr()
{
}

/* ------------------------------------------------------------------
|  io_nemo :                                                         
|  function called from a C program to perform I/O on a NEMO file.   
+------------------------------------------------------------------ */
int io_nemo(char * iofile,
            char * param, ...)
{
  va_list pa;      /* variable parameter list                 */


  t_ion_data * ion;

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
    init_io_one(maxbodies, read_one, save_one, set_history, &history_prog, MAXIO);
    first = FALSE;
  }

  /* allocate memory for a new data structure */
  ion = (t_ion_data *) malloc(sizeof(t_ion_data));
  if (!ion) {
    fprintf(stderr,"Unable to allocate memory of size [t_ion_data], aborting...\n");
    exit(1);
  }

  /* get the first parameter from the variable list */
  va_start(pa, param);
  p = param;    /* point to the first element of the variable list */

  while (*p) {

    field = (char *) get_field(&p); /* get the next element */
    switch (get_case(field)) {

    case 1 : N_io = 1;
      ion->nbody_p = va_arg(pa, int **);
      ion->nbody   = *(ion->nbody_p); /* match local and parameter pointer */
      break;

    case 2 : T_io = 1;
      ion->time_p  = va_arg(pa, char **);
      ion->timu    = *(ion->time_p); /* match local and parameter pointer */
      break;

    case 3 : M_io = 1;
      ion->mass_p  = va_arg(pa, char **);
      ion->mass    = *(ion->mass_p); /* match local and parameter pointer */
      break;

    case 4 : X_io = 1;
      ion->pos_p = va_arg(pa, char **);
      ion->pos   = *(ion->pos_p);
      break;

    case 5 : V_io = 1;
      ion->vel_p   = va_arg(pa, char **);
      ion->vel     = *(ion->vel_p); /* match local and parameter pointer */
      break; 

    case 6 : P_io = 1;
      ion->pot_p   = va_arg(pa, char **);
      ion->pot     = *(ion->pot_p); /* match local and parameter pointer */
      break;

    case 7 : A_io = 1;
      ion->acc_p   = va_arg(pa, char **);
      ion->acc     = *(ion->acc_p); /* match local and parameter pointer */
      break;

    case 8 : K_io = 1;
      ion->keys_p   = va_arg(pa, char **);
      ion->keys     = *(ion->keys_p); /* match local and parameter pointer */
      break;

    case 10 : EPS_io = 1;
      ion->eps_p   = va_arg(pa, char **);
      ion->eps     = *(ion->eps_p); /* match local and parameter pointer */
      break;

    case 11 : B_io = 1;
      ion->bits_p   = va_arg(pa, int **);
      ion->bits     = *(ion->bits_p); /* match local and parameter pointer */
      break;

    case 12 : AUX_io = 1;
      ion->aux_p   = va_arg(pa, char **);
      ion->aux     = *(ion->aux_p); /* match local and parameter pointer */
      break;

    case 13 : D_io = 1;
      ion->dens_p   = va_arg(pa, char **);
      ion->dens     = *(ion->dens_p); /* match local and parameter pointer */
      break;

    case 57 : ST_io = 1;
      ion->selt     = va_arg(pa, char *);
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
      ion->SelectionString123=va_arg(pa, char*);
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
  else {
    if (io_op) { /* operation on snaphot is reading */
      code = get_data_select(iofile, rtype, io_in, read_one, instr, MAXIO,ion);
      /* reajust pointers */
      /*reajust_ptr();*/
      if (N_io)    *(ion->nbody_p)    = ion->nbody;
      if (T_io)    *(ion->time_p )    = ion->timu;
      if (M_io)    *(ion->mass_p )    = ion->mass;
      if (X_io)    *(ion->pos_p  )    = ion->pos;
      if (V_io)    *(ion->vel_p  )    = ion->vel;
      if (XV_io)   *(ion->phase_p)    = ion->phase;
      if (AUX_io)  *(ion->aux_p)      = ion->aux;
      if (D_io)    *(ion->dens_p)     = ion->dens;
      if (P_io)    *(ion->pot_p  )    = ion->pot;
      if (A_io)    *(ion->acc_p  )    = ion->acc; 
      if (K_io)    *(ion->keys_p )    = ion->keys; 
      if (EPS_io)  *(ion->eps_p  )    = ion->eps; 
      if (B_io)    *(ion->bits_p )    = ion->bits; 
    }
    else { 
      if (!N_io) {
	fprintf(stderr,"Parameter error ## [io_nemo] param : \"%s\"\n",
		param);
	fprintf(stderr,
		"You must specify \"nbody\" "
		"in the field parameter for SAVE operation \n");
	exit(1);
      }
      code = put_data_select(iofile, rtype, io_out, save_one, outstr, MAXIO,ion);
    }
  }
  free ((t_ion_data *) ion);
  return code;
}
/* ------------------------------------------------------------------
|  close_io_nemo :                                                   
|  Close the opening snapshot.                                       
|  Return 0 if the file was not open, otherwise 1                    
+------------------------------------------------------------------ */
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
    maxbodies[no_io]= 0;
    set_history[no_io] = FALSE;
    free((char *) io_in[no_io]);
    code=1;
  }
  /* reset_history();  */
  return code;
}
/* ------------------------------------------------------------------
|  End of io_nemo.c                                                  
+------------------------------------------------------------------ */
