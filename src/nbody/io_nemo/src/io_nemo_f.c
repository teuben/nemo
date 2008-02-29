/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2006                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Perform I/O operations on NEMO structure data from a Fortran     
|  program.                                                         
| ==================================================================
| history                                                           
|                *            *             *                       
| 15-Jun-95	 V1.0 : created                                  JCL
| 21-Jun-95	 V1.10: bugs fixes                               JCL
| 12-Dec-95	 V1.11: possibility to close file                JCL
| 11-Mar-96	 V1.12: acceleration I/O added                   JCL
| 04-Apr-97	 V1.13: generic real format                      JCL
| 07-Apr-97	 V1.14: manual created                           JCL
| 19-Jul-02	 V1.20: io_nemo/io_nemo_f unified                JCL
| 18-Mar-04	 V1.21: bugs fixed, softening added              JCL
| 03-Mar-05	 V1.30: code cleaning, valgrind mem/leak safe    JCL
| 24-Apr-06      V1.31: memory leak fixed                        JCL
| 19-Jun-06      V1.32: happy gfortran                           JCL
| 29-May-07      V1.42: handle snapshot with different #bodies   JCL
+----------------------------------------------------------------  */

#ifdef ABSOFT
#  define IO_NEMO_F io_nemo_f__
#  define CLOSE_IO_NEMO_F close_io_nemo_f__
#else
#  define IO_NEMO_F io_nemo_f_
#  define CLOSE_IO_NEMO_F close_io_nemo_f_
#endif
/* ----------------------------------------------------------------
|  Include files                                                   
+---------------------------------------------------------------- */
#include <stdinc.h>
#include <stdarg.h>
#include <history.h>

/* ----------------------------------------------------------------
|  Local include files                                             
+---------------------------------------------------------------- */
#include "io_init.h" 
#include "check_file.h"
#include "parameters.h"
#include "io_get_put_f.h"
#include "io_nemo_f.h"
#include "io_nemo_tools.h"
#include "maxio.h"

/* extern variables */
#include "flags_data.h"

/* ----------------------------------------------------------------
|  Shared variables                                                
+---------------------------------------------------------------- */

/* variables for reading */
extern char * io_in[MAXIO]; 
extern FILE * instr[MAXIO];
extern bool   read_one[MAXIO];

/* variables for writing */
extern char * io_out[MAXIO]; 
extern FILE * outstr[MAXIO];
extern bool   save_one[MAXIO];

/* variables used to keep track of history (EXPORTED) */
bool   set_history_f[MAXIO];
char * history_prog_f=NULL;
char * hist_file_f;        /* history file name */

/* variables to store max #bodies per file */
int   maxbodies[MAXIO];

/* io_nemo_data_f (EXPORTED) */
char 
  * pos_f  ,   /* position           */
  * vel_f  ,   /* velocity           */
  * phase_f,   /* phase              */
  * pot_f  ,   /* potential          */
  * acc_f  ,   /* acceleration       */
  * mass_f ,   /* mass               */
  * eps_f  ,   /* softening          */
  * keys_f ,   /* keys               */
  * aux_f  ,   /* aux                */
  * dens_f ,   /* dens               */
  * timu_f ,   /* time steps         */
  * selt_f ,   /* selected time      */
  * selp_f ;   /* selected particles */
int   
  * nbody_f;   /* nbody              */

/* ----------------------------------------------------------------
|  io_nemo_f_ :                                                    
|  Function called from a FORTRAN program to perform I/O operations
|  on NEMO datas structure.                                        
|                                                                  
+---------------------------------------------------------------- */
int IO_NEMO_F(char * iofile,
	      int  * lg, 
	      int  * size_array, 
	      char * param, ...)
{
  va_list pa;
  char * p, * field;
  int 
    code;            /* return value according to I/O performed */
  bool io_op=TRUE;   /* TRUE -> READ,   FALSE -> SAVE */
  int rtype;         /* rtype = 1 (float), or = 2 (double) */
  static bool first=TRUE;	

  /* correct the FORTRAN string */
  iofile = (char *) f_ch_to_c(iofile,*lg);

  /* initialize selections flags */
  init_flag_io();
	
  /* initialize I/O variables */
  if (first) {
    init_io_one(maxbodies,read_one, save_one, set_history_f, &history_prog_f, MAXIO);
    first = FALSE;
  }
	
  param = set_eos(param,'\\');
  param = set_eos(param,'#');

  /* initialize unamed parameters */
  va_start(pa, param);
  p = param;
	
  while (*p) {
    field = (char *) get_field(&p);
    /*fprintf(stderr,"Field = [%s]\n",field);*/
    switch (get_case(field)) {
		
    case 1  : N_io = 1;
      nbody_f= va_arg(pa, int *);
      break;
		
    case 2  : T_io = 1;
      timu_f = va_arg(pa, char *);
      break;
			
    case 3  : M_io = 1;
      mass_f  = va_arg(pa, char *);
      break;
			
    case 4  : X_io = 1;
      pos_f  = va_arg(pa, char *);
      break;
			
    case 5  : V_io = 1;
      vel_f  = va_arg(pa, char *);
      break;
			
    case 6  : P_io = 1;
      pot_f  = va_arg(pa, char *);
      break;
			
    case 7  : A_io = 1;
      acc_f  = va_arg(pa, char *);
      break;
			
    case 8  : K_io = 1;
      keys_f = va_arg(pa, char *);
      break;

    case 9  : XV_io = 1;
      phase_f= va_arg(pa, char *);
      break;

    case 10  : EPS_io = 1;
      eps_f= va_arg(pa, char *);
      break;

    case 12  : AUX_io = 1;
      aux_f = va_arg(pa, char *);
      break;

    case 13  : D_io = 1;
      dens_f = va_arg(pa, char *);
      break;

    case 57 : ST_io = 1;
      selt_f = va_arg(pa, char *);
      break;
			
    case 50 : F_dim = 1;
      break;
			
    case 51 : F_dim = 0;
      break; 
			
    case 52 : io_op = FALSE;  /* save the snapshot */
      break;
			
    case 53 : io_op = TRUE;   /* read the snapshot */
      break;
			
    case 54 : rtype = 1;      /* float dimension */
      break;
			
    case 55 : rtype = 2;      /* double dimension */
      break;
			
    case 56 : I_io = 1;       /* print information */
      break;

    case 58 : SP_io=1;	      /* selected particles*/
      selp_f   = va_arg(pa, char *);
      break;
			
    default  : fprintf(stderr,
		       "Parameter error ## [io_nemo_f] \"%s\" unknown\n",
		       field);
    exit(1);
    break;
    }
    free(field);
  }
	
  va_end(pa);
	
  if (io_op) { /* operation on snaphot is reading */
    code = get_data_select_f(iofile, size_array, rtype, io_in, read_one, instr, MAXIO);
  }
  else {
    if (!N_io) {
      fprintf(stderr,"Parameter error ## [io_nemo_f] param : \"%s\"\n",
	      param);
      fprintf(stderr,
	      "You must specify \"nbody\" in the field parameter for SAVE operation \n");
      exit(1);
    }
    code = put_data_select_f(iofile, size_array, rtype, io_out, save_one, outstr, MAXIO);
  }
  return code;
}

/* ----------------------------------------------------------------
|  close_io_nemo_f_ :                                              
|  Close the opening snapshot.                                     
|  Return 0 if the file was not open, otherwise 1                  
+---------------------------------------------------------------- */
int CLOSE_IO_NEMO_F(char * iofile,int * lg)
{
  int no_io,code;

  /* correct the FORTRAN string */
  iofile = (char *) f_ch_to_c(iofile,*lg);
	
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
    free((char *) io_in[no_io]);
    code=1;
  }

  reset_history();
  return code;
}
/* ----------------------------------------------------------------
|  End of io_nemo_f.c                                              
+---------------------------------------------------------------- */
