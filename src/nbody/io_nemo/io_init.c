/* -------------------------------------------------------------- *\
|* io_init.c :	JC LAMBERT	 09-Feb-98            V1.5
|*
|* V1.0  : created
|*
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Nemo include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <string.h>
#include <history.h>

/* -------------------------------------------------------------- *\
|* Local include files
\* -------------------------------------------------------------- */
#include "maxio.h"
#include "file_data.h"
#include "flags_data.h"

#include "io_nemo.h"

/* -------------------------------------------------------------- *\ 
|* init_io_one :
|* Initilalize the NEMO engine and some variables.
\* -------------------------------------------------------------- */
void init_io_one()
{ 
  string defv[] = { "none=none","VERSION=1.3",NULL };
  string argv[] = { "IO_NEMO",NULL };
  int i;
  string * histo;

  initparam(argv,defv); 
  
  for (i=0; i< MAXIO; i++)
     { 
       read_one[i] = FALSE;
       save_one[i] = FALSE;
       set_history[i] = FALSE;
     }

  /* get command line history */
  histo = (char **) ask_history();
  history_prog = (char *) allocate_pointer(history_prog,
					   sizeof(char)*strlen(histo[0]));

  /* save command line history */
  strcpy(history_prog, histo[0]);

  /*fprintf(stderr,"history prog : <%s> \n", history_prog); */
}

/* -------------------------------------------------------------- *\ 
|* init_flag_io : 
|* Initialise les drapeaux de recuperation de parametres
\* -------------------------------------------------------------- */
void init_flag_io()
{
  N_io= T_io= M_io= X_io= V_io= P_io= A_io= S_io= R_io = C_io = H_io=0;
  ST_io=0;
  SP_io=0;
}

/* -------------------------------------------------------------- *\
|* End of io_init.c
\* -------------------------------------------------------------- */
