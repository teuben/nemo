/* ================================================================
|  Copyright Jean-Charles LAMBERT - 2005                           
|  e-mail:   Jean-Charles.Lambert@oamp.fr                          
|  address:  Dynamique des galaxies                                
|            Laboratoire d'Astrophysique de Marseille              
|            2, place Le Verrier                                   
|            13248 Marseille Cedex 4, France                       
|            CNRS U.M.R 6110                                       
| =================================================================
|* Keep track of open/close files                                  
+---------------------------------------------------------------- */

/* ----------------------------------------------------------------
|  Nemo include files                                              
+---------------------------------------------------------------- */
#include <stdinc.h>
#include <history.h>

/* ----------------------------------------------------------------
|  Local include files                                             
+---------------------------------------------------------------- */
#include "check_file.h"

/* ----------------------------------------------------------------
|  get_history_input_file :                                        
|  get the history stored in the 'history_file'                    
+---------------------------------------------------------------- */
int get_history_input_file(char * hist_file)
{
  stream instr;
  
  instr=stropen(hist_file,"r");
  get_history(instr);
  strclose(instr);
  return 0;
}
/* -----------------------------------------------------------------
|  get_old_file :                                                   
|  Check if the filename is already open. If yes return the logical 
|  name otherwise '-1'.                                             
+----------------------------------------------------------------- */
int get_old_file(char * curr_file,
		 char * io_file[],
		 bool   io_one[],
		 FILE * io_str[],
		 int MAXIO)
{ 
  int i;
  
  for (i=0; i < MAXIO; i++) {
    if (io_one[i]) {
      if (!strcmp(curr_file, io_file[i]))
	return i;
    }
  }
  return -1;
}
/* -----------------------------------------------------------------
|  get_new_file :                                                   
|  Open a file either for reading or either for writing.            
+----------------------------------------------------------------- */
int get_new_file(char * curr_file,
		 char * io_file[],
		 bool   io_one[],
		 FILE * io_str[],
		 char * mode,
		 int MAXIO)
{ 
  int i;

  for (i=0; i < MAXIO; i++) {
    if (!io_one[i]) {
      io_file[i] = (char *) malloc(strlen(curr_file)*sizeof(char)+1);
      
      if (io_file[i]==NULL) {
	fprintf(stderr,"Memory error ## [get_new_file]\n");
	fprintf(stderr,"Impossible to allocate memory\n");
	exit(1);
      }
      
      strcpy(io_file[i],curr_file);
      io_str[i] = stropen(io_file[i],mode);

      if (io_str[i]==NULL) {
	fprintf(stderr,"I/O error ## [get_new_file]\n");
	fprintf(stderr,"File \"%s\" open error in mode \"%s\"\n",
		io_file[i],mode);
	exit(1);
      }
      return i;
    }
  }
  fprintf(stderr,"Error!! ## MAXIO number ## [get_new_file]\n");
  fprintf(stderr,"number MAXIO=(%d) reached, too much FILES open\n",MAXIO);
  exit(1);
}
/* ----------------------------------------------------------------
|  End of check_file.c                                             
+---------------------------------------------------------------- */
