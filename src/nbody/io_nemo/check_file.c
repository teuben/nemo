/* -------------------------------------------------------------- *\
|* check_file.c :	JC LAMBERT	 20-Oct-97	V1.0
|*
|* V1.0  : created
|*
|* 
|*
|* Materiel     : Sparc STATION
|* OS           : Solaris 2.x
|* Langage      : C
|* Version Nemo : 2.0 Juillet 1994 
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Nemo include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <history.h>
#include <string.h>

/* -------------------------------------------------------------- *\
|* Local include files
\* -------------------------------------------------------------- */
#include "maxio.h"
#include "file_data.h"


/* -------------------------------------------------------------- *\ 
|* get_history_input_file :
|* get the history stored in the 'history_file'
\* -------------------------------------------------------------- */
int get_history_input_file()
{
  stream instr;
  
  instr=stropen(hist_file,"r");
  get_history(instr);
  strclose(instr);
}
/* -------------------------------------------------------------- *\ 
|* get_old_file :
|* Check if the filename is already open. If yes return the logical
|* name otherwise '-1'.
\* -------------------------------------------------------------- */
int get_old_file(char * curr_file,
		 char * io_file[],
		 bool   io_one[],
		 FILE * io_str[])
{ 
  int i;
  
  for (i=0; i < MAXIO; i++)
    if (io_one[i])
      {
        if (!strcmp(curr_file, io_file[i]))
            return(i);
      }

  
  return -1;
}
/* -------------------------------------------------------------- *\ 
|* get_new_file :
|* Open a file either for reading or either for writing.
\* -------------------------------------------------------------- */
int get_new_file(char * curr_file,
		 char * io_file[],
		 bool   io_one[],
		 FILE * io_str[],
		 char * mode)
{ 
  int i;

  for (i=0; i < MAXIO; i++)
    if (!io_one[i])
      { 

	io_file[i] = (char *) malloc(strlen(curr_file)*sizeof(char)+1);

        if (io_file[i]==NULL)
          { 
	    fprintf(stderr,"Memory error ## [get_new_file]\n");
            fprintf(stderr,"Impossible to allocate memory\n");
            exit(1);
          }

        strcpy(io_file[i],curr_file);

#ifdef DEBUG
	fprintf(stderr,"io_file[%d] = <%s> io_str[%d] = %x, mode=[%s]\n",
		i,io_file[i],i,io_str[i],mode);
#endif
        io_str[i] = stropen(io_file[i],mode);

#ifdef DEBUG
	fprintf(stderr,"io_file[%d] = <%s> io_str[%d] = %x, mode=[%s]\n",
		i,io_file[i],i,io_str[i],mode);
#endif

        if (io_str[i]==NULL)
          { 
	    fprintf(stderr,"I/O error ## [get_new_file]\n");
            fprintf(stderr,"File \"%s\" open error in mode \"%s\"\n",
                    io_file[i],mode);
            exit(1);
          }

        return i;
      }

   fprintf(stderr,"Error!! ## MAXIO number ## [get_new_file]\n");
   fprintf(stderr,"number MAXIO=(%d) reached, too much FILES open\n",MAXIO);
   exit(1);
}


/* -------------------------------------------------------------- *\
|* End of check_file.c
\* -------------------------------------------------------------- */
