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
#ifndef CHECK_FILE_H
#define CHECK_FILE_H
#ifdef __cplusplus
extern "C" {
#endif

int get_history_input_file(char *);

int get_old_file(char * curr_file ,
		 char * io_file[] ,
		 bool   io_one[]  ,
		 FILE * io_str[]  ,
		 int);

int get_new_file(char * curr_file ,
		 char * io_file[] ,
		 bool   io_one[]  ,
		 FILE * io_str[]  ,
		 char * mode,
		 int);

#ifdef __cplusplus
}
#endif

#endif /* CHECK_FILE_H */
/* ----------------------------------------------------------------
|  End of check_file.h                                             
+---------------------------------------------------------------- */
