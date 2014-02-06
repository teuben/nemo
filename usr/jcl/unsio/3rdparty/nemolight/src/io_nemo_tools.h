/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Functions prototypes                                             
+----------------------------------------------------------------  */

#ifndef IO_NEMO_TOOLS_H
#define IO_NEMO_TOOLS_H
#ifdef __cplusplus
extern "C" {
#endif

float   char2float(char * ,int);

double  char2double(char * ,int );

char *  allocate_pointer(char *, int );

char *  f_ch_to_c(char * ,int );

bool ** chk_select(int * ,int ,int ,string select_pts[]);

char *  get_selected(char *);
  char *  set_eos(char *,char);

#ifdef __cplusplus
}
#endif

#endif /* IO_NEMO_TOOLS_H */
/* ----------------------------------------------------------------
|  End of io_nemo_tools.h                                          
+---------------------------------------------------------------- */
