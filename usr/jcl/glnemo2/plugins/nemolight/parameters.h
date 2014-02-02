/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@lam.fr                           
|  address:  Centre de donneeS Astrophysique de Marseille (CeSAM)   
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 7326                                        
| ==================================================================
|* Functions prototypes                                             
+----------------------------------------------------------------  */

#ifndef PARAMETERS_H
#define PARAMETERS_H
#ifdef __cplusplus
extern "C" {
#endif

int get_case(char * );

char * get_field(char **);

int chk_parameters(bool ,int , int );


#ifdef __cplusplus
}
#endif

#endif /* PARAMETERS_H */
/* ----------------------------------------------------------------
|  End of parameters.h                                             
+---------------------------------------------------------------- */
