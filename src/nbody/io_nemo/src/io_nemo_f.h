/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Prototypes                                                       
+---------------------------------------------------------------- */

#ifndef IO_NEMO_F_H
#define IO_NEMO_F_H
#ifdef __cplusplus
extern "C" {
#endif

int IO_NEMO_F(char * , int *, int *, char * , ...);

int CLOSE_IO_NEMO_F(char * , int *);

#ifdef __cplusplus
}
#endif
#endif /* IO_NEMO_F_H */
/* ----------------------------------------------------------------
|  End of io_nemo_f.h                                              
+---------------------------------------------------------------- */
