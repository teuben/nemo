/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Get/Put selected Data from io_nemo                               
+----------------------------------------------------------------- */
#ifndef IO_GET_PUT_H
#define IO_GET_PUT_H
#include "io_nemo.h"

#ifdef __cplusplus
extern "C" {
#endif

int put_data_select(char * ,
		    int    ,
		    char * a[],
		    bool *    ,
                    FILE * c[],
		    int,
		    t_ion_data * ion);

int get_data_select(char * ,
		    int    ,
		    char * a[],
		    bool *    ,
                    FILE * c[],
		    int       ,
		    t_ion_data * ion);


#ifdef __cplusplus
}
#endif
#endif /* IO_GET_PUT_H */
/* ----------------------------------------------------------------
|  End of io_get_put.h                                             
+---------------------------------------------------------------- */
