/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@lam.fr                           
|  address:  Centre de donneeS Astrophysique de Marseille (CeSAM)   
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 7326                                        
| ==================================================================
|* Get/Put selected Data from io_nemo_f                             
| ==================================================================
| 03-Mar-05 :                                                       
+----------------------------------------------------------------- */
#ifndef IO_GET_PUT_F_H
#define IO_GET_PUT_F_H
#ifdef __cplusplus
extern "C" {
#endif

int put_data_select_f(char * ,
		      int  * ,
		      int    ,
		      char * a[],
		      bool *    ,
		      FILE * c[],
		      int     );

int get_data_select_f(char * ,
		      int  * , 
		      int    ,
		      char * a[],
		      bool *    ,
		      FILE * c[],
		      int     );


#ifdef __cplusplus
}
#endif

#endif  /* IO_GET_PUT_F_H */
/* ----------------------------------------------------------------
|  End of io_get_put_f.h                                           
+---------------------------------------------------------------- */
