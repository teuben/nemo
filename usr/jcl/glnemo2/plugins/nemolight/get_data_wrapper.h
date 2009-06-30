/* ================================================================
|  Copyright Jean-Charles LAMBERT - 2008                           
|  e-mail:   Jean-Charles.Lambert@oamp.fr                          
|  address:  Dynamique des galaxies                                
|            Laboratoire d'Astrophysique de Marseille              
|            2, place Le Verrier                                   
|            13248 Marseille Cedex 4, France                       
|            CNRS U.M.R 6110                                       
| =================================================================
|* Wrapper of basic NEMO procedure, it makes me life easiest :)    
+---------------------------------------------------------------- */
#ifndef GET_DATA_WRAPPER_H
#define GET_DATA_WRAPPER_H
#ifdef __cplusplus
extern "C" {
#endif

int get_data_gen(stream, char *, char *, int, int, int, int, char **);

int get_data_time(stream, char *, int, char **);

int get_data_nbody(stream, char *, int, int **);

int get_data_mass(stream, char *, int, int, char **);

int get_data_keys(stream, char *, int, int, char **);

int get_data_phase(stream, char *, int, int, char **, int);

int get_data_pos(stream, char *, int, int, char **, int);

int get_data_vel(stream, char *, int, int, char **, int);

int get_data_pot(stream, char *, int, int, char **);

int get_data_acc(stream, char *, int, int, char **, int);

int get_data_eps(stream, char *, int, int, char **);

int get_data_aux(stream, char *, int, int, char **);

int get_data_dens(stream, char *, int, int, char **);
#ifdef __cplusplus
}
#endif
#endif /* GET_DATA_WRAPPER_H */
/* ----------------------------------------------------------------
|  End of [get_dat_nemo.h]                                         
+---------------------------------------------------------------- */ 
