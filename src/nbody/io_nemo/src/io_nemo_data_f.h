/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Shared data imported from "io_nemo_f.c"                          
+----------------------------------------------------------------- */

#ifndef IO_NEMO_DATA_F_H
#define IO_NEMO_DATA_F_H
/* snapshot variables */
 extern char * pos_f,
             * vel_f,
             * phase_f,
             * pot_f,
             * acc_f,
             * mass_f,
             * eps_f,
             * keys_f,
             * aux_f,
             * dens_f,
             * timu_f,
             * selt_f,
             * selp_f;

 extern int  * nbody_f;

#endif /* IO_NEMO_DATA_F_H */
/* -----------------------------------------------------------------
|  End of io_nemo_data_f.h                                          
+----------------------------------------------------------------- */
