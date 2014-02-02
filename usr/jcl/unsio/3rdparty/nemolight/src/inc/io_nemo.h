/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@lam.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS UMR 7326                            
| ==================================================================
|* Functions prototypes                                             
+----------------------------------------------------------------  */

#ifndef IO_NEMO_H
#define IO_NEMO_H
typedef struct {
  char         
    * pos,   ** pos_p,      /* position           */
    * vel,   ** vel_p,      /* velocity           */
    * phase, ** phase_p,    /* phase              */
    * pot,   ** pot_p,      /* potential          */
    * acc,   ** acc_p,      /* acceleration       */
    * mass,  ** mass_p,     /* mass               */
    * aux ,  ** aux_p,      /* aux                */
    * dens,  ** dens_p,     /* density            */
    * keys,  ** keys_p,     /* keys               */
    * eps,   ** eps_p,      /* softening          */
    * timu,  ** time_p,     /* time steps         */
    * selt,  ** selt_p,     /* selected time      */  
    *SelectionString123;    /* selected particles */

  int 
    * nbody, ** nbody_p,    /* nbody              */
    * bits,  ** bits_p;     /* control bits       */  
} t_ion_data;

#ifdef __cplusplus
extern "C" {
#endif

void reajust_ptr();

int io_nemo(char * , char * , ...);

int close_io_nemo(char * );

#ifdef __cplusplus
}
#endif

#endif /* IO_NEMO_H */
/* ----------------------------------------------------------------
|  End of io_nemo.h                                                
+---------------------------------------------------------------- */
