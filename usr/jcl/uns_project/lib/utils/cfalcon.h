// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/
#ifndef CFALCON_H
#define CFALCON_H

namespace jclut {
class cfalcon
{
public:
    cfalcon();
    static bool addGravity(const int nbody,
                           const float * pos, const float * mass, 
                           float * acc, float * phi,
                           const float eps, 
                           const float G=1.0,
                           const float theta=0.6,
                           const int kernel_type=1,
                           const int ncrit=6);
    
    static bool addGravity2(const int nbody,
                             const float * pos, const float * mass, 
                             const int nbody_tp,
                             const float * pos_tp,
                             float * acc, float * phi,
                             const bool selfp,
                             const float eps,
                             const float G=1.0,
                             const float theta=0.6,
                             const int kernel_type=1,
                             const int ncrit=6);

private:
};
}
#endif // CFALCON_H
