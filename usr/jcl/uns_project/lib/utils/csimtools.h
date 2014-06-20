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
#ifndef CSIMTOOLS_H
#define CSIMTOOLS_H
#include "uns.h"
using namespace uns;
namespace jclut {
  
  class CSimtools
  {
  public:
    CSimtools(CunsIn *);
    
    
  private:
    CunsIn * unsin;
    
    void loadCod();
  };
}
#endif // CSIMTOOLS_H
