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

#include <iostream>
#include <fstream>                                    // C++ file I/O
#include <cmath>
#include <sstream>
#include "csimtools.h"

using namespace jclut;
using namespace uns;
// ============================================================================
// constructor
CSimtools::CSimtools(CunsIn * _uns)
{
  unsin = _uns;
}
// ============================================================================
// loadCod
void CSimtools::loadCod()
{
  std::ifstream         // File Handler
      fd;               // manipfile file desc

  std::string codfile=unsin->snapshot->getSimDir()+
                      "ANALYSIS/"+
                      unsin->snapshot->getFileName()+".";
}

// ============================================================================

// ============================================================================

// ============================================================================
