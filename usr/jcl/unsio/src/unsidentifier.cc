// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2014
//           Centre de donneeS Astrophysiques de Marseille (CeSAM)              
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Aix Marseille Universite, CNRS, LAM 
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS UMR 7326                                       
// ============================================================================

#include "unsidentifier.h"

namespace uns {
// ============================================================================
// getIdentUNSV                                                                
  int CunsIdentifier::getUnsvIndex(const int ident,CunsIdentifierVector * unsv)
{
  for (CunsIdentifierVector::iterator it=unsv->begin(); it!= unsv->end(); it++) {
    if (ident == it->ident) {
      return (it-unsv->begin()); // return the valid identifier
    }
  }
  return -1;
}
}
