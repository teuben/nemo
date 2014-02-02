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
#ifndef UNSIDENTIFIER_H
#define UNSIDENTIFIER_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/

#include <vector>
#include "uns.h"

namespace uns {

  // class CunsIdentifier                                 
  // store identifier index and pointer to UNS object     
  //                                                      
  class CunsIdentifier;
  typedef std::vector <CunsIdentifier> CunsIdentifierVector;
 
  class CunsIdentifier {
  public:
    CunsIdentifier() {};
    const CunsIdentifier  & operator=(const CunsIdentifier & m) {
      ident = m.ident;
      obj   = m.obj;
      return *this;
    };
    int ident;
    //CunsIn * obj;
    void * obj;
    static int getUnsvIndex(const int,CunsIdentifierVector *);
  };
 
}

#endif //UNSIDENTIFIER_H
