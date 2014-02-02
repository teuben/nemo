// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef GLNEMOUSERSELECTION_H
#define GLNEMOUSERSELECTION_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#include <string>
#include <vector>
#include "componentrange.h"
#include "particlesobject.h"

namespace glnemo {

class UserSelection{
// this class permit to the user to check if the components/range she/he has
// selected match with the components returned by the snapshot's plugin.    
//                                                                          
// if comp is empty, only range is taken in account                         
public:
    UserSelection();
    ~UserSelection();
    bool setSelection(const std::string,const ComponentRangeVector *,
                      ParticlesObjectVector *);
    const int * getIndexesTab() const { return indexes; };
    int         getNSel()       const { return nsel   ; };
private:
  
  std::string select;           // input range (console | GUI)
  bool parse();                 // parse selection
  std::string parseString(std::string&);
  int parseComponent(const std::string);
  bool checkComponent(const std::string);
  int isRange(const std::string);
  int isComponent(const std::string);
  std::string out_range;            // range selected             
  int * indexes;                    // particles' indexes selected
  int nbody;                        // #bodies max                
  int nsel;                         // #bodies selected           
  const ComponentRangeVector * crv; // vector of component range
  void fillIndexes(const int, const int, const int);
  void storeObject(const int, const int, const int);
  ParticlesObjectVector * pov, * povcpy2;
  int min,max;
  void findMinMax(const int, const int);
};
}
#endif
