// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
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
#ifndef GLNEMOCOMPONENTRANGE_H
#define GLNEMOCOMPONENTRANGE_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#include <vector>
#include <string>

namespace glnemo {

class ComponentRange;

typedef std::vector <ComponentRange> ComponentRangeVector;

class ComponentRange{
public:
    ComponentRange();
    ComponentRange(const ComponentRange&);  // copy constructor
    ~ComponentRange();

  void setData(const int _f,const  int _l, const std::string _type="");
  void setType(std::string _type) { type = _type; };
  static int getIndexMatchType(const ComponentRangeVector * crv, const std::string type);
  static void list(const ComponentRangeVector * crv);
   std::string range,   // "0:99999"
               type;    // "all", "halo", "disk", "bulge", "gas", "other"....
   int first,last,n;
   
  private:
    int computeN();
    void buildRange();
};

}

#endif
