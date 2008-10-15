// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#ifndef GLNEMOCOMPONENTRANGE_H
#define GLNEMOCOMPONENTRANGE_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
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
  static int print(const ComponentRangeVector * crv, std::string select);
   std::string range,   // "0:99999"
               type;    // "all", "halo", "disk", "bulge", "gas", "other"....
   int first,last,n;
   
  private:
    int computeN();
    void buildRange();
};

}

#endif
