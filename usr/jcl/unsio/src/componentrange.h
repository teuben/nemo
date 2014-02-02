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
#ifndef UNSCOMPONENTRANGE_H
#define UNSCOMPONENTRANGE_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#include <vector>
#include <string>

namespace uns {

class ComponentRange;

typedef std::vector <ComponentRange> ComponentRangeVector;

class ComponentRange{
public:
    ComponentRange();
    ComponentRange(const ComponentRange&);  // copy constructor
    ~ComponentRange();

  void setData(const int _f,const  int _l, const std::string _type="");
  void setType(std::string _type) { type = _type; }
  void setPosition(const int _pos) { position = _pos;}
  static int getIndexMatchType(const ComponentRangeVector * crv, const std::string type, 
                               int & offset,bool crvuser=false);
  static void list(const ComponentRangeVector * crv);
  static int print(const ComponentRangeVector * crv, std::string select);
  static bool compareFirst(const ComponentRange &a ,const ComponentRange &b) {
    return a.first < b.first;
  };
  static bool comparePos(const ComponentRange &a ,const ComponentRange &b) {
    return a.position < b.position;
  };
  
  std::string range,   // "0:99999"
               type;    // "all", "halo", "disk", "bulge", "gas", "other"....
   int first,last,n;
   int position;
  private:
    int computeN();
    void buildRange();
};

}

#endif
