// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#ifndef UNSUSERSELECTION_H
#define UNSUSERSELECTION_H
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#include <string>
#include <vector>
#include "componentrange.h"

namespace uns {
#define HaloBit     (1  <<  0)
#define DiskBit     (1  <<  1)
#define BulgeBit    (1  <<  3)
#define StarsBit    (1  <<  4)
#define GasBit      (1  <<  5)
  struct indexes_tab { int i,p; };
  typedef struct indexes_tab t_indexes_tab;

class UserSelection{
// this class allows to the user to check if the components/range she/he has
// selected match with the components returned by the snapshot's plugin.    
//                                                                          
public:
  UserSelection();
  //const UserSelection& operator=(const UserSelection& m);
  ~UserSelection();
  bool setSelection(const std::string,const ComponentRangeVector *);
  const t_indexes_tab * getIndexesTab() const { return indx; }
  int   getNSel()       const { return nsel   ; }
  static std::string parseString(std::string&);
  ComponentRangeVector * getCrvFromSelection() { return &crvsel;}
  
private:
  static int comparePos(const void * a, const void * b) {
    t_indexes_tab * aa = (t_indexes_tab *) a;
    t_indexes_tab * bb = (t_indexes_tab *) b;
    return (aa->p - bb->p);
  }
  std::string select;           // input range (console | GUI)
  bool parse();                 // parse selection
  int parseComponent(const std::string);
  bool checkComponent(const std::string);
  int isComponent(const std::string);
  int isRange(const std::string);
  std::string out_range;            // range selected             
  int nbody;                        // #bodies max                
  int nsel;                         // #bodies selected           
  int pos;                          // current component position 
  const ComponentRangeVector * crv; // vector of component range  
  ComponentRangeVector crvsel;      // crv for selected components
  void fillIndexes(const std::string, const int, const int, const int, int);
  void crvResize(ComponentRangeVector&);
  int crvPermut(ComponentRange&, const int, const int, int&);
  t_indexes_tab * indx;
  int min,max;
  void findMinMax(const int, const int);
};
}
#endif
