// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2013                                       
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
class  ParticlesObject;
typedef std::vector <ParticlesObject> ParticlesObjectVector; 

class UserSelection{
// this class allows to the user to check if the components/range she/he has
// selected match with the components returned by the snapshot's plugin.    
//                                                                          
public:
  UserSelection();
  //const UserSelection& operator=(const UserSelection& m);
  ~UserSelection();
  bool setSelection(const std::string,const ComponentRangeVector *, bool nodata=false);
  const t_indexes_tab * getIndexesTab() const { return indx; }
  int   getNSel()       const { return nsel   ; }
  static std::string parseString(std::string&);
  ComponentRangeVector * getCrvFromSelection() { return &crvsel;}
  int compBits() { return comp_bits; }
  std::vector <int> selectOrder() {
    return select_order;
  }
  void setCrv(ComponentRangeVector _crv) {
    crvsel = _crv;
  }

private:
  ParticlesObjectVector pov;
  static int comparePos(const void * a, const void * b) {
    t_indexes_tab * aa = (t_indexes_tab *) a;
    t_indexes_tab * bb = (t_indexes_tab *) b;
    return (aa->p - bb->p);
  }
  std::vector <int> select_order;
  std::string select;           // input range (console | GUI)
  bool parse();                 // parse selection
  int parseComponent(const std::string);
  bool checkComponent(const std::string);
  int isComponent(const std::string);
  int isRange(const std::string);
  bool nodata; // no data before setSelection
  std::string out_range;            // range selected             
  int nbody;                        // #bodies max                
  int nsel;                         // #bodies selected           
  int pos;                          // current component position 
  int comp_bits;                    // components bits requested
  const ComponentRangeVector * crv; // vector of component range  
  ComponentRangeVector crvsel;      // crv for selected components
  void fillIndexes(const std::string, const int, const int, const int, int);
  void crvResize(ComponentRangeVector&);
  int crvPermut(ComponentRange&, const int, const int, int&);
  t_indexes_tab * indx;
  int min,max;
  void findMinMax(const int, const int);
};
//
// class ParticlesObject
//
// this class is used to re-ordering the object according to
// user selection
class ParticlesObject{
  public:
  ParticlesObject() {
    npart=0;    
    step=first=last=-1;
    pos=-1;
  }

  const ParticlesObject& operator=(const ParticlesObject&m) {
    npart = m.npart;
    first = m.first;
    last  = m.last;
    pos   = m.pos;
    return *this;
  }

  int npart;       // #particles in the object
  int first;       // index of the first particle
  int last;        // index of the last particle
  int step;        // incremental step between particles.
  int pos;         // position of the object
};

}
#endif
