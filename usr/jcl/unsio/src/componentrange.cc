// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2013                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#include "componentrange.h"
#include <sstream>
#include <iostream>
#include "assert.h"
namespace uns {

// ============================================================================
// constructor                                                                 
ComponentRange::ComponentRange()
{
  n=0;
  position=first=last=-1;
  range="";
  type="";
}
// ============================================================================
// copy constructor                                                            
ComponentRange::ComponentRange(const ComponentRange&m)
{
  n        = m.n;
  first    = m.first;
  last     = m.last;
  range    = m.range;
  type     = m.type;
  position = m.position;
}
// ============================================================================
ComponentRange::~ComponentRange()
{
}
// ============================================================================
void ComponentRange::setData(const int _f, const int _l, const std::string _t)
{
  first=_f;
  last = _l;
  setType(_t);
  computeN();
  buildRange();
}
// ============================================================================
int ComponentRange::computeN()
{
  n=last-first+1;
  return n;
}
// ============================================================================
void ComponentRange::buildRange()
{
  std::ostringstream stm1,stm2;
  stm1 << first;
  stm2 << last;
  range = stm1.str()+":"+stm2.str();
}
// ============================================================================
// getIndexMatchType                                                           
// return the index vector matching the type of object in the component range  
// vector. Return also the object offset in the index tab                      
int ComponentRange::getIndexMatchType(const ComponentRangeVector * crv, const std::string type,
                                      int &offset, bool crvuser)
{
  int status=-1;
  offset=0;
  assert(crv);
  for (unsigned int i=0; i<crv->size()&&status==-1; i++) {
    //std::cerr << "i="<<i<<" crv.type =" <<(*crv)[i].type<<" type="<<type<<"\n"; 
    if ((*crv)[i].type == type)  status=i;
    else {
      if (i>0 || crvuser) { // skip first which match "all" component
                            // or User's crv is processing "no all"
	offset+=(*crv)[i].n;
      }
    }
  }
  return status;
}
// ============================================================================
void ComponentRange::list(const ComponentRangeVector * crv)
{
  std::cerr << "ComponentRange::list size"<<crv->size()<<"\n";
  for (unsigned int i=0; i<crv->size(); i++) {
    std::cerr << "-----------------------------------------------------------\n";
    std::cerr << "Component #"<<i<<"\n";
    std::cerr << "type  :"<<(*crv)[i].type<<"\n";
    std::cerr << "range :"<<(*crv)[i].range<<"\n";
    std::cerr << "nbody :"<<(*crv)[i].n<<"\n";
  }
}
// ============================================================================
int ComponentRange::print(const ComponentRangeVector * crv, std::string select)
{
  for (unsigned int i=0; i<crv->size(); i++) {
    if ((*crv)[i].type == select) return i;
  }
  return -1;
}
} // namespace
// ============================================================================

