// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

#include <iostream>
#include <sstream>
#include "userselection.h"
#include <stdlib.h>
#include "assert.h"
namespace glnemo {

// ============================================================================
// constructor                                                                 
UserSelection::UserSelection()
{
  nbody   = 0;
  nsel    = 0;
  indx    = NULL;
  pos     = 0;
}
// ============================================================================
// destructor                                                                  
UserSelection::~UserSelection()
{
  if (indx) {
    delete [] indx;
  }
}
// ============================================================================
// setSelection:                                                               
// according to User's selection (string select) and the differents components 
// (ranges, components) inside the snapshot, the POV is filled up.             
bool UserSelection::setSelection(const std::string _sel,
                                 const ComponentRangeVector * _crv)
{
  select =_sel; // cop selection              
  crv    = _crv;   // link component range vector
  assert(crv);                   // must not be NULL
  assert((*crv)[0].type=="all"); // first entry must be "all"    
  nbody = (*crv)[0].n;           // #bodies max in the snapshot  
  if (indx) 
    delete [] indx;
  indx = new t_indexes_tab[nbody];
  for (int i=0;i<nbody;i++) {
    indx[i].i=-1;    // reset indexes
    indx[i].p=10000; // set position to high number
  }
  nsel = 0;

  bool status=parse();
  if (status  ) { // we force here
    qsort(indx,nbody,sizeof(t_indexes_tab),UserSelection::comparePos);
  }
  return status;
}
// ============================================================================
// parse:                                                                      
bool UserSelection::parse()
{
  bool status=true;
  std::string current_s,next_s;
  next_s = select;
  while ((current_s=parseString(next_s)) != "" && (status)) {  // || 1 force parsing
    status=checkComponent(current_s);
  }
  return status;
}
// ============================================================================
// checkComponent                                                              
// find out the component stored in the string                                 
bool UserSelection::checkComponent(const std::string comp)
{
  bool status=true;
  int fail;
  if((fail=isComponent(comp)))
    status=false;
  return status;
}
// ============================================================================
// isComponent                                                                 
// return true is the component is component type                             
int UserSelection::isComponent(const std::string comp)
{
  int status;
  // Regular expression => all|halo|disk ......
  const char *  rx[] = {"all","halo","disk","bulge","stars","gas","bndry",NULL};
  int i=0;
  int match=-1;
  while (rx[i] && match==-1) {
    if (rx[i]) {
      std::string tmp=rx[i];
      if (tmp == comp) match=i;
    }
    i++;
  }
  if (match == -1) { // not match
    status=1;        // misformated
  }
  else {
    int first,last,step=1;
    // get component's type
    std::string type=rx[match];
    int icrv=ComponentRange::getIndexMatchType(crv,type);
    if (icrv != -1 ) {
      assert(icrv<(int)crv->size());
      first=(*crv)[icrv].first;
      last =(*crv)[icrv].last;
      assert(last>=first);
      int npart=last-first+1; // #part
      assert(npart<=nbody);
      status=0;
      fillIndexes(first,last,step,pos); // fill indexes array
      pos++;
    }
    else {
      status=4;      // type does not exist
    }
  }
  return status;
}
// ============================================================================
// fillIndexes                                                                 
// fill indexes array with new range of particles                              
  void UserSelection::fillIndexes(const int first, const int last, const int step, int pos)
{
  int npart=(last-first+1)/step;
  assert(npart<=nbody);
  for (int i=first; i<=last; i+=step) {
    if (indx[i].i==-1) nsel++; // one more particles
    indx[i].i = i;             // set new particles 
    indx[i].p = pos;
    assert(nsel<=nbody);
  }

}
// ============================================================================
// parseString                                                                 
// return the string at the position after the next ',' otherwise ""           
std::string UserSelection::parseString(std::string & next_string)
{
  std::string return_string;
  std::string::size_type coma=next_string.find(",",0);  // try to find ","
  if (coma != std::string::npos) { // found ","
    return_string = next_string.substr(0,coma);
    next_string   = next_string.substr(coma+1,next_string.length());
  } else {                         // not found
    return_string = next_string;
    next_string = "";
  }
  return return_string;
}
} //namespace glnemo

