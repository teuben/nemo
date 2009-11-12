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
#include <vector>
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
  pos = 0; // current component selection
  select =_sel; // cop selection              
  crv    = _crv;   // link component range vector
  assert(crv);                   // must not be NULL
  assert((*crv)[0].type=="all"); // first entry must be "all"    
  nbody = (*crv)[0].n;           // #bodies max in the snapshot  
  std::cerr << "nbody = "<<nbody<<"\n";
  if (indx) 
    delete [] indx;
  indx = new t_indexes_tab[nbody];
  
  for (int i=0;i<nbody;i++) {
    indx[i].i=-1;    // reset indexes
    indx[i].p=10000; // set position to high number
  }
  nsel = 0;

  bool status=parse();
  if (status || 1 ) { // we force here
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
  while ((current_s=parseString(next_s)) != "") {// && (status)) {  // || 1 force parsing
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
  if ((fail=isRange(comp)))
    if((fail=isComponent(comp)))
      status=false;
  return status;
}
// ============================================================================
// isRange                                                                     
// return true is the component is a range of particles                        
int UserSelection::isRange(const std::string comp)
{
  int status;
  std::vector<int> store;
  int ppos=0;
  bool stop=false;
  int cpt=0;
  while (! stop) {
    size_t found = comp.find(':',ppos);
    if (found!=std::string::npos) {
      if (found > (size_t) (ppos)) {
	cpt++;
	std::string str=comp.substr(ppos,found-ppos);
	std::istringstream ss(str);
	int val;
	ss>>val;
	store.push_back(val);
      }
      ppos=found+1; //
    } else { // no more ":"
      if (cpt>0) {
	std::string str=comp.substr(ppos);
	std::istringstream ss(str);
	int val;
	ss>>val;
	store.push_back(val);
      }
      stop=true;
    }
  }
  for (std::vector<int>::iterator it=store.begin(); it!=store.end(); it++) {
    //std::cerr << " range i ="<< it-store.begin() << "\n";
    //std::cerr << " value   =" << *it << "\n";
  }
  if (cpt>0) { // one ":" has been found
    int step=1;
    int first=store[0];
    int last=first;
    if (store.size()>1) {
      last = store[1];
    }
    if (store.size()>2) {
      step = store[2];
    }
    assert(last>=first);
    int npart=last-first+1; // #part
    assert(npart<=nbody);
    fillIndexes(first,last,step,pos); // fill indexes array
    pos++;
    status=0;
  }
  else {
    status=1;
  }


  return status;
}
// ============================================================================
// isComponent                                                                 
// return true is the component is component type                              
int UserSelection::isComponent(const std::string comp)
{
  int status;
  // Regular expression => all|halo|disk ......
  const char *  rx[] = {"all","halo","disk","bulge","stars","gas","bndry","halo2",NULL};
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
    int offset;
    int icrv=ComponentRange::getIndexMatchType(crv,type,offset);
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

