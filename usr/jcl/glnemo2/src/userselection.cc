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
#include <iostream>
#include <sstream>
#include <QRegExp>
#include <QString>
#include <QTime>
#include "userselection.h"
#include <algorithm>
#include "assert.h"
namespace glnemo {

// ============================================================================
// constructor                                                                 
UserSelection::UserSelection()
{
  nbody   = 0;
  nsel    = 0;
  indexes = NULL;
  povcpy2 = new ParticlesObjectVector();
}
// ============================================================================
// destructor                                                                  
UserSelection::~UserSelection()
{
  if (indexes) {
    delete [] indexes;
  }
  povcpy2->clear();
}
// ============================================================================
// setSelection:                                                               
// according to User's selection (string select) and the differents components 
// (ranges, components) inside the snapshot, the POV is filled up.             
bool UserSelection::setSelection(std::string _sel,
                                 const ComponentRangeVector * _crv,
                                 ParticlesObjectVector * _pov)
{
  select =_sel; // cop selection              
  crv    = _crv;   // link component range vector
  pov    = _pov;   // link Particle Object Vector
  if (select=="all" && crv->size()>1) { // replace "all" by list of components
    select=(*crv)[1].type;
    for (unsigned int i=2;i<crv->size();i++) {
      select += ","+(*crv)[i].type;
    }
    std::cerr << "New selected component=" << select << "\n";
  }
  assert(crv);                   // must not be NULL
  assert((*crv)[0].type=="all"); // first entry must be "all"    
  nbody = (*crv)[0].n;           // #bodies max in the snapshot  
  if (indexes) 
     delete [] indexes;
  if (pov->size()) { // copy previous objects
     povcpy2->clear();
     ParticlesObject::nobj=0;
     *povcpy2 = *pov;
  }
  pov->clear();                   // reset particles object vector
  ParticlesObject::nobj=0;
  indexes = new int[nbody];
  for (int i=0;i<nbody;i++) indexes[i]=-1; // reset indexes
  nsel = 0;
  min = max = -1;

  bool status=parse();
  if (status || 1 ) { // we force here
    QTime tbench;
    tbench.restart();
    // ascending sort according to the 'first' element
    int nobj=ParticlesObject::nobj;
    std::sort(pov->begin(),pov->end(),ParticlesObject::compareFirst);
    ParticlesObject::nobj=nobj;
    
    int new_max=max;
    int next_first=0;
    for (unsigned int i=0; i<pov->size(); i++) {
      new_max=(*pov)[i].resizeRange(min,new_max,next_first); // resize range
      (*pov)[i].buildIndexList();                            // rebuild index list
    }
    // ascending sort according to the 'pos' element
    nobj=ParticlesObject::nobj;
    std::sort(pov->begin(),pov->end(),ParticlesObject::comparePos);
    ParticlesObject::nobj=nobj;

    for (unsigned int i=0; i<pov->size(); i++) {
      if (i < povcpy2->size()) {  // old object exist
        (*pov)[i].copyProperties((*povcpy2)[i]);
      }
    }
    std::cerr << "UserSelection::setSelection, Time elapsed for sorting="<< 
        tbench.elapsed()/1000 << "\n";
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
  while ((current_s=parseString(next_s)) != "" && (status||1)) {  // || 1 force parsing
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
// return 0 is the component is type of range
int UserSelection::isRange(const std::string comp)
{
  int status;
  // Regular expression => first:last:step
  QRegExp rx("^(\\d{1,})((:)(\\d{1,})){,1}((:)(\\d{1,})){,1}$");
  int match=rx.indexIn(QString(comp.c_str()));
  if (match == -1) { // not match
    status=1;        // misformated
  }
  else {
    int first,last,step=1;
    // get first
    std::istringstream iss((rx.cap(1)).toStdString());
    iss >> first;
    // get last
    if (rx.captureCount()>4) {
      std::istringstream  iss((rx.cap(4)).toStdString());
      iss >> last;
      // get step
      if (rx.captureCount()>=7) {
        std::istringstream  iss((rx.cap(7)).toStdString());
        iss >> step;
      }
    }
    else {
      last=first;
    }
    if (last>=first) {
      int npart=last-first+1; // #part 
      if (npart <= nbody) {   // valid range
        status=0;
        fillIndexes(first,last,step); // fill indexes array
        storeObject(first,last,step); // store pov
      }
      else {
        status=3; // range too big
      }
    }
    else {
      status=2;   // range misformated
    }
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
  QRegExp rx("^(all|halo|disk|disc|bulge|stars|gas|gaz|bndry|other(\\d{,}))$");
  int match=rx.indexIn(QString(comp.c_str()));
  if (match == -1) { // not match
    status=1;        // misformated
  }
  else {
    int first,last,step=1;
    // get component's type
    std::string type=(rx.cap(1)).toStdString();
    int icrv=ComponentRange::getIndexMatchType(crv,type);
    if (icrv != -1 ) {
      assert(icrv<(int)crv->size());
      first=(*crv)[icrv].first;
      last =(*crv)[icrv].last;
      assert(last>=first);
      int npart=last-first+1; // #part
      assert(npart<=nbody);
      status=0;
      fillIndexes(first,last,step); // fill indexes array
      storeObject(first,last,step); // store pov
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
void UserSelection::fillIndexes(const int first, const int last, const int step)
{
  int npart=(last-first+1)/step;
  assert(npart<=nbody);
  for (int i=first; i<=last; i+=step) {
    if (indexes[i]==-1) nsel++; // one more particles
    indexes[i] = i;             // set new particles 
    assert(nsel<=nbody);
  }
}
// ============================================================================
// storeObject                                                                 
// create a new Particle Object                                                
void UserSelection::storeObject(const int first, const int last, const int step)
{
  int npart=(last-first+1)/step;
  assert(npart<=nbody);
  int nlast=first+npart-1;        // resize last in case step > 1
  // step =1  by default now
  findMinMax(first,nlast); // find min and max for later processing
  ParticlesObject * po = new ParticlesObject(npart,first,nlast); // new object
  pov->push_back(*po);                           // store in vector
  delete po;
}
// ============================================================================
// findMinMax                                                                  
//                                                                             
void UserSelection::findMinMax(const int first, const int last)
{
  if (min == -1 ) min=first;
  if (max == -1 ) max=last;
  if (first<min) {
    min = first;
    max = last;
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

