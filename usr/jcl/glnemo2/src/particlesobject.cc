// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "glwindow.h"
#include "particlesobject.h"
#include <iostream>
#include <assert.h>

namespace glnemo {
int ParticlesObject::nobj=0;
long long int ParticlesObject::cpt=0;
// ============================================================================
// constructor                                                                 
//ParticlesObject::ParticlesObject()
//{
//  index_tab    =  NULL;
//}
// ============================================================================
// constructor                                                                 
ParticlesObject::ParticlesObject(const ObjFrom _of, const std::string _name)
{
  init(_of,_name);
}
// ============================================================================
// constructor                                                                 
ParticlesObject::ParticlesObject( const int _npart, const int _first,
                                  const int _last, const int _step,
                                  const ObjFrom _of, const std::string _name)
{
  init(_of,_name);
  npart = _npart;
  first = _first;
  last  = _last;
  step  = _step;
}

// ============================================================================
// copy data object                                                             
void ParticlesObject::copyDataObject(const ParticlesObject&m, const bool garbage)
{
  obj_from     = m.obj_from;
  obj_name     = m.obj_name;
  npart        = m.npart;
  first        = m.first;
  last         = m.last;
  step         = m.step;
  
  copyProperties(m);
  
  // loop on orbits_max
  OrbitsVector oo = m.ov;
  for (OrbitsVector::iterator oit =oo.begin();oit!=oo.end() ; oit++) {  
    OrbitsList ol;
    // loop on orbit_history
    for (OrbitsList::iterator oil=(*oit).begin(); oil != (*oit).end(); oil++){
      Orbits * ob = new Orbits(*oil); // get each orbits
      ol.push_back(*ob);              // add in new list
      delete ob;
    }
    ov.push_back(ol);  // insert new list in new object
  }
  if (m.index_tab) {
    if (garbage && index_tab)
      delete [] index_tab;
    cpt +=npart;
    index_tab    = new int[npart];
    memcpy(index_tab,m.index_tab,sizeof(int)*npart);
    //for (int i=0;i<npart;i++) index_tab[i] = m.index_tab[i];
  }
  else index_tab=NULL;
  
}

// ============================================================================
// copy onstructor                                                             
ParticlesObject::ParticlesObject(const ParticlesObject&m)
{
#if 1
  copyDataObject(m);
#else
  obj_from     = m.obj_from;
  obj_name     = m.obj_name;
  npart        = m.npart;
  first        = m.first;
  last         = m.last;
  step         = m.step;
  visible      = m.visible;
  part         = m.part;
  part_size    = m.part_size;
  part_alpha   = m.part_alpha;
  gaz          = m.gaz;
  gaz_size     = m.gaz_size;
  gaz_alpha    = m.gaz_alpha;
  gaz_size_max = m.gaz_size_max;
  gaz_rotate   = m.gaz_rotate;
  gaz_glsl     = m.gaz_glsl;
  texture_index=m.texture_index;
  vel          = m.vel;
  vel_size     = m.vel_size;
  vel_alpha    = m.vel_alpha;
  vel_factor   = m.vel_factor;
  vel_size_max = m.vel_size_max;
  color        = m.color;
  pos          = m.pos;
  orbits       = m.orbits;
  o_record     = m.o_record;
  orbits_max   = m.orbits_max;
  orbits_history=m.orbits_history;
  orbits_animate=m.orbits_animate;
  min_phys   = m.min_phys;
  max_phys   = m.max_phys;
  min_percen_phys = m.min_percen_phys;
  max_percen_phys = m.max_percen_phys;
  has_physic      = m.has_physic;
  OrbitsVector oo = m.ov;
  // loop on orbits_max
  for (OrbitsVector::iterator oit =oo.begin();oit!=oo.end() ; oit++) {
    OrbitsList ol;
    // loop on orbit_history
    for (OrbitsList::iterator oil=(*oit).begin(); oil != (*oit).end(); oil++){
      Orbits * ob = new Orbits(*oil); // get each orbits
      ol.push_back(*ob);              // add in new list
      delete ob;
    }
    ov.push_back(ol);  // insert new list in new object
  }
  if (m.index_tab) {
    cpt +=npart;
    index_tab    = new int[npart];
    memcpy(index_tab,m.index_tab,sizeof(int)*npart);
    //for (int i=0;i<npart;i++) index_tab[i] = m.index_tab[i];
  }
  else index_tab=NULL;
#endif
}
// ============================================================================
// copy constructor                                                            
const ParticlesObject& ParticlesObject::operator=(const ParticlesObject&m)
{
#if 1
  copyDataObject(m,true);
#else
  obj_from     = m.obj_from;
  obj_name     = m.obj_name;
  npart        = m.npart;
  first        = m.first;
  last         = m.last;
  step         = m.step;
  visible      = m.visible;
  part         = m.part;
  part_size    = m.part_size;
  part_alpha   = m.part_alpha;
  gaz          = m.gaz;
  gaz_size     = m.gaz_size;
  gaz_alpha    = m.gaz_alpha;
  gaz_size_max = m.gaz_size_max;
  gaz_rotate   = m.gaz_rotate;
  gaz_glsl     = m.gaz_glsl;
  texture_index= m.texture_index;
  vel          = m.vel;
  vel_size     = m.vel_size;
  vel_alpha    = m.vel_alpha;
  vel_factor   = m.vel_factor;
  vel_size_max = m.vel_size_max;  
  color        = m.color;
  pos          = m.pos;
  orbits       = m.orbits;
  o_record     = m.o_record;
  orbits_max   = m.orbits_max;
  orbits_history=m.orbits_history;
  orbits_animate=m.orbits_animate;
  min_phys   = m.min_phys;
  max_phys   = m.max_phys;
  min_percen_phys = m.min_percen_phys;
  max_percen_phys = m.max_percen_phys;
  has_physic      = m.has_physic;
  OrbitsVector oo = m.ov;
  // loop on orbits_max
  for (OrbitsVector::iterator oit =oo.begin();oit!=oo.end() ; oit++) {
    OrbitsList ol;
    // loop on orbit_history
    for (OrbitsList::iterator oil=(*oit).begin(); oil != (*oit).end(); oil++){
      Orbits * ob = new Orbits(*oil); // get each orbits
      ol.push_back(*ob);              // add in new list
      delete ob;
    }
    ov.push_back(ol);  // insert new list in new object
  }
  if (m.index_tab) {
    if (index_tab)
      delete [] index_tab;
    cpt += npart;
    index_tab    = new int[npart];
    memcpy(index_tab,m.index_tab,sizeof(int)*npart);
    //for (int i=0;i<npart;i++) index_tab[i] = m.index_tab[i];
  }
  else index_tab=NULL;
#endif
  return *this;
}
// ============================================================================
// copyProperties                                                             
void ParticlesObject::copyProperties(const ParticlesObject&m)
{
  visible      = m.visible;
  part         = m.part;
  part_size    = m.part_size;
  part_alpha   = m.part_alpha;
  gaz          = m.gaz;
  gaz_size     = m.gaz_size;
  gaz_alpha    = m.gaz_alpha;
  gaz_size_max = m.gaz_size_max;
  gaz_rotate   = m.gaz_rotate;
  gaz_glsl     = m.gaz_glsl;
  texture_index=m.texture_index;
  vel          = m.vel;
  vel_size     = m.vel_size;
  vel_alpha    = m.vel_alpha;
  vel_factor   = m.vel_factor;
  vel_size_max = m.vel_size_max;
  color        = m.color;
  pos          = m.pos; // 02/sep/2011 added
  orbits       = m.orbits;
  o_record     = m.o_record;
  orbits_max   = m.orbits_max;
  orbits_history=m.orbits_history;
  orbits_animate=m.orbits_animate;
  min_phys   = m.min_phys;
  max_phys   = m.max_phys;
  min_percen_phys = m.min_percen_phys;
  max_percen_phys = m.max_percen_phys;
  has_physic      = m.has_physic;
  //ol           = m.ol;
  //pos          = m.pos;
}
// ============================================================================
// copyVVkeepProperties                                                        
void ParticlesObject::copyVVkeepProperties(ParticlesObjectVector& src,ParticlesObjectVector& dest, const int nsel)
{
  for (unsigned int i=0; i<src.size(); i++) {
    if (i<dest.size()) { // object already exist
      ParticlesObject::nobj--;
      ParticlesObject * podest = new ParticlesObject();
      *podest = dest[i]; // copy dest po
      dest[i] = src[i];  // src to dest
      dest[i].copyProperties(*podest); // get back properties
      delete (podest);
    }
    else {               // new object
      dest.push_back(src[i]);
    }
  }
  // check if all the object fit in memory
  for (unsigned int i=src.size(); i<dest.size(); i++) {
    if (dest[i].last>nsel) {      // object out of range
      dest.erase(dest.begin()+i); // remove object      
    }
  }
}
// ============================================================================
// copyVVProperties
void ParticlesObject::backupVVProperties(ParticlesObjectVector& src,ParticlesObjectVector& dest, const int nobj)
{
  for (int i=0; i<(int)src.size()&&i<nobj; i++) {
    if (i<(int)dest.size()) { // object already exist
      dest[i].copyProperties(src[i]); // copy properties only
    }
    else {               // new object
      ParticlesObject * podest = new ParticlesObject();
      dest.push_back(*podest);
      delete podest;
      dest[i].copyProperties(src[i]); // copy properties only
    }
  }
}
// ============================================================================
// clearOrbitsVectorPOV                                                        
void ParticlesObject::clearOrbitsVectorPOV(ParticlesObjectVector& pov)
{
  for (ParticlesObjectVector::iterator pvit=pov.begin(); pvit!=pov.end(); pvit++) {
    //std::cerr << "Object size (before) =" << (*pvit).ov.size() << "\n";
    for (OrbitsVector::iterator ovit=(*pvit).ov.begin(); ovit!=(*pvit).ov.end();ovit++) {
      (*ovit).clear();
    }
  }
}
// ============================================================================
// initOrbitsVectorPOV                                                        
void ParticlesObject::initOrbitsVectorPOV(ParticlesObjectVector& pov)
{
  for (ParticlesObjectVector::iterator pvit=pov.begin(); pvit!=pov.end(); pvit++) {
    //std::cerr << "Object size (before) =" << (*pvit).ov.size() << "\n";
    //std::cerr << "pvit->npart = " << pvit->npart << " orbit max ="<<pvit->getOrbitsMax()<<"\n";
    pvit->setOrbitsMax(std::min(pvit->npart,pvit->getOrbitsMax()));
  }
}
// ============================================================================
// destructor                                                                  
ParticlesObject::~ParticlesObject()
{
  if (index_tab) {
    cpt -= npart;
    delete [] index_tab;
  }
  int nn=0;
  for (OrbitsVector::iterator oit =ov.begin();oit!=ov.end() ; oit++) {    
    (*oit).clear();  
    nn++;
  }
  ov.clear();
//   if (ol && ol->size())
//     ol->clear(); // Orbits List

}
// ============================================================================
// init                                                                        
void ParticlesObject::init(const ObjFrom _of, const std::string _name)
{
  obj_from = _of;
  obj_name = _name;
  pos = nobj;
  nobj++;
  npart        =  0;
  first        = -1;
  last         = -1;
  step         =  1;
  visible      =  true;
  part         =  false;
  part_size    =  1;
  part_alpha   =  255;
  gaz          =  true;
  gaz_size     =  1.;
  if (!GLWindow::GLSL_support) {
    gaz_size   =  0.1;
    gaz        =  false;
    part       =  true;
  }
  gaz_alpha    =  255;
  gaz_size_max =  1.0;
  gaz_rotate   = true;
  gaz_glsl     = true;
  texture_index= 0;
  vel          =  false;
  vel_size_max =  4.;
  vel_size     =  1;
  vel_alpha    =  255;
  vel_factor   =  1.0;
  index_tab    =  NULL;
  freed        = false;
  orbits       = false;
  o_record     = false;
  orbits_max   = 20;
  orbits_history=50;
  orbits_animate=false;
  min_phys   = -1.;//0.;
  max_phys   = -1.;//10000000000.;
  min_percen_phys=0;
  max_percen_phys=99;
  has_physic = false;
  setColor();

}
// ============================================================================
// buildIndexList                                                              
void ParticlesObject:: buildIndexList(std::vector<int> & indexes)
{
  npart = indexes.size();
  if (index_tab) {
    delete [] index_tab;
    cpt -= npart;
  }
  first = 0;
  last  = npart-1;
  step  = 1;
  index_tab    = new int[npart];
  cpt += npart;
  for (int i=0; i<npart; i++) {
    index_tab[i] = indexes[i];
  }
}
// ============================================================================
// buildIndexList                                                              
void ParticlesObject::buildIndexList( const int _npart, const int _first,
                                      const int _last, const int _step)
{
  if (index_tab) {
    delete [] index_tab;
    cpt -= npart;
  }
  npart = _npart;
  first = _first;
  last  = _last;
  step  = _step;
  index_tab    = new int[npart];
  cpt += npart;
  for (int i=0; i<npart; i++) {
    index_tab[i] = first+i;
  }
}
// ============================================================================
// buildIndexList                                                              
void ParticlesObject::buildIndexList()
{
  buildIndexList(npart,first,last,step);
}
// ============================================================================
// resizeRange                                                                 
int ParticlesObject::resizeRange(const int min,const int max,int &next_first)
{
  int ret;
  if (last > max) ret=last;  // new max
  else ret=max;              // keep max

  if (first <= max ) {       // first <= max
    first -= min;            // left move
    last=first+npart-1;      // resize last
  }
  else {                     // first > max
    first = next_first;      // move to the right of latest object
    last=first+npart-1;      // resize last
  }
  next_first=last+1;
  return ret;
}
// ============================================================================
// setColor                                                                    
void ParticlesObject::setColor()
{
  int modulo_col[5][3] = {
                            { 255, 75, 39 },
                            { 255,255, 255 },
                            { 114,214, 32 },
                            {  58, 61,214 },
                            { 214, 47,197 }
  };
  //{ 214,214, 52 },
  //std::cerr << "NOBJ = "<<nobj<<"\n";
  color = QColor(modulo_col[(nobj-1)%5][0],
                 modulo_col[(nobj-1)%5][1],
                 modulo_col[(nobj-1)%5][2]);
}
// ============================================================================
// addOrbits                                                                   
void ParticlesObject::addOrbits(const ParticlesData * p_data)
{
  if (isOrbitsRecording()) {
    if ((int)ov.size() != orbits_max) {
      std::cerr << "resize orbits vector\n";
      ov.resize(orbits_max);
    }
    assert((int)ov.size()==orbits_max);
    int no=0;
    for (OrbitsVector::iterator oit=ov.begin(); oit!=ov.end(); oit++) {
      int size_current_list = (*oit).size();  // size of the current list
      Orbits * moreorbits = new Orbits(p_data,index_tab[no]); // generate one history
      if (size_current_list < orbits_history) { // max history *not* reached
        //std::cerr << "add new orbits\n";
        (*oit).push_back(*moreorbits);          // add new orbits record at the end
        delete moreorbits;
      }
      else {                                    // max history reached
        //std::cerr << "shift orbits\n";
        (*oit).pop_front();                     // remove first element
        (*oit).push_back(*moreorbits);          // add new orbits record at the end
        delete moreorbits;
      }
      no++; // next orbit
    }
  }
}
}
