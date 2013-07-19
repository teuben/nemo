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

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#include "snapshotramses.h"
#include "camr.h"
#include "cpart.h"
#include <limits>
#include "uns.h"

namespace uns {
// ============================================================================
// constructor : CSnapshotRamsesIn
CSnapshotRamsesIn::CSnapshotRamsesIn(const std::string _name,
                                     const std::string _comp,
                                     const std::string _time,
                                     const bool verb):
  CSnapshotInterfaceIn(_name, _comp, _time, verb)
{
  first_loc=true;
  particles = new CParticles();
  valid=false;
  part = new ramses::CPart(filename,verbose);
  amr  = new ramses::CAmr(filename,verbose);
  if (part->isValid() || amr->isValid()) {
    valid=true;
    interface_type = "Ramses";
    file_structure = "component";
    interface_index= 2;

    // >> !! NEEDS TO BE FIXED
    // all
    uns::ComponentRange cr;
    cr.setData(0,0);
    cr.setType("all");
    crv.clear();
    crv.push_back(cr);
    // << !! NEEDS TO BE FIXED

  }
}
// ============================================================================
// desstructor : ~CSnapshotRamsesIn
CSnapshotRamsesIn::~CSnapshotRamsesIn()
{
  delete amr;
  delete part;
  delete particles;
}
// ============================================================================
// PURE virtual FUNCTION delared in CSnapshotInterfaceIn
// ============================================================================

// ============================================================================
// getSnapshotRange
ComponentRangeVector * CSnapshotRamsesIn::getSnapshotRange()
{

  if (valid && crv.size()) {
    //crv = getCRV();
    //ComponentRange::list(&crv);
    if (first) {
      first       = false;
      crv_first   = crv;
      //nbody_first = getNtotal();
      //std::cerr << "CSnapshotGadgetIn::getSnapshotRange() = " << nbody_first << "\n";
      //time_first  = getTime();
    }
  }
  return &crv;
}
// ============================================================================
// nextFrame
int CSnapshotRamsesIn::nextFrame(uns::UserSelection &user_select)
{
  int status=0;
  assert(valid==true);
  if (first_loc) {
    first_loc = false;
    if (1 /*checkRangeTime(getTime())*/) {
      // check component bits selected
      user_select.setSelection(getSelectPart(),&crv,true);
      unsigned int comp_bits=user_select.compBits();

      // set boundaries
      float x[8];

      x[0]=x[2]=x[4]=std::numeric_limits<float>::min();
      x[1]=x[3]=x[5]=std::numeric_limits<float>::max();
      x[6]=(float )0; // level min
      x[7]=0; // nlevelmax

      if ((comp_bits&HALO_BIT || comp_bits&STARS_BIT) && part->isValid()) {
        part->setBoundary(x);
        part->loadData(particles,req_bits,comp_bits);
      }
      if (comp_bits&GAS_BIT && amr->isValid()) {
        std::cerr << "in gas\n";
        amr->setBoundary(x);
        amr->loadData(particles,req_bits);
      }
      std::cerr << "ntot   = "<< particles->ntot <<"\n";
      std::cerr << "ngas   = "<< particles->ngas <<"\n";
      std::cerr << "ndm    = "<< particles->ndm <<"\n";
      std::cerr << "nstars = "<< particles->nstars <<"\n";
      std::cerr << "Box len=" << amr->getHeader()->boxlen << "\n";
      std::cerr << "Start reordering...\n";
      reorderParticles(user_select);
      std::cerr << "Stop reordering...\n";
      status = 1;
    }
  }
  return status;
}
// ============================================================================
// close operation
int CSnapshotRamsesIn::close()
{
  return 1;
}
// ============================================================================
// reorderParticles
// here we re order particles according to user's component selection
// select_order is a vector containing indexes of components sorted
// according to user request
int CSnapshotRamsesIn::reorderParticles(uns::UserSelection &user_select)
{
  std::cerr <<"Nbody particles loaded="<<particles->ntot<<"\n";
  std::vector <int> offset_comp(6,-1); // init 6 offsets with value =-1
  std::vector <int> npart_comp(6,0);   // initialise #part per component to ZERO
  char * comp[] = { (char*) "gas",(char*) "halo",(char*) "disk",(char*) "bulge",(char*) "stars",(char*) "bndry",(char*) "all" };
  // get ordering
  std::vector <int> select_order = user_select.selectOrder();

  // component range
  uns::ComponentRange cr;
  crv.clear();

  // set #particles per component
  npart_comp[0] = particles->ngas;
  npart_comp[1] = particles->ndm;
  npart_comp[4] = particles->nstars;

  assert(select_order.size()<7);
  bool first=true;
  // according to user selection
  // we reformat offset of each components
  for (unsigned int i=0;i<select_order.size(); i++) {
    assert(select_order[i]<7);
    if (first) { // first in the order
      if (npart_comp[select_order[i]] > 0 ) { // part exist for the component
        offset_comp[select_order[i]]=0; // offset 0
        first=false;

      }
    } else { // remaining sorted component
      // current offset of the component
      assert(offset_comp[select_order[i-1]]!=-1);
      offset_comp[select_order[i]]=offset_comp[select_order[i-1]]+ // previous offset
                                   npart_comp[select_order[i-1]];  // size of previous component
    }
    if (!first && npart_comp[select_order[i]]>0) { // if first && #npart > 0
      cr.setData(offset_comp[select_order[i]],     // first part
                 offset_comp[select_order[i]]+     // last part
                 npart_comp[select_order[i]]-1,
                 comp[select_order[i]]);           // component name
      crv.push_back(cr);
    }
  }

  if (select_part=="all") { // if "all" selected
    uns::ComponentRangeVector crvall;
    uns::ComponentRange cr;
    // all
    cr.setData(0,particles->ntot-1);
    cr.setType("all");
    crvall.push_back(cr);
    user_select.setCrv(crvall);  // we force CRV with "all" only because select_part="all"
  } else {                  // if not "all" selected
    // set the new real user select
    user_select.setCrv(crv);  // we rebuild CRV with the user selection
  }
  if (verbose)
    uns::ComponentRange::list(&crv);

  // all
  cr.setData(0,particles->ntot-1);
  cr.setType("all");
  uns::ComponentRangeVector::iterator it;
  it = crv.begin();
  // insert "all" at the beginning of componentRangeVector
  crv.insert(it,cr);

  for (unsigned int i=0;i<offset_comp.size(); i++) {
    std::cerr << "i="<<i<<"["<<comp[i]<<"] offset="<<offset_comp[i]<<"\n";
  }
  if (verbose)
    uns::ComponentRange::list(&crv);

  // particles reordering
  if (particles->ntot) { // exist particles to reorder
    std::vector <float> pos,vel,mass,metal;
    std::vector <int> id;
    if (particles->pos.size()>0)
      pos.resize(particles->pos.size());      // resize new pos   vector
    if (particles->vel.size()>0)
      vel.resize(particles->vel.size());      // resize new pos   vector
    if (particles->mass.size()>0)
      mass.resize(particles->mass.size());    // resize new mass  vector
    if (particles->metal.size()>0)
      metal.resize(particles->metal.size());  // resize new metal vector
    if (particles->id.size()>0)
      id.resize(particles->id.size());     // resize new id vector

    for (int i=0; i<particles->ntot; i++) {

      bool found=false;

      int icomp=particles->indexes[i]; // integer component
      assert(icomp==0 ||icomp==1 || icomp==4); // gas || halo || stars only
      int istart=offset_comp[icomp]; // index start in the new pos array

      // positions
      if (particles->pos.size()>0) {
        assert((istart*3)+2<(int)particles->pos.size());
        found=true;
        pos[istart*3+0] = particles->pos[i*3+0]; // x
        pos[istart*3+1] = particles->pos[i*3+1]; // y
        pos[istart*3+2] = particles->pos[i*3+2]; // z
      }

      // velocities
      if (particles->vel.size()>0) {
        assert((istart*3)+2<(int)particles->vel.size());
        found=true;
        vel[istart*3+0] = particles->vel[i*3+0]; // x
        vel[istart*3+1] = particles->vel[i*3+1]; // y
        vel[istart*3+2] = particles->vel[i*3+2]; // z
      }

      // masses
      if (particles->mass.size()>0) {
        assert(istart<(int)particles->mass.size());
        found=true;
        mass[istart]    = particles->mass[i];
      }

      // id
      if (particles->id.size()>0) {
        assert(istart<(int)particles->id.size());
        found=true;
        id[istart]    = particles->id[i];
      }
      // metal
      if (particles->metal.size()>0) { // && (icomp==0 || icomp==4)) { // metal for gas or stars
        if (!(istart<(int)particles->metal.size())) {
          std::cerr << " istart ="<<istart<< " metal.size ="<< particles->metal.size() << "\n";
        }
        assert(istart<(int)particles->metal.size());
        found=true;
        metal[istart]    = particles->metal[i];
      }

      if (found) { // found particles
        offset_comp[icomp]++; //  update offset
      }
    }
    // copy back arrays
    particles->pos   = pos;
    particles->vel   = vel;
    particles->mass  = mass;
    particles->metal = metal;
    particles->id    = id;
    //std::cerr << "metal.size() ="<<particles->metal.size() <<"\n";
  }
  return 1;
}
// ============================================================================
// getData
// return requested array according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string comp, std::string name, int *n,float **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;

  int nbody,first,last;
  bool status=getRangeSelect(comp.c_str(),&nbody,&first,&last,false); // find components ranges
  if (!status && comp=="all") { // retreive all particles selected by the user
    status=1;
    first=0;
    nbody=particles->ntot;
  }
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Nbody :
    if (status) {
      *data = NULL;
      *n = nbody;
    } else {
      ok = false;
    }
    break;
  case uns::Nsel   :
    if (status) {
      *n    = nbody;
    } else {
      ok=false;
    }
  case uns::Pos   :
    if (status && particles->pos.size()>0) {
      *data = &particles->pos[first*3];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Vel  :
    if (status && particles->vel.size()>0) {
      *data = &particles->vel[first*3];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Mass  :
    if (status && particles->mass.size()) {
      *data = &particles->mass[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;

  case uns::Rho :
    if (status && comp=="gas" && particles->rho.size()>0) {
      *data = &particles->rho[0];
      *n=particles->rho.size();
    } else {
      ok=false;
    }
    break;
  case uns::U :
      ok=false;
    break;
  case uns::Hsml :
    if (status && comp=="gas" && particles->hsml.size()>0) {
      *data = &particles->hsml[0];
      *n=particles->hsml.size();
    } else {
      ok=false;
    }
    break;
  case uns::Temp :
    if (status && comp=="gas" && particles->temp.size()>0) {
      *data = &particles->temp[0];
      *n=particles->temp.size();
    } else {
      ok=false;
    }
    break;
  case uns::Age :
    if (status && comp=="stars" && particles->age.size()>0) {
      *data = &particles->age[0];
      *n=particles->age.size();
    } else {
      ok=false;
    }
    break;
  case uns::Metal :
    if (status && particles->metal.size()>0) {
      *data = &particles->metal[first];
      *n    = nbody;
    } else {
      ok=false;
    }
    break;

  default: ok=false;
  }
  if (ok && !*data &&
      (CunsOut::s_mapStringValues[name]!=uns::Nbody &&
       CunsOut::s_mapStringValues[name]!=uns::Nsel)) ok = false; // not ok because array is NULL
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotGadgetIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotGadgetIn::getData Value ["<<name<<"] for component <"<<comp<<"> does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData
// return requested array according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string name,int *n,float **data)
{
  bool ok=true;
  *data = NULL;
  *n = 0;
  int first=0;
  bool status=true;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Nsel   :
    if (status) {
      *n    = particles->ntot;;
    } else {
      ok=false;
    }
    break;
  case uns::Pos   :
    if (particles->pos.size()>0) {
      *data = &particles->pos[first*3];
      *n    = particles->pos.size()/3;
    } else {
      ok=false;
    }
    break;
  case uns::Vel  :
    if (particles->vel.size()>0) {
      *data = &particles->vel[first*3];
      *n    = particles->vel.size()/3;
    } else {
      ok=false;
    }
    break;
  case uns::Mass  :
    if (particles->mass.size()) {
      *data = &particles->mass[first];
      *n    = particles->mass.size();
    } else {
      ok=false;
    }
    break;
  case uns::Rho :
    if (particles->rho.size()>0) {
      *data = &particles->rho[0];
      *n=particles->rho.size();
    } else {
      ok=false;
    }
    break;
  case uns::U :
      ok=false;
    break;
  case uns::Hsml :
    if (particles->hsml.size()>0) {
      *data = &particles->hsml[0];
      *n=particles->hsml.size();
    } else {
      ok=false;
    }
    break;
  case uns::Temp :
    if (particles->temp.size()>0) {
      *data = &particles->temp[0];
      *n=particles->temp.size();
    } else {
      ok=false;
    }
    break;
  case uns::Age :
    if (particles->age.size()>0) {
      *data = &particles->age[0];
      *n=particles->age.size();
    } else {
      ok=false;
    }
    break;

  default: ok=false;
  }
  if (ok && !*data &&
      (CunsOut::s_mapStringValues[name]!=uns::Nbody &&
       CunsOut::s_mapStringValues[name]!=uns::Nsel)) ok = false; // not ok because array is NULL
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotGadgetIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotGadgetIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData
// return requested float according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string name,float *data)
{
  bool ok=true;
  *data=0.0;
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Time   :
    *data = amr->getHeader()->time; // find time
    break;

  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotGadgetIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    }  else {
      std::cerr << "**WARNING** CSnapshotGadgetIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData
// return requested array according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string comp,const std::string name,int *n, int **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;

  int nbody,first,last;
  bool status=getRangeSelect(comp.c_str(),&nbody,&first,&last,false); // find components ranges
  if (!status && comp=="all") { // retreive all particles selected by the user
    status=1;
    first=0;
    nbody=particles->ntot;
  }
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Id :
    if (status && particles->id.size()>0) {
      *data = &particles->id[first];
      *n = nbody;
    } else {
      ok = false;
    }
    break;
  case uns::Nbody :
    if (status) {
      *data = NULL;
      *n = nbody;
    } else {
      ok = false;
    }
    break;
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotGadgetIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotGadgetIn::getData Value ["<<name<<"] for component <"<<comp<<"> does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData
// return requested array according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string name,int *n, int **data)
{
  bool ok=false;
  return ok;
}
// ============================================================================
// getData
// return requested int according 'name' selection
bool CSnapshotRamsesIn::getData(const std::string name,int *data)
{
  bool ok=true;
  *data=0;
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Nsel   :
    *data = particles->ntot;
    break;
  case uns::Ngas   :
    *data = particles->ngas;
    break;
  case uns::Nhalo   :
    *data = particles->ndm;
    break;
  case uns::Nstars   :
    *data = particles->nstars;
    break;
  default: ok=false;
  }
  if (ok && !*data) ok = false; // not ok because array is NULL
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotGadgetIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotGadgetIn::getData Value ["<<name<<"] does not exist or empty\n";
    }
  }
  return ok;
}
// ============================================================================
// Normal function
// ============================================================================

} // end of namespace uns
