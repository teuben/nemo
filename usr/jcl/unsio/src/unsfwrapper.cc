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

#include "unsidentifier.h"
#include "userselection.h"
#include "ctools.h"
#include <iostream>
#include <assert.h>
#include <cstdlib>

#define DEBUG 0
#include "unsdebug.h"

namespace uns {
  extern "C" {
    // load protocol
    int uns_init_(const char *,const char *, const char * );
    int uns_load_(const int *);
    int uns_close_(const int *);
    // load data
    int uns_get_nbody_( const int * id, int   * nbody            );
    int uns_get_time_ ( const int * id, float * timu             );
    int uns_get_pos_  ( const int * id, float * pos   , int * size);
    int uns_get_vel_  ( const int * id, float * vel   , int * size);
    int uns_get_mass_ ( const int * id, float * mass  , int * size);
    int uns_get_age_  ( const int * id , float * age  , int * size);
    int uns_get_metal_( const int * id , float * metal, int * size);
    int uns_get_metal_gas_( const int * id , float * metal, int * size);
    int uns_get_metal_stars_( const int * id , float * metal, int * size);
    int uns_get_range_( const int * id, const char * comp ,
			int * nbody, int * first, int * last);
    int uns_get_eps_  ( const int * id, const char * comp, float * eps);
    int uns_sim_type_ ( const int * id);
    void uns_sim_dir_ ( const int * id, char * simdir, int len);
    int uns_get_u_    ( const int * id , float * u  , int * size);
    int uns_get_temp_ ( const int * id , float * temp  , int * size);
    int uns_get_rho_  ( const int * id , float * rho  , int * size);
    int uns_get_hsml_ ( const int * id , float * hsml , int * size);
    int uns_get_cod_  ( const int * id, const char * select, const float * time,
		       float * tcod);
    int uns_get_array_f_( const int * id, const char * comp, const char * tag, float * array , int * size);
    int uns_get_array_i_( const int * id, const char * comp, const char * tag, int   * array , int * size);
    int uns_get_value_f_ ( const int * id, const char * tag, float *  data);
    int uns_get_value_i_ ( const int * id, const char * tag, int   *  data);
    
    // save protocol
    int uns_save_init_(const char *,const char *);
    // save data
    int uns_set_nbody_  ( const int * id, int   * nbody            );
    int uns_set_time_   ( const int * id, float * timu             );
    int uns_set_array_f_( const int * id, const char * comp, const char * tag, float * array , int * size);
    int uns_set_array_i_( const int * id, const char * comp, const char * tag, int   * array , int * size);
    int uns_set_value_f_ ( const int * id, const char * tag, float *  data);
    int uns_set_value_i_ ( const int * id, const char * tag, int   *  data);
    int uns_set_pos_    ( const int * id, float * pos   , int * size);
    int uns_set_vel_    ( const int * id, float * vel   , int * size);
    int uns_set_mass_   ( const int * id, float * mass   , int * size);
    int uns_set_ncomp_  ( const int * id, const char * comp, int * nbody);
    int uns_save_       ( const int * id);
  }

void getCrv(const int ident);
int getUnsvIndex(const int ident);
void checkFArray(const int size, const int nbody);
// ----------------------------------------------------------------------------
// C++ / Fortran interface implementation                                      
// ----------------------------------------------------------------------------
  uns::CunsIdentifierVector unsv;

static bool first=true;
static int ident=10;
ComponentRangeVector * crv;
uns::UserSelection user_select;

// ----------------------------------------------------------------------------
// uns_init_:                                                                  
// initialyse UNS object with a new snapshot                                   
// Input :                                                                     
// simname     -> snapshot's name or simulation's name                         
// select_comp -> selected components or range                                 
// select_time -> selected time                                                
// Output :                                                                    
int uns_init_(const char * simname,const char * select_comp, const char * select_time )
{
  first=false;
  int ret=0;
  
  uns::CunsIn * uns = new uns::CunsIn(simname,select_comp,select_time,false);
  bool valid = uns->isValid();
  
  if (valid) {
    uns::CunsIdentifier * unsi = new  uns::CunsIdentifier();
    unsi->ident = ident++;
    unsi->obj   = uns;
    unsv.push_back(*unsi);
    ret = unsi->ident;
  }
  return ret;
}
// ----------------------------------------------------------------------------
// uns_load                                                                    
// load the next time steps of the snapshot with the idendity "ident           
// return value:                                                               
// -1 snapshot was not open                                                    
//  0 end of snapshot reached                                                  
//  1 sucessfull                                                               
int uns_load_(const int * ident)
{
  PRINT("uns_load ident requested : "<< *ident << "\n";)
  // look for "ident" exist in the already opened snapshots
  // return the index in the vector of -1 if false         
  int index=uns::CunsIdentifier::getUnsvIndex(*ident,&unsv);
  PRINT("index in UNSV ="<< index << "\n";)

  if (index >= 0) {
    uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
    index = snap->nextFrame();
  }

  return index;
}
// ----------------------------------------------------------------------------
// uns_close                                                                    
// close the UNS file with the identifier ident 
int uns_close_(const int * ident)
{
  PRINT("uns_close ident requested : "<< *ident << "\n";)
  // look for "ident" exist in the already opened snapshots
  // return the index in the vector of -1 if false         
  int index=uns::CunsIdentifier::getUnsvIndex(*ident,&unsv);
  PRINT("index in UNSV ="<< index << "\n";)

  if (index >= 0) {
    uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
    snap->close();
    delete ((CunsIn *)unsv[index].obj);
  }

  return index;
}
// ----------------------------------------------------------------------------
// uns_get_array_f:                                                                
// return component's array belonging to the simulation with the identifier ident         
int uns_get_array_f_( const int * id, const char * _comp, const char * _tag, float * array , int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  std::string tag  = tools::Ctools::fixFortran(_tag);
  std::string comp = tools::Ctools::fixFortran(_comp);  
  float *data;
  int nbody;
  bool ok=snap->getData(comp,tag,&nbody,&data);
  int status=0;
  if (ok) {
    int dim=1;
    if (tag=="pos" || tag == "vel" || tag == "acc") dim=3;
    checkFArray(*size*dim,nbody*dim);
    memcpy(array,data,sizeof(float)*nbody*dim);
    status=nbody;
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_get_array_i:                                                                
// return component's array belonging to the simulation with the identifier ident         
int uns_get_array_i_( const int * id, const char * _comp, const char * _tag, int * array , int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  std::string tag  = tools::Ctools::fixFortran(_tag);
  std::string comp = tools::Ctools::fixFortran(_comp);  
  float *data;
  int nbody;
  bool ok=snap->getData(comp,tag,&nbody,&data);
  int status=0;
  if (ok) {
    checkFArray(*size,nbody);
    memcpy(array,data,sizeof(float)*nbody);
    status=nbody;
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_get_value_f_:                                                               
// return tag's value belonging to the simulation with the identifier ident        
int uns_get_value_f_( const int * id, const char * _tag, float   *  data)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  std::string tag = tools::Ctools::fixFortran(_tag);
  int status=0;
  bool ok=snap->getData(tag,data);
  if (ok) {
    status=1;
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_get_value_i_:                                                               
// return tag's value belonging to the simulation with the identifier ident        
int uns_get_value_i_( const int * id, const char * _tag, int   *  data)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  std::string tag = tools::Ctools::fixFortran(_tag);
  int status=0;
  bool ok=snap->getData(tag,data);
  if (ok) {
    status=1;
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getNbody:                                                               
// return #bodies belonging to the simulation with the identifier ident        
int uns_get_nbody_( const int * id, int   * nbody)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  *nbody = snap->getNSel();//user_select.getNSel(); 
  return *nbody;
}
// ----------------------------------------------------------------------------
// uns_getTime:                                                                
// return time belonging to the simulation with the identifier ident           
int uns_get_time_( const int * id, float  * time)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  snap->getData("time",time);
  return 1;
}
// ----------------------------------------------------------------------------
// uns_getPos:                                                                 
// return positions belonging to the simulation with the identifier ident      
int uns_get_pos_( const int * id, float  * pos, int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  float *data;
  int nbody;
  bool ok=snap->getData("pos",&nbody,&data);
  assert(ok==true);
  checkFArray(*size,nbody);
  memcpy(pos,data,sizeof(float)*3*nbody);
  return 1;
}
// ----------------------------------------------------------------------------
// uns_getVel:                                                                 
// return velocities belonging to the simulation with the identifier ident     
int uns_get_vel_( const int * id, float  * vel, int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  float *data;
  int nbody;
  bool ok=snap->getData("vel",&nbody,&data);
  assert(ok==true);
  checkFArray(*size,nbody);
  memcpy(vel,data,sizeof(float)*3*nbody);
  return 1;
}
// ----------------------------------------------------------------------------
// uns_getMass:                                                                
// return masses belonging to the simulation with the identifier ident         
int uns_get_mass_( const int * id, float  * mass, int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  float *data;
  int nbody;
  bool ok=snap->getData("mass",&nbody,&data);
  assert(ok==true);
  checkFArray(*size,nbody);
  memcpy(mass,data,sizeof(float)*nbody);
  return 1;
}
// ----------------------------------------------------------------------------
// uns_getAge:                                                                 
// give ages of the stars belonging to the simulation with the identifier ident
// return #stars                                                               
int uns_get_age_( const int * id, float  * age, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("age",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(age,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getMetal:                                                               
// give metallicity of the gas+stars belonging to the simulation with the      
// identifier ident                                                            
// return #(gas+stars)
int uns_get_metal_( const int * id, float  * metal, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("metal",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(metal,data,sizeof(float)*nbody);
  }

  return status;
}
// ----------------------------------------------------------------------------
// uns_getMetalGas:                                                            
// give metallicity of the gas belonging to the simulation with the            
// identifier ident                                                            
// return #(gas+stars)
int uns_get_metal_gas_( const int * id, float  * metal, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("gas_metal",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(metal,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getMetalStars:                                                          
// give metallicity of the stars belonging to the simulation with the          
// identifier ident                                                            
// return #(gas+stars)
int uns_get_metal_stars_( const int * id, float  * metal, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("stars_metal",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(metal,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getU:                                                                   
// give Internal energy of the gas belonging to the simulation with the        
// identifier ident                                                            
// return #(gas)
int uns_get_u_( const int * id, float  * u, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("u",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(u,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getTemp:                                                                
// give temperature of the gas belonging to the simulation with the            
// identifier ident                                                            
// return #(gas)                                                               
int uns_get_temp_( const int * id, float  * temp, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("temp",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(temp,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getRho:                                                                 
// give density of the gas belonging to the simulation with the                
// identifier ident                                                            
// return #(gas)                                                               
int uns_get_rho_( const int * id, float  * rho, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("rho",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(rho,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_getHsml:                                                                 
// give Hsml of the gas belonging to the simulation with the                
// identifier ident                                                            
// return #(gas)                                                               
int uns_get_hsml_( const int * id, float  * hsml, int * size)
{
  int status=0;
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  // populate user_select object according to the selection
  // and the object
  int nbody;
  float * data;
  bool ok=snap->getData("hsml",&nbody,&data);
  if (ok) {
    status=nbody;
    checkFArray(*size,nbody);
    memcpy(hsml,data,sizeof(float)*nbody);
  }
  return status;
}
// ----------------------------------------------------------------------------
// uns_get_range:                                                              
// return nbody, first and last particles of the component selected and        
// belonging to the simulation with the identifier ident                       
int uns_get_range_( const int * id, const char * _comp ,
			int * nbody, int * first, int * last)
{
  int index=getUnsvIndex(*id);
  std::string comp = tools::Ctools::fixFortran(_comp);
  int status=((CunsIn *)unsv[index].obj)->snapshot->getRangeSelect(comp.c_str(),nbody,first,last,true);
  return status;
}
// ----------------------------------------------------------------------------
// uns_get_eps:                                                                
// return softening of the component selected and                              
// belonging to the simulation with the identifier ident                       
int uns_get_eps_( const int * id, const char * _comp ,
		    float * eps)
{
  int index=getUnsvIndex(*id);
  std::string comp = tools::Ctools::fixFortran(_comp);
  *eps=((CunsIn *)unsv[index].obj)->snapshot->getEps(comp);
  return (((*eps)>0.)?1:0);
}
// ----------------------------------------------------------------------------
// uns_get_cod:                                                                
// return time+cod of the components selected and                              
// belonging to the simulation with the identifier ident                       
int uns_get_cod_ ( const int * id, const char * _select, const float * time,
		   float * tcod)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceIn * snap = ((CunsIn *)unsv[index].obj)->snapshot;
  std::string select = tools::Ctools::fixFortran(_select);
  int status = snap->getCod(select,*time,tcod);
  return status;
}

// ----------------------------------------------------------------------------
// uns_sim_type:                                                               
// return simulation type (nemo|gadget etc....) of the requested simulation    
// and belonging to the simulation with the identifier ident                   
int uns_sim_type_(const int * id)
{
  int index=getUnsvIndex(*id);
  return(((CunsIn *)unsv[index].obj)->snapshot->getInterfaceIndex());
}
// ----------------------------------------------------------------------------
// uns_sim_dir:                                                               
// return simulation dirname of the requested simulation    
// and belonging to the simulation with the identifier ident                   
void uns_sim_dir_(const int * id, char * simdir, int lenstring)
{
  int index=getUnsvIndex(*id);
  std::string dir= ((CunsIn *)unsv[index].obj)->snapshot->getSimDir();
  assert(dir.length() <= (unsigned int) lenstring);
  strcpy(simdir,dir.c_str());
  for (int i=strlen(simdir); i<lenstring; i++) {
    simdir[i] = ' ';
  }
}
// ----------------------------------------------------------------------------
// getCrv:                                                                     
// return the Component Range Vector belonging to the UNS's object "index"     
void getCrv(int index) 
{
  assert(index<(int)unsv.size());
  crv = ((CunsIn *)unsv[index].obj)->snapshot->getSnapshotRange();
  //ComponentRange::list(crv);
}
// ----------------------------------------------------------------------------
// getUnsvIndex:                                                               
// return the UNS's object index of the ident simulation                       
int getUnsvIndex(const int ident) 
{
  // look for "ident" exist in the already opened snapshots
  // return the index in the vector of -1 if false         
  int index=uns::CunsIdentifier::getUnsvIndex(ident,&unsv);
  if (index < 0) {
    std::cerr << "\n\nIdentifier #["<<ident<<"] does not exist\n\n";
    std::cerr << "Aborting..........\n\n";
    std::exit(1);
  }
  return index;
}
// ----------------------------------------------------------------------------
// checkFArray                                                                 
// check if fortran array are big enough to store data                         
void checkFArray(const int size, const int nbody)
{
  if (nbody > size) {
    std::cerr << "Your fortran array size["<<size<<"] is too small to handle\n"
	      << "all the nbody["<<nbody<<"] particles.\n"
	      << "\nprogram aborted....\n";
    std::exit(1);
  }
}

// ----------------------------------------------------------------------------
//             S A V I N G         O P E R A T I O N
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// uns_save_init_:                                                                  
// initialyse UNS object for saving with a new snapshot                                   
// Input :                                                                     
// simname     -> snapshot's output file name
// simtype     -> simulation type
// Output :                                                                    
int uns_save_init_(const char * simname,const char * simtype)
{
  int ret=0;
  
  uns::CunsOut * uns = new uns::CunsOut(simname,simtype,false);
    
  uns::CunsIdentifier * unso = new  uns::CunsIdentifier();
  unso->ident = ident++;
  unso->obj   = uns;
  unsv.push_back(*unso);
  ret = unso->ident;

  return ret;
}

// ----------------------------------------------------------------------------
// uns_set_pos_:                                                                 
// set positions belonging to the simulation with the identifier ident      
int uns_set_pos_( const int * id, float  * pos, int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;

  // save positions
  int status=snap->setData("pos",*size,pos,false);
  return status;
}
// ----------------------------------------------------------------------------
// uns_set_time_:                                                                 
// set time belonging to the simulation with the identifier ident      
int uns_set_time_( const int * id, float  * time)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;

  // save positions
  int status=snap->setData("time",*time);
  return status;
}
// ----------------------------------------------------------------------------
// uns_set_array_f_:                                                                 
// set comp's float array belonging to the simulation with the identifier ident      
int uns_set_array_f_( const int * id, const char * _comp, const char * _array, float * data , int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;

  std::string comp  = tools::Ctools::fixFortran(_comp);
  std::string array = tools::Ctools::fixFortran(_array);
    
  // set data
  int status=snap->setData(comp,array,*size,data,true);  
  return status;
}
// ----------------------------------------------------------------------------
// uns_set_array_i_:                                                                 
// set comp's integer array belonging to the simulation with the identifier ident      
int uns_set_array_i_( const int * id, const char * _comp, const char * _array, int * data , int * size)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;

  std::string comp  = tools::Ctools::fixFortran(_comp);
  std::string array = tools::Ctools::fixFortran(_array);
  
  // set data
  int status=snap->setData(comp,array,*size,data,true);
  return status;
}
// ----------------------------------------------------------------------------
// uns_set_value_f_:                                                                 
// set comp's float data belonging to the simulation with the identifier ident      
int uns_set_value_f_( const int * id, const char * tag, float * data)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;
  std::string comp  = tools::Ctools::fixFortran(tag);
  // set data
  int status=snap->setData(comp,*data);
  return status;
}
// ----------------------------------------------------------------------------
// uns_set_value_i_:                                                                 
// set comp's integer data belonging to the simulation with the identifier ident      
int uns_set_value_i_( const int * id, const char * tag, int * data)
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;

  std::string comp  = tools::Ctools::fixFortran(tag);  
  // set data
  int status=snap->setData(comp,*data);
  return status;
}
// ----------------------------------------------------------------------------
// uns_save_:                                                                 
// save data belonging to the simulation with the identifier ident      
int uns_save_     ( const int * id) 
{
  int index=getUnsvIndex(*id);
  uns::CSnapshotInterfaceOut * snap = ((CunsOut *)unsv[index].obj)->snapshot;
  
  snap->save();
  return 1;
}
} // namespace uns
// ----------------------------------------------------------------------------
// End of file                                                                 
// ---------------------------------------------------------------------------- 
