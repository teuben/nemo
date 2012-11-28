// ============================================================================
// Copyright Jean-Charles LAMBERT - 2008-2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/* 
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */

#include <iostream>
#include <cstdlib>
#include <assert.h>
#include <algorithm>
#include <cctype> // for toupper()
#include "snapshotnemo.h"
#include "ctools.h"
#include "uns.h"

#define DEBUG 0
#include "unsdebug.h"

namespace uns {

// ============================================================================
// Constructor                                                                 
CSnapshotNemoIn::CSnapshotNemoIn(const std::string _name,
			     const std::string _comp, 
			     const std::string _time,
			     const bool verb)
  :CSnapshotInterfaceIn(_name, _comp, _time, verb)
{
  const char * defv[] = { "none=none","VERSION=XXX",NULL };
  const char * argv[] = { "CSnapshotNemoIn",NULL };

  interface_type="Nemo";
  file_structure = "range";
  interface_index=0;
  first_stream=false;
  nemobits=NULL;
  ionbody =NULL;
  iotime  =NULL;
  iopos   =NULL;
  iovel   =NULL;
  iomass  =NULL;
  ioacc   =NULL;
  iopot   =NULL;
  ioaux   =NULL;
  iorho   =NULL;
  iokeys  =NULL;
  
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  acc    = NULL;
  pot    = NULL;
  rho    = NULL;
  aux    = NULL;
  keys   = NULL;
  last_nbody=0;
  last_nemobits=-1;
  initparam(const_cast<char**>(argv),const_cast<char**>(defv));
  valid=isValidNemo();
}

// ============================================================================
// Destructor                                                                  
CSnapshotNemoIn::~CSnapshotNemoIn()
{

  if (ionbody) free ((int   *) ionbody);
  if (iotime ) free ((float *) iotime );
  if (iopos  ) free ((float *) iopos  );
  if (iovel  ) free ((float *) iovel  );
  if (iomass ) free ((float *) iomass );
  if (iorho  ) free ((float *) iorho  );
  if (ioaux  ) free ((float *) ioaux  );
  if (ioacc  ) free ((float *) ioacc  );
  if (iopot  ) free ((float *) iopot  );
  if (iokeys ) free ((int   *) iokeys );
  
  if (pos    ) delete [] pos;
  if (vel    ) delete [] vel;
  if (mass   ) delete [] mass;
  if (rho    ) delete [] rho;
  if (aux    ) delete [] aux;
  if (acc    ) delete [] acc;
  if (pot    ) delete [] pot;
  if (keys   ) delete [] keys;
  
  if (valid) close();
}
// ============================================================================
// isValidNemo()                                                               
// return true if it's a NEMO snapshot.                                        
bool CSnapshotNemoIn::isValidNemo()
{
  bool status;
  valid=true;
  
  if (filename == "-") { // we assume here that "-"
    status=true;         // is standard input and a NEMO stream...
    first_stream=true;
    //select_time="all";
    std::string force_select = "all"; 
    status_ionemo=io_nemo(filename.c_str(),"float,read,sp,n,pos,vel,mass,dens,aux,acc,pot,key,t,st,b",
            force_select.c_str(),&ionbody,&iopos,&iovel,&iomass,&iorho,&ioaux,&ioacc,&iopot,&iokeys,
            &iotime, select_time.c_str(),&nemobits);

    full_nbody = *ionbody;        
  }
  else { // Normal file
    stream str=stropen(filename.c_str(),(char *) "r"); // open NEMO file for reading
    if ( ! str )  status = false; // failed to open
    if (qsf(str)) status = true;  // it's a structured binary file (NEMO)
    else          status = false; // it's not                            
    strclose(str);
    if (status)  {                // it's a NEMO snapshot
      int * ptr=NULL;      // get the full nbody
      if (io_nemo(filename.c_str(),"float,read,n,t,b",&ptr,&iotime,&nemobits) != 0) {
        io_nemo(filename.c_str(),"close");
      } else {
      }
      assert(ptr);
      full_nbody=*ptr;
      free((int *) ptr);
    }
  }
  valid=status;
  if (valid) {
    if ( ! ( *nemobits & TimeBit)) { // no TimeBit
      time_first = 0.0;
    }
    else time_first = *iotime;
  }
  return status;
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * CSnapshotNemoIn::getSnapshotRange()
{
  crv.clear();
  if (valid) {
    ComponentRange * cr = new ComponentRange();
    cr->setData(0,full_nbody-1);
    cr->setType("all");
    crv.push_back(*cr);
    //ComponentRange::list(&crv);
    delete cr;
    if (first) {
      first       = false;
      crv_first   = crv;
      nbody_first = full_nbody;
      //time_first  = 0.0;
    }
  }
  return &crv;
}
// ============================================================================
// nextFrame()                                                                 
int CSnapshotNemoIn::nextFrame(uns::UserSelection &user_select)
{
  const uns::t_indexes_tab * index_tab=user_select.getIndexesTab();
  const int nsel=user_select.getNSel();
  
  int status;  // io_nemo status
  std::string force_select = "all"; 
  if (! first_stream) { // normal file or second reading (stream)
    status=io_nemo(filename.c_str(),"float,read,sp,n,pos,vel,mass,dens,aux,acc,pot,key,t,st,b",
                   force_select.c_str(),&ionbody,&iopos,&iovel,&iomass,&iorho,&ioaux,&ioacc,&iopot,&iokeys,
		   &iotime, select_time.c_str(),&nemobits);
  } else { // "-" first stream, no need to read
    first_stream=false;
    status=status_ionemo; // status read the first time cf : isValid()
  }
  if (status == -1) {  // Bad nemobits
    PRINT("io_nemo status="<<status<<"\n";);
  } else {
    PRINT("io_nemo status="<<status<<" time="<<*iotime<<"\n";);
  }
  if ( status != 0) {
    if (status == -1) {  // Bad nemobits
      if ( ! ( *nemobits & TimeBit)) { // no TimeBit
        if ( ! iotime )
          iotime = (float *) malloc(sizeof(float));
        std::cerr << "CSnapshotNemoIn::nextFrame => Forcing time to [0.0]\n";
        *(iotime) = 0.0;
      }
    }
    if (status != -2) { // NEMO snapshot must have particles TAG
      if (*ionbody > last_nbody || (last_nemobits>0 && last_nemobits!=(*nemobits))) { // we must resize arrays
	if (pos) delete [] pos;
	if ( *nemobits & PosBit && req_bits&POS_BIT)  pos = new float[*ionbody*3];
	else pos=NULL;
	if (vel) delete [] vel;
	if ( *nemobits & VelBit && req_bits&VEL_BIT)  vel = new float[*ionbody*3];
	else vel=NULL;
	if (mass) delete [] mass;
	if ( *nemobits & MassBit && req_bits&MASS_BIT)  mass = new float[*ionbody];
	else mass=NULL;
	if (rho) delete [] rho;
	if ( *nemobits & DensBit && req_bits&RHO_BIT)  rho = new float[*ionbody];
	else rho=NULL;
	if (acc) delete [] acc;
	if ( *nemobits & AccelerationBit && req_bits&ACC_BIT)  acc = new float[*ionbody*3];
	else acc=NULL;
	if (aux) delete [] aux;
	if ( *nemobits & AuxBit && req_bits&AUX_BIT)  aux = new float[*ionbody];
	else aux=NULL;
	if (pot) delete [] pot;
	if ( *nemobits & PotentialBit && req_bits&POT_BIT)  pot = new float[*ionbody];
	else pot=NULL;
	if (keys) delete [] keys;
        if ( *nemobits & KeyBit && (req_bits&KEYS_BIT || req_bits&ID_BIT))  keys = new int[*ionbody];
        else keys=NULL;
      }
      last_nbody    = *ionbody;  // save nbody
      last_nemobits = *nemobits; // save nemobits
      int cpt=0;
      // fill array according to selection
      for (int i=0; i<*ionbody; i++) {
	int idx=index_tab[i].i;
	if (idx!=-1) { // index selected
	  for (int j=0; j<3; j++) {
	    if ( *nemobits & PosBit && req_bits&POS_BIT) pos[cpt*3+j] = iopos[idx*3+j];
	    if ( *nemobits & VelBit && req_bits&VEL_BIT) vel[cpt*3+j] = iovel[idx*3+j];
            if ( *nemobits & AccelerationBit && req_bits&ACC_BIT) acc[cpt*3+j] = ioacc[idx*3+j];
	  }
	  if ( *nemobits & MassBit && req_bits&MASS_BIT) mass[cpt] = iomass[cpt];
          if ( *nemobits & DensBit && req_bits&RHO_BIT) rho[cpt]  = iorho[cpt];
          if ( *nemobits & AuxBit  && req_bits&AUX_BIT) aux[cpt]  = ioaux[cpt];
          if ( *nemobits & PotentialBit && req_bits&POT_BIT ) pot[cpt]  = iopot[cpt];
          if ( *nemobits & KeyBit  && (req_bits&KEYS_BIT || req_bits&ID_BIT)) keys[cpt]  = iokeys[cpt];
	  cpt++;
	}
      }
      assert(nsel==cpt);
    }
  }
  if (nsel) ; // remove compiler warning
  if (verbose) std::cerr << "CSnapshotNemoIn::nextFrame status = " << status << "\n";
  if (status == -1) status=1;
  return status;

}
// ============================================================================
// close()                                                                     
int CSnapshotNemoIn::close()
{
  int status=0;
  if (valid) {
    status = io_nemo(filename.c_str(),"close");
    end_of_data = false;
  }
  return status;
}
// ============================================================================
// checkBits                                                                   
void CSnapshotNemoIn::checkBits(std::string comp,const int bits)
{
  if ( ! ( *nemobits & bits)) { // no XXXbits
    std::cerr << "You have requested the component ["<<comp<<"] which is missing\n"
	      << " in the file. Aborting program.....\n\n";
    std::exit(1);
  }
}
// ============================================================================
// getData                                                               
// return requested array according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string name,int *n,float **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;
  
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Pos   :
    *data = getPos();
    *n    = getNSel();
    break;
  case uns::Vel  :
    *data = getVel();
    *n    = getNSel();
    break;
  case uns::Mass  :
    *data = getMass();
    *n    = getNSel();
    break;
  case uns::Acc :
    *data = getAcc();
    *n    = getNSel();
    break;        
  case uns::Rho :
    *data = getRho();
    *n    = getNSel();
    break;
  case uns::Pot :
    *data = getPot();
    *n    = getNSel();
    break;    
  case uns::Aux :
  case uns::Hsml :
    *data = getAux();
    *n    = getNSel();
    break;
    
  default: ok=false;
  }
  if (*data == NULL) {
    ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData                                                               
// return requested array according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string comp,const std::string name,int *n,float **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;
  
  int nbody,first,last;
  bool status=getRangeSelect(comp.c_str(),&nbody,&first,&last,false); // find components ranges
  if (!status && comp=="all") { // retreive all particles selected by the user
    status=1;
    first=0;
    nbody=getNSel();
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
    if (status && getPos()) {
      *data = &getPos()[first*3];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Vel  :
    if (status && getVel()) {
      *data = &getVel()[first*3];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Acc  :
    if (status && getAcc()) {
      *data = &getAcc()[first*3];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;    
  case uns::Pot  :
    if (status && getPot()) {
      *data = &getPot()[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;    
  case uns::Mass  :
    if (status && getMass()) {
      *data = &getMass()[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Rho :
    if (status && getRho()) {
      *data = &getRho()[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
  case uns::Aux :
  case uns::Hsml :
    if (status && getAux()) {
      *data = &getAux()[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
    
  default: ok=false;
  }

  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] for component <"<<comp<<"> does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData                                                               
// return requested float according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string name,float *data)
{
  bool ok=true;
  *data=0.0;  
  switch(CunsOut::s_mapStringValues[name]) {
    case uns::Time   :
     *data = getTime();   
     break;
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData                                                               
// return requested array according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string name,int *n,int **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;
    
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Keys :
  case uns::Id   :
    *data = getKeys();
    *n    = getNSel();
    break;
  default: ok=false;
  }
  if (*data == NULL) {
    ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData                                                               
// return requested array according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string comp,const std::string name,int *n,int **data)
{
  bool ok=true;
  *data=NULL;
  *n = 0;
  
  int nbody,first,last;
  bool status=getRangeSelect(comp.c_str(),&nbody,&first,&last,false); // find components ranges
  if (!status && comp=="all") { // retreive all particles selected by the user
    status=1;
    first=0;
    nbody=getNSel();
  }
  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Id   :
    if (status && getKeys()) {
      *data = &getKeys()[first];
      *n    = nbody;//getNSel();
    } else {
      ok=false;
    }
    break;
    
    default: ok=false;
    }

  if (verbose) {
    if (ok) {      
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================
// getData                                                               
// return requested int according 'name' selection                                               
bool CSnapshotNemoIn::getData(const std::string name, int *data)
{
  bool ok=true;
  *data=0;  
  switch(CunsOut::s_mapStringValues[name]) {
    case uns::Nsel   :
    *data = getNSel();   
     break;
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoIn::getData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoIn::getData Value ["<<name<<"] does not exist...\n";
    }
  }
  return ok;
}
// ============================================================================

//
// CSnapshotNemoOut CLASS implementation
//


// ----------------------------------------------------------------------------
// constructor
CSnapshotNemoOut::CSnapshotNemoOut(const std::string _n, const std::string _t, const bool _v):CSnapshotInterfaceOut(_n, _t, _v)
{
  if (simtype != "nemo") {
    std::cerr << "CSnapshotNemoOut::CSnapshotNemoOut Unkwown file type : ["<<simtype<<"]\n"
              << "aborting .....\n";
    std::exit(1);
  }
  //
  interface_type = "Nemo";
  file_structure = "range"; // "range" like file
  // RAZ
  mass   = NULL;
  pos    = NULL;
  vel    = NULL;
  aux    = NULL;
  acc    = NULL;
  pot    = NULL;
  rho    = NULL;
  keys   = NULL;
  
  ptrIsAlloc["mass" ]=false; 
  ptrIsAlloc["pos"  ]=false; 
  ptrIsAlloc["vel"  ]=false; 
  ptrIsAlloc["pot"  ]=false; 
  ptrIsAlloc["acc"  ]=false; 
  ptrIsAlloc["aux"  ]=false; 
  ptrIsAlloc["keys" ]=false; 
  ptrIsAlloc["rho"  ]=false; 
  ptrIsAlloc["id"   ]=false; 
  
  nbody = -1;
  bits = 0;
  is_saved=false;
  is_closed=false;
}
// ----------------------------------------------------------------------------
// desstructor
CSnapshotNemoOut::~CSnapshotNemoOut()
{
  if (mass && ptrIsAlloc["mass"]) delete [] mass;
  if (pos  && ptrIsAlloc["pos" ]) delete [] pos;
  if (vel  && ptrIsAlloc["vel" ]) delete [] vel;
  if (pot  && ptrIsAlloc["pot" ]) delete [] pot;
  if (acc  && ptrIsAlloc["acc" ]) delete [] acc;
  if (aux  && ptrIsAlloc["aux" ]) delete [] aux;
  if ((keys && ptrIsAlloc["keys"]) ||
      (keys && ptrIsAlloc["id"  ])  ) delete [] keys;
  if (rho  && ptrIsAlloc["rho" ]) delete [] rho;
  close();
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setHeader(void * h)
{
  if (h) {;} // remove compiler warning
  return 0;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setNbody(const int _n)
{
  nbody = _n;
  return nbody;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name, float data) 
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Time : 
    status = 1;
    time = data;
    bits |= TimeBit;
    break;
  default: ok=false;
  }

  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name, const int n ,int * data,const bool _addr)
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Keys :
  case uns::Id:
    //status = setKeys(n, data, _addr);
    status = setArray(n,1,data,&keys,name.c_str(),KeyBit,_addr);
    break;
  default: ok=false;
  }

  if (verbose) {
    if (ok) {      
      std::cerr << "CCSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}

// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name,std::string array,  const int n ,float * data,const bool _addr)
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::All  :
    status = setData(array,n,data,_addr);
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name,std::string array,  const int n ,int * data,const bool _addr)
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::All  :
    status = setData(array,n,data,_addr);
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name,  const int n ,float * data,const bool _addr)
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::Pos  :
    status = setArray(n,3,data,&pos,name.c_str(),PosBit,_addr);
    break;
  case uns::Vel  :
    status = setArray(n,3,data,&vel,name.c_str(),VelBit,_addr);
    break;
  case uns::Mass : 
    status = setArray(n,1,data,&mass,name.c_str(),MassBit,_addr);
    break;
  case uns::Acc :
    status = setArray(n,3,data,&acc,name.c_str(),AccelerationBit,_addr);
    break;
  case uns::Pot :
    status = setArray(n,1,data,&pot,name.c_str(),PotentialBit,_addr);
    break;
  case uns::Aux :
  case uns::Hsml:
    status = setArray(n,1,data,&aux,name.c_str(),AuxBit,_addr);
    break;
  case uns::Rho :
    status = setArray(n,1,data,&rho,name.c_str(),DensBit,_addr);
    break;
  default: ok=false;
  }
  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}
// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::setData(std::string name, const int n ,float * data, float * data1, float * data2, const bool _addr)
{
  bool ok=true;
  int status=0;

  switch(CunsOut::s_mapStringValues[name]) {
  case uns::All :
    status = setArray(n,1,data ,&mass,"mass",MassBit,_addr);
    status = setArray(n,3,data1,&pos ,"pos" ,PosBit ,_addr);
    status = setArray(n,3,data2,&vel ,"vel" ,VelBit ,_addr);
    break;
  default: ok=false;
  }

  if (verbose) {
    if (ok) {
      std::cerr << "CSnapshotNemoOut::setData name["<<name<<"]=" << CunsOut::s_mapStringValues[name] << "\n";
    } else {
      std::cerr << "**WARNING** CSnapshotNemoOut::setData Value ["<<name<<"] does not exist....\n";
      //std::exit(1);
    }
  }
  return status;
}
// ============================================================================
// setArray (float)
int CSnapshotNemoOut::setArray(const int _n, const int dim, float * src, float ** dest, const char * name, const int tbits, const bool addr)
{
  if (nbody<0) {
    nbody = _n;
  } else {
    assert(nbody==_n);
  }
  
  if (addr) { // map address
    *dest = src;
  }
  else {
    ptrIsAlloc[name]=true;
    if (*dest)  delete [] (*dest);
    *dest = new float[_n*dim];
    memcpy(*dest,src,sizeof(float)*_n*dim);
  }
  bits |= tbits;
  return 1;
}
// ============================================================================
// setArray (int)
int CSnapshotNemoOut::setArray(const int _n, const int dim, int * src, int ** dest, const char * name, const int tbits, const bool addr)
{
  if (addr) { // map address
    *dest = src;
  }
  else {
    ptrIsAlloc[name]=true;
    if (*dest)  delete [] (*dest);
    *dest = new int[_n*dim];
    memcpy(*dest,src,sizeof(int)*_n*dim);
  }
  bits |= tbits;
  return 1;
}

// ----------------------------------------------------------------------------
//
int CSnapshotNemoOut::save()
{
  int   * n = &nbody;
  float * t = &time;
  int   * b = &bits;
  int status=io_nemo(simname.c_str(),"float,save,n,t,x,v,m,p,a,aux,k,dens,b",
                     &n,&t,&pos,&vel,&mass,&pot,&acc,&aux,&keys,&rho,&b);
  if (status!=0) {
    is_saved=true;
  }
  return status;
}
// ============================================================================
// close()                                                                     
int CSnapshotNemoOut::close()
{
  int status=0;
  if (is_saved && !is_closed) {
    is_closed = true;
    status = io_nemo(simname.c_str(),"close");    
  }
  return status;
}
// ============================================================================
// moveToCom()                                                                     
std::vector<double> CSnapshotNemoOut::moveToCom()
{
  std::vector<double> com(6,0.);
  double masstot=0.0;
  // compute center of mass
  for (int i=0; i<nbody;i++) {
    float massi;
    if (mass) massi = mass[i];
    else      massi = 1.0;
    masstot+=massi;
    if (pos) {
      com[0] +=(pos[i*3  ]*massi);
      com[1] +=(pos[i*3+1]*massi);
      com[2] +=(pos[i*3+2]*massi);
    }
    if (vel) {
      com[3] +=(vel[i*3  ]*massi);
      com[4] +=(vel[i*3+1]*massi);
      com[5] +=(vel[i*3+2]*massi);
    }
  }
  if (!mass) {
    std::cerr << "CSnapshotNemoOut::moveToCom => No mass in the snapshot, we assum mass=1.0 for each particles...\n";
  }
  // shift to center of mass
  for (int i=0; i<nbody;i++) {
    if (pos) {
      pos[i*3+0] -= (com[0]/masstot);
      pos[i*3+1] -= (com[1]/masstot);
      pos[i*3+2] -= (com[2]/masstot);
    }
    if (vel) {
      vel[i*3+0] -= (com[3]/masstot);
      vel[i*3+1] -= (com[4]/masstot);
      vel[i*3+2] -= (com[5]/masstot);
    }
  }  
  return com;
}
    
} // end of uns namespace

//
