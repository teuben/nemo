// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// SnapshotData Class implementation                                           
//                                                                             
//                                                                             
// ============================================================================

extern "C" { // Must include NEMO's header before calling [using namespace std]
             // to avoid confict on 'string' in <stdinc.h>
#include <stdinc.h>
#include <filestruct.h>
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>
#include <snapshot/snapshot.h>
}
#include "snapshot_data.h"
#include <iostream>
#include <qmessagebox.h>

#define LOCAL_DEBUG 0
#include "print_debug.h"

using namespace std;

// ============================================================================
// Constructor                                                                 
SnapshotData::SnapshotData(const char * filename, const char * selp, 
                           const char * selt):VirtualData()
{
  nemo_file   = new char[strlen(filename)+1];
  strcpy(nemo_file,filename);
  select_part = selp;
  select_time = selt;
  
  // make a copy of the select_string bc it's overwriten during io_nemo
  sel2 = new char[strlen(selp)+1];
  strcpy(sel2,selp);
  
  pos   = NULL;
  PRINT_D std::cerr << "SnapshotData pos="<<pos <<"\n";
  nbody = NULL;
  timu  = NULL;
  nemobits = NULL;
  is_open = FALSE;
  is_parsed = FALSE;
  is_loading_thread = FALSE;
  is_new_data_loaded = FALSE;
}
// ============================================================================
// Destructor                                                                  
SnapshotData::~SnapshotData()
{
  if (pos) free((float *) pos);
  if (nbody) free((int *) nbody);
  delete [] sel2;
  close();
  delete [] nemo_file;
  if (nemobits) delete nemobits;
}
// ============================================================================
// SnapshotData::isValidData()                                                 
// return true if it's a NEMO snapshot. Standard input '-' is assumed to be a  
// valid NEMO snapshot.                                                        
bool SnapshotData::isValidData()
{
  bool status;
  
  if ( ! strcmp(nemo_file,"-")) {
    return true; // we assume here that - is standard input and a NEMO stream...
  }
  stream str=stropen(nemo_file,"r"); // open NEMO file for reading
  if ( ! str ) {
    status=FALSE;
  }
  if (qsf(str)) { // is it a structured binary file (NEMO?)
    status = TRUE;
  }
  else {
    status = FALSE;
  }
  strclose(str);
  if (status)  { // it's a true snapshot
    // get the full nbody
    int * ptr=&full_nbody;
    if (io_nemo(nemo_file,"read,n",&ptr) > 0) {
      io_nemo(nemo_file,"close");
    } 
  }
  return status;
}
// ============================================================================
// SnapshotData::close()                                                       
// close nemo snapshot                                                         
int SnapshotData::close()
{
  int status=0;
  if (is_open) {
    status = io_nemo(nemo_file,"close");
  }
  return status;
}

// ============================================================================
// SnapshotData::reload()                                                      
// load the snapshot from the beginning                                        
int SnapshotData::reload(ParticlesSelectVector * psv)
{
  int ret=0;
  if (close()) {
    is_open            = FALSE;
    is_parsed          = FALSE;
    is_loading_thread  = FALSE;
    is_new_data_loaded = FALSE;
    is_end_of_data     = FALSE;
    ret = loadPos(psv);
  }
  else {
    ret = -2;
  }
  return ret;
}
// ============================================================================
// SnapshotData::loadPos()                                                     
// Read snapshot positions, send them to the glbox                             
int SnapshotData::loadPos(ParticlesSelectVector * psv)
{
  PRINT_D std::cerr << "loadPos pos="<<pos <<" nemo_file ["<<nemo_file<<"]\n";
  PRINT_D std::cerr << "loadPos select_part="<<select_part<<" nbody="<<nbody<<"\n";
  // load via 'io_nemo'
  int status=io_nemo(nemo_file,"float,read,sp,n,pos,t,st,b",
             select_part,&nbody,&pos, &timu, select_time,&nemobits);  
  if (status != 0) {
    if ( status == -1) {  // Bad nemobits
      if ( ! ( *nemobits & TimeBit)) { // no TimeBit
        if ( ! timu ) {   
          timu = (float *) malloc(sizeof(float));          
        }
        std::cerr << "Forcing time to [0.0]\n";
        *timu = 0.0;
      }
    } 
    PRINT_D std::cerr << "time = " << *timu << "\n";        
    is_open=TRUE;
    is_new_data_loaded = TRUE;
    computeCooMax();
    float * cmax = getCooMax();
    PRINT_D std::cerr << "COOmax = " <<cmax[0]
                      <<" "<<cmax[1]<<" "<<cmax[2]<<"\n";
    if (! is_parsed) {
      is_parsed = TRUE;
      if ( ! strcmp(nemo_file,"-")) {
          full_nbody=*nbody;
      }
      // particles selected from a RANGE
      int nobject=VirtualData::fillParticleRange(psv,full_nbody,sel2);
      if (nobject) ; // do nothing (remove compiler warning)
    }
    if (! is_loading_thread) {
      // a running Thread MUST not emit data to the GLBox     
      // because it's not possible to share the OpenGl Display
      // it crashs the aplication                             
      //emit newTime(*timu);             // send to animation   
      emit loadedData(nbody,pos,psv);  // send to glbox       
      is_new_data_loaded = FALSE;
    }    
  } else {
    is_end_of_data = TRUE;
  }
  return status;
}
// ============================================================================
// SnapshotData::uploadGlData()                                                
// Upload Data to GLbox                                                        
void SnapshotData::uploadGlData(ParticlesSelectVector * psv)
{
  if (is_new_data_loaded) {
    //emit newTime(*timu);               // send to animation   
    emit loadedData(nbody,pos,psv);    // send to glbox       
    is_new_data_loaded = FALSE;
  }
}
// ============================================================================
// SnapshotData::endOfDataMessage()                                            
// return end of data message                                                  
QString SnapshotData::endOfDataMessage()
{
  QString nemof=nemo_file;
  QString message="Snapshot ["+nemof+"] end of snapshot reached!";
  return message;
}
// ============================================================================
