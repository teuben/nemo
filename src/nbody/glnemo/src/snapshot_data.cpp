// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
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
			   const char * selt, const bool load_vel_,
                           pthread_mutex_t * _mutex_data):VirtualData()
{
  nemo_file   = new char[strlen(filename)+1];
  strcpy(nemo_file,filename);
  select_part = selp;
  select_time = selt;
  load_vel    = load_vel_;
  mutex_data  = _mutex_data;
  
  // make a copy of the select_string bc it's overwriten during io_nemo
  sel2 = new char[strlen(selp)+1];
  strcpy(sel2,selp);
  
  part_data->pos   = NULL;
  PRINT_D std::cerr << "SnapshotData pos="<<part_data->pos <<"\n";
  part_data->nbody = NULL;
  part_data->vel   = NULL;
  part_data->timu  = NULL;
  part_data->nemobits = NULL;
  part_data->vel_norm = NULL;
  is_open = FALSE;
  is_parsed = FALSE;
  is_loading_thread = FALSE;
  is_new_data_loaded = FALSE;
}
// ============================================================================
// Destructor                                                                  
SnapshotData::~SnapshotData()
{
  delete [] sel2;
  close();
  delete [] nemo_file;
  delete part_data;
}
// ============================================================================
// SnapshotData::isValidData()                                                 
// return true if it's a NEMO snapshot. Standard input '-' is assumed to be a  
// valid NEMO snapshot.                                                        
bool SnapshotData::isValidData()
{
  bool status;
  
  if ( ! strcmp(nemo_file,(char *)"-")) {
    return true; // we assume here that - is standard input and a NEMO stream...
  }
  stream str=stropen(nemo_file,(char *)"r"); // open NEMO file for reading
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
int SnapshotData::reload(ParticlesSelectVector * psv, const bool load_vel)
{
  int ret=0;
  if (close()) {
    is_open            = FALSE;
    is_parsed          = FALSE;
    is_loading_thread  = FALSE;
    is_new_data_loaded = FALSE;
    is_end_of_data     = FALSE;
    ret = loadPos(psv, load_vel);
  }
  else {
    ret = -2;
  }
  return ret;
}
// ============================================================================
// SnapshotData::loadPos()                                                     
// Read snapshot positions, send them to the glbox                             
int SnapshotData::loadPos(ParticlesSelectVector * psv, const bool load_vel)
{
  
  //                                                         
  // we must protect data access during multithreaded loading
  //                                                         
  pthread_mutex_lock(mutex_data);
  PRINT_D std::cerr << "loadPos pos="<<part_data->pos <<" nemo_file ["<<nemo_file<<"]\n";
  PRINT_D std::cerr << "loadPos select_part="<<select_part<<"\n"; 
  
  // load via 'io_nemo'
  int status;
  if (load_vel) {
    status=io_nemo(nemo_file,"float,read,sp,n,pos,vel,t,st,b",
		   select_part,&part_data->nbody,&part_data->pos,&part_data->vel,
                   &part_data->timu, select_time,&part_data->nemobits);
  }
  else {
    status=io_nemo(nemo_file,"float,read,sp,n,pos,t,st,b",
		  select_part,&part_data->nbody,&part_data->pos, &part_data->timu,
                  select_time,&part_data->nemobits);
  }
  is_open=TRUE;
  if (status != 0) {
    if ( status == -1) {  // Bad nemobits
      if ( ! ( *part_data->nemobits & TimeBit)) { // no TimeBit
        if ( ! part_data->timu ) {   
          part_data->timu = (float *) malloc(sizeof(float));          
        }
        std::cerr << "Forcing time to [0.0]\n";
        *(part_data->timu) = 0.0;
      }
    }
    PRINT_D std::cerr << "(0)time = " << *(part_data->timu) << "\n";
    // intialyse tree_depht array
    part_data->allocTree();

    PRINT_D std::cerr << "time = " << *(part_data->timu) << "\n";        
    is_new_data_loaded = TRUE;
    computeCooMax();
    float * cmax = getCooMax();
    PRINT_D std::cerr << "COOmax = " <<cmax[0]
                      <<" "<<cmax[1]<<" "<<cmax[2]<<"\n";
    full_nbody=*(part_data->nbody);
    PRINT_D std::cerr << "part_data->nbody = " << *(part_data->nbody) << "\n";
    if (! is_parsed) {
      is_parsed = TRUE;
      if ( ! strcmp(nemo_file,"-")) {
          full_nbody=*(part_data->nbody);
      }
      // particles selected from a RANGE
      //!!!int nobject=VirtualData::fillParticleRange(psv,full_nbody,sel2);
      //!!!if (nobject) ; // do nothing (remove compiler warning)
    }
    ParticlesSelectVector psv2=*psv;
    int nb_select2=VirtualParticlesSelect::nb_select;
    VirtualParticlesSelect::nb_select=0;
    psv->clear();
    int nobject=VirtualData::fillParticleRange(psv,full_nbody,sel2);
    if (nobject) {;} // do nothing (remove compiler warning)
    // copy back color and visibility
    for (unsigned int i=0; i<psv2.size();i++) {
      if (i<psv->size()) {
          (*psv)[i].vps->col = psv2[i].vps->col;
          (*psv)[i].vps->is_visible = psv2[i].vps->is_visible;
      }
    }
    VirtualParticlesSelect::nb_select=nb_select2;
    //pthread_mutex_lock(mutex_data);
    // compute velocity vector norm
    part_data->computeVelNorm();
    // we can now unlock the data
    pthread_mutex_unlock(mutex_data);
    
    if (! is_loading_thread) {
      // a running Thread MUST not emit data to the GLBox     
      // because it's not possible to share the OpenGl Display
      // it crashs the aplication                             
      //emit newTime(*timu);             // send to animation   
      emit loadedData(part_data, psv);  // send to glbox
      //is_new_data_loaded = FALSE;
    }    
  } else {
    // we can now unlock the data
    pthread_mutex_unlock(mutex_data);
    is_end_of_data = TRUE;
  }
  return status;
}
// ============================================================================
// SnapshotData::uploadGlData()                                                
// Upload Data to GLbox                                                        
void SnapshotData::uploadGlData(ParticlesSelectVector * psv)
{
  if (1||is_new_data_loaded) {
    //emit newTime(*timu);               // send to animation 
    emit loadedData(part_data, psv);    // send to glbox
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
