// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
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
}
#include "snapshot_data.h"
#include <iostream>

#define LOCAL_DEBUG 0
#include "print_debug.h"

using namespace std;

// ----------------------------------------------------------------------------
// SnapshotData => Constructor
// Initialyze data

SnapshotData::SnapshotData(const char * filename, const char * selp, const char * selt):VirtualData()
{
  nemo_file   = new char[strlen(filename)+1];
  strcpy(nemo_file,filename);
  select_part = selp;
  select_time = selt;
  
  // make a copy of the select_string bc it's overwriten during io_nemo
  sel2 = new char[strlen(selp)+1];
  strcpy(sel2,selp);
  pos   = NULL;
  std::cerr << "SnapshotData pos="<<pos <<"\n";
  nbody = NULL;
  timu  = NULL;
  is_open = FALSE;
  is_parsed = FALSE;
  is_loading_thread = FALSE;
  is_new_data_loaded = FALSE;
}
SnapshotData::~SnapshotData()
{
  delete [] pos;
  delete nbody;
  delete sel2;
  close();
  delete nemo_file;
  //ParticlesRange::nb_select=0;  // Reset 
}
// ----------------------------------------------------------------------------
// isNemo
// return true if it's a NEMO snapshot
bool SnapshotData::isValidData()
{
  bool status;
  
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
// ----------------------------------------------------------------------------
// close
// close nemo snapshot
int SnapshotData::close()
{
  int status=0;
  if (is_open) {
    std::cerr << ">> BUG around this call :Trying to io nemo close\n";
    status = io_nemo(nemo_file,"close");
    std::cerr << "<< BUG around this call :Trying to io nemo close\n";
  }
  return status;
}

// ----------------------------------------------------------------------------
// restart:
// read the snapshot from the beginnig
int SnapshotData::reload(ParticlesRangeVector * prv)
{
  int ret=0;
  if (close()) {
    is_open            = FALSE;
    is_parsed          = FALSE;
    is_loading_thread  = FALSE;
    is_new_data_loaded = FALSE;
    is_end_of_data     = FALSE;
    ret = loadPos(prv);
  }
  else {
    ret = -2;
  }
  return ret;
}
// ----------------------------------------------------------------------------
// loadData
// Read nemo snapshot
int SnapshotData::loadPos(ParticlesRangeVector * prv)
{
  std::cerr << "loadPos pos="<<pos <<" nemo_file ["<<nemo_file<<"]\n";
  std::cerr << "loadPos select_part="<<select_part<<" nbody="<<nbody<<"\n";
  int status=io_nemo(nemo_file,"info,float,read,sp,n,pos,t,st",
             select_part,&nbody,&pos, &timu, select_time);
  std::cerr << "time = " << *timu << "\n";      
  if (status > 0) {
    is_open=TRUE;
    is_new_data_loaded = TRUE;
    computeCooMax();
    float * cmax = getCooMax();
    PRINT_D std::cerr << "COOmax = " <<cmax[0]<<" "<<cmax[1]<<" "<<cmax[2]<<"\n";
    if (! is_parsed) {
      is_parsed = TRUE;
      //int nobject=VirtualData::fillParticleRange(prv,*nbody,sel2);
      int nobject=VirtualData::fillParticleRange(prv,full_nbody,sel2);
      if (nobject) ; // do nothing (remove compiler warning)
    }
    if (! is_loading_thread) {
      // a running Thread MUST not emit data to the GLBox
      // because it's not possible to share the OpenGl Display
      // it crashs the aplication
      emit loadedData(nbody,pos,prv);
      is_new_data_loaded = FALSE;
    }
    
  } else {
    is_end_of_data = TRUE;
    //QString message="End of NEMO snapshot";
    //cerr << "Iam HERE !!!\n";
    //emit messageLoad(&message);
  }
  return status;
}
// ----------------------------------------------------------------------------
// Upload Data to GLbox
// 
void SnapshotData::uploadGlData(ParticlesRangeVector * prv)
{
  if (is_new_data_loaded) {
    emit loadedData(nbody,pos,prv);
    is_new_data_loaded = FALSE;
  }
}

// ----------------------------------------------------------------------------
// spread all the selected particles in separated data structure 
#if 0
int SnapshotData::fillParticleRange(ParticlesRangeVector * prv)
{
  ParticlesRange * pr;
  char * s = sel2;

  while (s) {
    pr = new ParticlesRange();
    s=pr->parseString(s,*nbody,prv);
    if (s) 
      cerr << " >>>> s sring = ["<< s << "]\n";
    prv->push_back(*pr);
    cerr << "In globwin, prv->size() = " << prv->size() << "\n";	  
    delete pr;
  }
  for (int i=0; i< (int) prv->size(); i++) {
    cerr << " - - - - - - - - - - - \n";
    cerr << i << "\n";
    (*prv)[i].printRange();
  }
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //  exit(1); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  return prv->size();
}
#endif
// ----------------------------------------------------------------------------
//
QString SnapshotData::endOfDataMessage()
{
  QString nemof=nemo_file;
  QString message="Snapshot ["+nemof+"] end of snapshot reached!";
  return message;
}
//
