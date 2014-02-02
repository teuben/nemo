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
#include "loadingthread.h"
#include <assert.h>

namespace glnemo {

// ============================================================================
// constructor                                                                 
LoadingThread::LoadingThread(SnapshotInterface * si,
                             UserSelection * _user_select,
                             ParticlesObjectVector * _pov,
                             const std::string _select, QMutex * _mutex,
                             GlobalOptions * _so)
                             //:
    //current_data(si),user_select(_user_select)
{
  current_data    = si;
  user_select     = _user_select;
  pov             = _pov;
  select          = _select;
  mutex_data      = _mutex;
  store_options   = _so;
  valid_new_frame = false;
}

// ============================================================================
// Destructor                                                                  
LoadingThread::~LoadingThread()
{
}
// ============================================================================
// run()                                                                       
// this function will run in parallel in a thread                              
void LoadingThread::run()
{
  mutex_data->lock(); // lock data access
  // get snapshot component ranges
  ComponentRangeVector * crv = current_data->getSnapshotRange();
  //ComponentRange::list(crv);
  if (crv) {
    assert(crv);
    assert(crv->size());
    user_select->setSelection(select,crv,pov);
    if (current_data->nextFrame(user_select->getIndexesTab(),user_select->getNSel())) {
      valid_new_frame = true;
      store_options->new_frame = true;
    }
  }
  mutex_data->unlock(); // unlock data access
}

}
