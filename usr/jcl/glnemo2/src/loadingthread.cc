// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "loadingthread.h"
#include <assert.h>

namespace glnemo {

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


LoadingThread::~LoadingThread()
{
}

void LoadingThread::run()
{
  // get snapshot component ranges
  ComponentRangeVector * crv = current_data->getSnapshotRange();
  //ComponentRange::list(crv);
  if (crv) {
    assert(crv);
    assert(crv->size());
    mutex_data->lock();
#if 0
    if (current_data->getInterfaceType() == "SnapshotList")
#endif
      user_select->setSelection(select,crv,pov);
    if (current_data->nextFrame(user_select->getIndexesTab(),user_select->getNSel())) {
      valid_new_frame = true;
      store_options->new_frame = true;
    }
    mutex_data->unlock();
  }
}

}
