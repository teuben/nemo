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
// AcquireDataThread class implementation                                      
//                                                                             
// Allow parallele data acquisition via a THREAD                               
// ============================================================================

//#include <qmessagebox.h>
#include "acquire_data_thread.h"

// ============================================================================
// Constructor
AcquireDataThread::AcquireDataThread(VirtualData * vd, 
                                     ParticlesRangeVector * pr_v)
{
  virtual_data = vd;
  prv          = pr_v;
  is_loaded    = false;
}
// ============================================================================
// destructor
AcquireDataThread::~AcquireDataThread()
{
}
// ============================================================================
// run method
void AcquireDataThread::run()
{
#if 1 
  if (! virtual_data->loadPos(prv) ) {
    QString message="End of snapshot Reached !";
    //QMessageBox::information( this,"Warning",message,"Ok");
  }
#else
  int i=0;
  while (1) {
    i++;
    std::cerr << "I am in the running thread [" << i << "]\n";
    sleep(1);
  }
#endif
  virtual_data->is_loading_thread = FALSE;
  is_loaded = true;
}
//
