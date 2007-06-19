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
// AcquireDataThread class definition                                          
//                                                                             
// Allow parallele data acquisition via a THREAD                               
// ============================================================================
#ifndef ACQUIRE_DATA_THREAD_H
#define ACQUIRE_DATA_THREAD_H
#include <qobject.h>
#include <qthread.h> 
#include "virtual_data.h"

class AcquireDataThread : public QThread
{
  //Q_OBJECT
  public:
  AcquireDataThread(VirtualData *, ParticlesSelectVector *, const bool );
  //AcquireDataThread();
  ~AcquireDataThread();
  
  void run();
  bool is_loaded;  
  private:
  VirtualData * virtual_data;
  ParticlesSelectVector * psv;
  bool load_vel;

};  // DO NOT FORGET ';'
#endif
//
