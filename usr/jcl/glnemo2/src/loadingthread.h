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
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
 */
#ifndef GLNEMOLOADINGTHREAD_H
#define GLNEMOLOADINGTHREAD_H
#include "snapshotinterface.h"
#include "userselection.h"
#include "particlesobject.h"
#include "globaloptions.h"

#include <QThread>
#include <QMutex>
namespace glnemo {


class LoadingThread : public QThread {
  Q_OBJECT
  public:
    LoadingThread(SnapshotInterface *,
                  UserSelection *,
                  ParticlesObjectVector *, const std::string,
                  QMutex *, GlobalOptions *);
    void run();
    ~LoadingThread();
    bool isValidLoading() { return valid_new_frame;};
  private:
    SnapshotInterface * current_data;
    UserSelection * user_select;
    GlobalOptions * store_options;
    bool valid_new_frame;
    ParticlesObjectVector * pov;
    std::string select;
    QMutex * mutex_data;
};

}

#endif
