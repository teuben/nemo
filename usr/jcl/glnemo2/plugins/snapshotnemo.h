// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOSNAPSHOTNEMO_H
#define GLNEMOSNAPSHOTNEMO_H
#include <QObject>
#include "snapshotinterface.h"
#include "globaloptions.h"

namespace glnemo {


class SnapshotNemo : public SnapshotInterface
{
     Q_OBJECT
     Q_INTERFACES(glnemo::SnapshotInterface)

public:
    SnapshotNemo();
    ~SnapshotNemo();
	SnapshotInterface * newObject(const std::string _filename, const int x=0);
    ComponentRangeVector * getSnapshotRange();
    int initLoading(GlobalOptions * so);
    bool isValidData();
    int nextFrame(const int * index_tab, const int nsel);
    int close();
    QString endOfDataMessage();
private:
    int full_nbody;
    bool valid;
    float * npos, *nvel, *nrho, *nrneib, *ntimu;
    int * nnemobits , *nnbody, *nid;
    bool first_stream;
    int status_ionemo;
};

}

#endif
