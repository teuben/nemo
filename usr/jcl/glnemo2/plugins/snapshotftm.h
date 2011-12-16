// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#ifndef GLNEMOSNAPSHOTFTM_H
#define GLNEMOSNAPSHOTFTM_H
#include <QObject>
#include "snapshotinterface.h"
#include "globaloptions.h"
#include "ftmio.h"
namespace glnemo {


class SnapshotFtm: public SnapshotInterface
{
     Q_OBJECT
     Q_INTERFACES(glnemo::SnapshotInterface)

public:
    SnapshotFtm();

    ~SnapshotFtm();
    
	SnapshotInterface * newObject(const std::string _filename, const int x=0);
    ComponentRangeVector * getSnapshotRange();
    bool isValidData();
    int initLoading(GlobalOptions * so);
    int nextFrame(const int * index_tab, const int nsel);
    int close();
    QString endOfDataMessage();

private:
    int full_nbody;
    bool valid;
    ftm::FtmIO * ftm_io;
};

}

#endif
