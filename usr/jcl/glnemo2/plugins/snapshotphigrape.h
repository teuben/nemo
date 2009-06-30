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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
*/
#ifndef SNAPSHOTPHIBAR_H
#define SNAPSHOTPHIBAR_H
#include <QObject>
#include "snapshotinterface.h"
#include "zlib.h"

namespace glnemo {

class SnapshotPhiGrape: public QObject,
                      public SnapshotInterface
{
    Q_OBJECT
    Q_INTERFACES(glnemo::SnapshotInterface);

public:
    SnapshotPhiGrape();

    ~SnapshotPhiGrape();
    
    SnapshotInterface * newObject(const std::string _filename);
    ComponentRangeVector * getSnapshotRange();
    bool isValidData();
    int initLoading(const bool _load_vel, const std::string _select_time);
    int nextFrame(const int * index_tab, const int nsel);
    int close();
    QString endOfDataMessage();
 private:
    bool valid;
    gzFile           file;               // file handle for compressed file
    bool detectHeader();
    int full_nbody, frame_number;
};

}
#endif // SNAPSHOTPHIBAR_H
