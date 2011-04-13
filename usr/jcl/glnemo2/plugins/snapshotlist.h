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
#ifndef GLNEMOSNAPSHOTLIST_H
#define GLNEMOSNAPSHOTLIST_H
#include <QObject>
#include <fstream>
#include "snapshotinterface.h"
#include "pluginsmanage.h"
#include "globaloptions.h"

namespace glnemo {


class SnapshotList: public QObject,
                     public SnapshotInterface
{
     Q_OBJECT
     Q_INTERFACES(glnemo::SnapshotInterface);
     
public:
    SnapshotList();
    ~SnapshotList();

	SnapshotInterface * newObject(const std::string _filename, const int x=0);
    ComponentRangeVector * getSnapshotRange();
    int initLoading(GlobalOptions * so);
    bool isValidData();
    int nextFrame(const int * index_tab, const int nsel);
    int close();
    QString endOfDataMessage();
private:
    std::ifstream fi;
    static const char * magic;
    bool valid;
    bool openFile();
    bool getLine(const bool force=false);
    std::string snapshot;
    SnapshotInterface * current_data;
    PluginsManage * plugins;
    std::string interface_type_ori;
    GlobalOptions * go;
    QString dirpath;
};

}

#endif
