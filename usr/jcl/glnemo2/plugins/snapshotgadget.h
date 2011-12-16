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
#ifndef GLNEMOSNAPSHOTGADGET_H
#define GLNEMOSNAPSHOTGADGET_H
#include <QObject>
#include "snapshotinterface.h"
#include "gadgetio.h"
#include "globaloptions.h"
namespace glnemo {


class SnapshotGadget: public SnapshotInterface{
     Q_OBJECT
     Q_INTERFACES(glnemo::SnapshotInterface)
public:
    SnapshotGadget();

    ~SnapshotGadget();
	SnapshotInterface * newObject(const std::string _filename, const int x=0);
    bool isValidData();
    ComponentRangeVector * getSnapshotRange();
    int initLoading(GlobalOptions * so);
    int nextFrame(const int * index_tab, const int nsel);
           // index_tab = array of selected indexes (size max=part_data->nbody)
           // nsel      = #particles selected in index_tab                     
           // particles not selected must have the value '-1'                  
    int close();
    QString endOfDataMessage();

private:
    bool valid;
    gadget::GadgetIO * gadget_io;
};

}

#endif
