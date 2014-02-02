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
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#ifndef GLNEMOSNAPSHOTRAMSES_H
#define GLNEMOSNAPSHOTRAMSES_H
#include <QObject>
#include "snapshotinterface.h"
#include "camr.h"
#include "cpart.h"
#include "globaloptions.h"

namespace glnemo {

class SnapshotRamses: public SnapshotInterface
{
  Q_OBJECT
  Q_INTERFACES(glnemo::SnapshotInterface)
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
  Q_PLUGIN_METADATA(IID "fr.glnemo2.ramsesPlugin" FILE "ramsesPlugin.json")
#endif

public:
    SnapshotRamses();

    ~SnapshotRamses();
    
    SnapshotInterface * newObject(const std::string _filename, const int x=0);
    ComponentRangeVector * getSnapshotRange();
    bool isValidData();
    int initLoading(GlobalOptions * so);
    int nextFrame(const int * index_tab, const int nsel);
    int close() { return 1;}
    QString endOfDataMessage();

private:
    ramses::CAmr * amr;
    ramses::CPart * part;
    int namr;
    int nstars;
    int ndm;
    bool valid;
    int full_nbody;
    GlobalOptions * go;
    bool take_gas, take_halo, take_stars;
};

}

#endif
