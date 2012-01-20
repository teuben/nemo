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
#ifndef GLNEMOPLUGINSMANAGE_H
#define GLNEMOPLUGINSMANAGE_H
#include <QObject>
#include <string>

namespace glnemo {

// IMPORTANT : it's mandatory to prototype the class which belong to
// the namespace INSIDE the namespace....                           
class SnapshotInterface;

class PluginsManage: public QObject
{
  Q_OBJECT
public:
    PluginsManage();
    ~PluginsManage();

    SnapshotInterface * getObject(const std::string filename);
};

}

#endif
