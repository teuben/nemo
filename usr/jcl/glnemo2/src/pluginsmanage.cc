// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
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
#include <iostream>
#include <QtGui>
#include "pluginsmanage.h"
#include "snapshotinterface.h"

namespace glnemo {

PluginsManage::PluginsManage()
{
}


PluginsManage::~PluginsManage()
{
}

SnapshotInterface * PluginsManage::getObject(const std::string filename )
{
  // first: try to load STATIC plugins
  foreach (QObject *plugin, QPluginLoader::staticInstances()) {
    //populateMenus(plugin);
    SnapshotInterface * iface = qobject_cast<SnapshotInterface *>(plugin);
    if (iface) {
      SnapshotInterface * iface1 = iface->newObject(filename);
      std::cerr << "Trying Interface:" << iface1->getInterfaceType() << "\n";
      if (iface1->isValidData()) {        
        return iface1;
      }
      else {
         delete iface1;
      }
    }
  }
  
#if 1
  // second: try to load DYNAMIC plugins
  QDir pluginsDir= QDir(qApp->applicationDirPath());
  pluginsDir.cd("../lib");

  foreach (QString fileName, pluginsDir.entryList(QDir::Files)) {
    QPluginLoader loader(pluginsDir.absoluteFilePath(fileName));
    QObject *plugin = loader.instance();
    if (plugin) {
      SnapshotInterface * iface = qobject_cast<SnapshotInterface *>(plugin);
      if (iface) {
        SnapshotInterface * iface1 = iface->newObject(filename);
        return iface1;
      }
    }
  }
#endif
  return NULL;
}

}
