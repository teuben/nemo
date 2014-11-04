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
#ifndef GLNEMOSNAPSHOTLIST_H
#define GLNEMOSNAPSHOTLIST_H
#include <QObject>
#include <fstream>
#include "snapshotinterface.h"
#include "pluginsmanage.h"
#include "globaloptions.h"

namespace glnemo {


class SnapshotList: public SnapshotInterface
{
  Q_OBJECT
  Q_INTERFACES(glnemo::SnapshotInterface)
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
  Q_PLUGIN_METADATA(IID "fr.glnemo2.listPlugin" FILE "listPlugin.json")
#endif
     
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
    int getNumberFrames() { return vector_file.size();}
    int getCurrentFrameIndex() { return current_file_index;} // >
    virtual void checkJumpFrame(const int _v=-1) {
      frame.lock();
      jump_frame = _v;
      if (jump_frame>=0 && jump_frame < (int) vector_file.size()) {
        end_of_data=false;
      }
      frame.unlock();
    }
    bool isEndOfData() {
      if ( (current_file_index>0 && current_file_index<((int) vector_file.size()-1) ) ||  // in the limit
           ( play_forward &&  current_file_index<(int) vector_file.size()-1)          ||  // forward
           (!play_forward &&  current_file_index>0))                                       // backward
        end_of_data=false;
      else {
        end_of_data=true;
        if (current_file_index>=(int)vector_file.size()) {
          current_file_index=vector_file.size()-1;
        }
      }
      return end_of_data;
    }
    std::vector<std::string> getVectorFile() { return vector_file; }

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
    std::vector<std::string> vector_file;
    int current_file_index;
    bool getNextFile();
};

}

#endif
