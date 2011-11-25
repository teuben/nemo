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
#include <QtGui> // Mandatory for plugins management
#include <QDir>
#include <sstream>
#include "snapshotlist.h"

namespace glnemo {
const char * SnapshotList::magic="#glnemo_file_list";
// ============================================================================
SnapshotList::SnapshotList():SnapshotInterface()
{
  valid=false;
  part_data = NULL;
  part_data = new ParticlesData();
  current_data = NULL;
  plugins = new PluginsManage();
  interface_type    ="List";
  interface_type_ori="List";
}

// ============================================================================
SnapshotList::~SnapshotList()
{
  close();
  if (current_data)
    delete current_data;
  if (part_data)
    delete part_data;
}
// ============================================================================
SnapshotInterface * SnapshotList::newObject(const std::string _filename, const int x)
{
  if (x) {;} // get rid of compiler warning
  filename = _filename;
  obj      = new SnapshotList();
  obj->setFileName(_filename);
  return obj;
}
// ============================================================================
ComponentRangeVector * SnapshotList::getSnapshotRange()
{
  bool stop=false;
  ComponentRangeVector * crv=NULL;
  
  int ipvs = part_data->getIpvs(); // backup ipvsr
  while (!stop && getLine()) {  
    if (current_data) {
      delete current_data;
    }
    current_data = plugins->getObject(snapshot);
    if (current_data) {
      stop=true;
      std::cerr << "SnapshotList::getSnapshotRange select_time=" << select_time << "\n";
      current_data->initLoading(go);      
      //loadNewData(select,select_time,load_vel);
      // get snapshot component ranges
      crv = current_data->getSnapshotRange();
      interface_type = interface_type_ori+" of "+current_data->getInterfaceType();
      //ComponentRange::list(crv);
      assert(crv);
      assert(crv->size());
      //user_select->setSelection(select,crv,&pov);
      current_data->part_data->setIpvs(ipvs); // copy ipvs from the previous shots
      // load from disk
      //current_data->nextFrame(user_select->getIndexesTab(),user_select->getNSel());
      if (first) {
	first       = false;
	crv_first   = current_data->crv_first;
	nbody_first = current_data->nbody_first;
	time_first  = current_data->time_first;
      }
    }
  }
  if (!stop) end_of_data = true; // EOF reached
  return crv;
}
// ============================================================================
int SnapshotList::initLoading(GlobalOptions * so)
{
  go = so;
  load_vel    = so->vel_req;
  select_time = so->select_time;
  std::cerr << "SnapshotList::initLoading select_time=" << select_time << "\n";
  isValidData();  // call it in case of reload
  return 1;
}
// ============================================================================
bool SnapshotList::isValidData()
{
  std::cerr << "current dir = " << ((QDir::current()).dirName()).toStdString() << "\n";
  QDir dir(QString(filename.c_str()));
  dirpath=QString(filename.c_str());          // absolute file name
  dirpath.remove(dir.dirName());              // dirname
  dir = QDir(dirpath);                        // dirpath
  if (dir.exists()) {                         // path exist
    bool b=dir.cd(dirpath);                   // cd to path
    std::cerr << "bool b = " << b << "\n";
    std::cerr << "current dir = " << ((QDir::current()).dirName()).toStdString() << "\n";
  }

  if (!valid)
    valid=openFile();

  return valid;
  
}
// ============================================================================
bool SnapshotList::openFile()
{
  bool status;
    // open file
  if (filename == "-")
    ;//fi = &std::cin;//fi.open(std::cin,std::ios::in);
  else
    fi.open(filename.c_str(),std::ios::in);
  if (! fi.is_open()) {
    std::cerr << "Unable to open file ["<<filename<<"] for reading, aborting...\n";
    status = false;
  }
  else {
    // read magic number
    std::string line;
    getline(fi,line);
    std::string header(magic);
    if (line == header) status=true;
    else {
#if 1
      fi.seekg(0, std::ios::beg); // go back to the beginning
      if (getLine(true)) { // read the first file
        SnapshotInterface * test_data = plugins->getObject(snapshot);
        if (test_data) { // it's a valid snaphot
          delete test_data;
          status = true;
          fi.seekg(0, std::ios::beg); // go back to the beginning
        }
      } 
      else {        
        status=false;
        fi.close();
      }
#else
      status=false;
        fi.close();
#endif
    }
  }
  return status;
}
// ============================================================================
int SnapshotList::nextFrame(const int * index_tab, const int nsel)
{
  int status=0;
  if (valid) {
    status = current_data->nextFrame(index_tab,nsel);
    if (status) {
      *part_data = *current_data->part_data;
    }
  }
  return status;
}
// ============================================================================
bool SnapshotList::getLine(const bool force)
{
  bool status=false,stop=false;
  if (valid || force) {
    while (!stop && ! fi.eof()) {           // while ! eof
      std::string line;
      getline(fi,line);
      if ( ! fi.eof()) {
        std::istringstream str(line);  // stream line
        std::string parse;
        // following loop parse each lines previously read   
        //
        int cpt=0;
        while (  str >> parse   &&              // something to read 
                parse[0] != '#' &&              // not commented out 
                parse[0] != '!' &&              // not commented out 
                parse[0] != '\n'                // not a blank line
                ) {
            cpt++;
            if (cpt==1) snapshot=parse;
        }
        if (cpt > 0 ) {
	  unsigned int i=0;
	  while(i<snapshot.length() && snapshot[i]==' ') i++; // search first non blank
	  if (i<snapshot.length() && snapshot[i]!='/')        // first char not a '/'  
	    snapshot = dirpath.toStdString() + snapshot;      // append to dirpath     
	  stop   = true; // we have a snapshot
	  status = true; // so we can stop reading          
        }
      }
      else { // end of file
        stop   = true;
        status = false;
      }
      std::cerr << "SnapshotList::getLine line="<<line<<"\n";
    }
  }
  else status=false;
  return status;

}
// ============================================================================
int SnapshotList::close()
{
  int status=1;
  if (valid) {
    fi.close();
    current_data->close();
    valid=false;
    end_of_data=false;
  }
  else status=0;
  return status;
}
// ============================================================================
QString SnapshotList::endOfDataMessage()
{
  QString message=tr("Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
} // end of namespace glnemo

// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotList                            

Q_EXPORT_PLUGIN2(listplugin, glnemo::SnapshotList);
//
