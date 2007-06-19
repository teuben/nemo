// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "filelistdata.h"
#include "snapshot_data.h"
#include <sstream>
const char * FileListData::magic="#glnemo_file_list";

// ============================================================================
// Constructor                                                                 
FileListData::FileListData(const char * _filename, const char * selp,
			   const char * selt, const bool load_vel_,
                           pthread_mutex_t * _mutex_data):VirtualData(),filename(_filename)
{
  //filename = _filename;
  select_part = selp;
  select_time = selt;
  load_vel    = load_vel_;
  mutex_data  = _mutex_data;
  virtual_data= NULL;
  is_open = FALSE;
  is_parsed = FALSE;
  is_loading_thread = FALSE;
  is_new_data_loaded = FALSE;
  openFile();
}
// ============================================================================
//
FileListData::~FileListData()
{
  if (virtual_data) delete virtual_data;
  close();
}
// ============================================================================
//
void FileListData::openFile()
{
  // open file
  fi.open(filename,std::ios::in);
  if (! fi.is_open()) {
    std::cerr << "Unable to open file ["<<filename<<"] for reading, aborting...\n";
    valid = false;
  }
  else {
    // read magic number
    std::string line;
    getline(fi,line);
    std::string header(magic);
    if (line == header) valid=true;
    else {
      valid=false;
      fi.close();
    }
  } 
}
// ============================================================================
//
int FileListData::close()
{
  int status = 1;
  if (valid) fi.close();
  else status = 0;
  return status;  
}
// ============================================================================
// getLine()
bool FileListData::getLine()
{
  bool status=false,stop=false;
  if (valid) {
    while (!stop && ! fi.eof()) {           // while ! eof
      std::string line;
      getline(fi,line);
      if ( ! fi.eof()) {
        std::istringstream str(line);  // stream line
        std::string parse;
        // following loop parse each lines previously read   
        //
        int cpt=0;
        snaprange="";
        while (  str >> parse   &&              // something to read 
                parse[0] != '#' &&              // not commented out 
                parse[0] != '!'                 // not commented out 
                ) {

            cpt++;
            if (cpt==1) snapname=parse;
            else {
              if (snaprange=="") snaprange=parse;
              else               snaprange+=","+parse;
            }
        }
        if (cpt > 0) {   // something has been read
          if (snaprange=="") snaprange="all";
          stop   = true; // we have a snapshot and a range
          status = true; // so we can stop reading
        } else {
          
        }
      }
      else { // end of file
        stop   = true;
        status = false;
      }
    }
  }
  else status=false;
  return status;
  
}
// ============================================================================
// isValidData()
bool FileListData::isValidData()
{
  return valid;
}
// ============================================================================
// reload()                                                      
// load the snapshot from the beginning                                        
int FileListData::reload(ParticlesSelectVector * psv, const bool load_vel)
{
  int ret=0;
  if (close()) {
    is_open            = FALSE;
    is_parsed          = FALSE;
    is_loading_thread  = FALSE;
    is_new_data_loaded = FALSE;
    is_end_of_data     = FALSE;
    openFile();
    ret = loadPos(psv, load_vel);
  }
  else {
    ret = -2;
  }
  return ret;
}
// ============================================================================
// getData()
//                                            
void FileListData::getData(const ParticlesData * _p_data,
                    ParticlesSelectVector  * psv)
{
  emit loadedData(_p_data,psv);
}
// ============================================================================
// loadPos()
// Read snapshot positions, send them to the glbox                             
int FileListData::loadPos(ParticlesSelectVector * psv, const bool load_vel)
{
  bool stop=false,status=false;

  while (!stop && getLine() ) {
    //std::cerr << "SNAPNAME ["<<snapname<<"] range <"<<snaprange<<">\n";
    // try to instantiate a snapshotdata class
    VirtualData * new_virtual_data = new SnapshotData(snapname.c_str(),
                                    snaprange.c_str(),select_time,load_vel,mutex_data);
    if (! new_virtual_data->isValidData()) {
      std::cerr << "File [" << snapname << "] is not a NEMO snapshot, aborting...\n";
      delete new_virtual_data;
      //exit(1);
    }
    else {
      pthread_mutex_lock(mutex_data);
      if (virtual_data) delete(virtual_data);
      virtual_data = new_virtual_data;
      pthread_mutex_unlock(mutex_data);
      virtual_data->is_loading_thread = is_loading_thread;
      // establish Signal connexion
      connect(virtual_data,SIGNAL(loadedData(const ParticlesData *,
      ParticlesSelectVector * )),
      this,SLOT(getData(const ParticlesData *, ParticlesSelectVector *)));
            
      psv->clear();   // clear particles range vectors
      VirtualParticlesSelect::nb_select = 0;
      
      if (virtual_data->loadPos(psv, load_vel)) {
        part_data = virtual_data->getParticlesData();
        status = true;
        stop   = true;
      } else {
        stop   = false;
      }
    }
    
  }
 if (!status) is_end_of_data = TRUE;
 return status;
  
}
// ============================================================================
// uploadGlData()
// Upload Data to GLbox                                                        
void FileListData::uploadGlData(ParticlesSelectVector * psv)
{
  virtual_data->uploadGlData( psv);
}
// ============================================================================
// endOfDataMessage()
// return end of data message                                                  
QString FileListData::endOfDataMessage()
{
  QString fname=filename;
  QString message="File ["+fname+"] no more files in the list";
  return message;
}
// ============================================================================
