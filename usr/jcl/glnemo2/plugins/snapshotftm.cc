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
#include <QtGui> // Mandatory for plugins management
#include "snapshotftm.h"
#include "particlesdata.h"
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
Q_PLUGIN_METADATA(IID "fr.glnemo2.ftmPlugin")
#endif

namespace glnemo {

using namespace ftm;

// ============================================================================
// Constructor                                                                 
SnapshotFtm::SnapshotFtm():SnapshotInterface()
{
  valid=false;
  part_data = new ParticlesData();
  ftm_io = NULL;
  interface_type="Ftm";
}
// ============================================================================
// Destructor                                                                  
SnapshotFtm::~SnapshotFtm()
{
  if (obj)       delete obj;
  if (part_data) delete part_data;
  if (valid) close();
  if (ftm_io)    delete ftm_io;
}
// ============================================================================
// newObject                                                                   
// instantiate a new object and return a pointer on it                         
SnapshotInterface * SnapshotFtm::newObject(const std::string _filename, int x)
{
  if (x) {;} // get rid of compiler warning
  filename = _filename;
  obj      = new SnapshotFtm();
  obj->setFileName(_filename);

  return obj;
}
// ============================================================================
// isValidData()                                                               
// return true if it's a Ftm snapshot.                                         
bool SnapshotFtm::isValidData()
{
  valid=false;
  if (isFileExist()) {
    ftm_io = new FtmIO(filename);
    int fail = ftm_io->open();
    if (!fail) valid=true;
  }
  return valid; 
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * SnapshotFtm::getSnapshotRange()
{
  crv.clear();
  if (valid) {
    ComponentRange * cr = new ComponentRange();
    cr->setData(0,ftm_io->getNtotal()-1);
    cr->setType("all");
    crv.push_back(*cr);
    const FtmComponent * comp = ftm_io->getComp();
    //std::cerr << "Halo.type = " << comp->halo.type << "\n";
    if ( comp->halo.type == "halo") crv.push_back(comp->halo);
    if ( comp->disk.type == "disk") crv.push_back(comp->disk);
    if ( comp->gas.type  == "gas" ) crv.push_back(comp->gas);
    //ComponentRange::list(&crv);
    if (first) {
      first       = false;
      crv_first   = crv;
      nbody_first = ftm_io->getNtotal();
      time_first  = ftm_io->getTime();
    }
    delete cr;
  }
  return &crv;
}
// ============================================================================
// initLoading()                                                               
int SnapshotFtm::initLoading(GlobalOptions * so)
{
  go = so;
  load_vel = so->vel_req;
  select_time = so->select_time;
  return 1;
}
// ============================================================================
// nextFrame()                                                                 
int SnapshotFtm::nextFrame(const int * index_tab, const int nsel)
{
  int status=0;
  load_vel = go->vel_req;
  if (valid) {
  if (ftm_io->read(part_data,index_tab,nsel,load_vel)) {
    status=1;
    part_data->computeVelNorm();
    if (part_data->rho) {
      part_data->rho->computeMinMax();
    }
    end_of_data=true; // only one frame from an ftm snapshot
  }
    else {
      end_of_data=true;
    }
  }
  return status;
}
// ============================================================================
// close()                                                                     
int SnapshotFtm::close()
{
  ftm_io->close();
  end_of_data=false;
  return 1;
}
// ============================================================================
// endendOfDataMessage()                                                       
QString SnapshotFtm::endOfDataMessage()
{
  QString message=tr("Ftm Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
}
// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotFtm                             

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(ftmplugin, glnemo::SnapshotFtm);
#endif
