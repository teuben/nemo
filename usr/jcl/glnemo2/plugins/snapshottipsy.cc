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
#include <sstream>
#include "snapshottipsy.h"
#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
Q_PLUGIN_METADATA(IID "fr.glnemo2.tipsyPlugin")
#endif

namespace glnemo {

SnapshotTipsy::SnapshotTipsy():SnapshotInterface()
{
  valid=false;
  part_data = NULL;
  part_data = new ParticlesData();
  tipsy_io = NULL;
  interface_type="Tipsy";
  take_gas = take_halo = take_stars = false;
}

// ============================================================================
// Destructor                                                                  
SnapshotTipsy::~SnapshotTipsy()
{
  if (obj)       delete obj;
  if (part_data) delete part_data;
  if (valid) close();
  if (tipsy_io) delete tipsy_io;
}
// ============================================================================
// newObject                                                                   
// instantiate a new object and return a pointer on it                         
SnapshotInterface * SnapshotTipsy::newObject(const std::string _filename, const int x)
{
  if (x) {;} // get rid of compiler warning
  filename = _filename;
  obj      = new SnapshotTipsy();
  std::cerr << "SnapshotTipsy:: " << _filename << "\n";
  obj->setFileName(_filename);

  return obj;
}
// ============================================================================
// isValidData()                                                               
// return true if it's a Tipsy snapshot.
bool SnapshotTipsy::isValidData()
{
  valid=false;
  tipsy_io = new tipsy::TipsyIO(filename);
  int fail = tipsy_io->open(filename);
  if (!fail) {
    valid=true;
  }
  return valid; 
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * SnapshotTipsy::getSnapshotRange()
{
  crv.clear();
  if (valid) {
    crv = tipsy_io->getCRV();
    ComponentRange::list(&crv);
    if (first) {
      first       = false;
      crv_first   = crv;
      nbody_first = tipsy_io->getNtotal();
      time_first  = tipsy_io->getTime();
    }
  }
  return &crv;
}
// ============================================================================
// initLoading()                                                               
int SnapshotTipsy::initLoading(GlobalOptions * so)
{
  load_vel = so->vel_req;
  select_time = so->select_time;
  std::cerr << "select_time ="<<select_time<<"\n";
  go = so;
  
  return 1;
}
// ============================================================================
// nextFrame()                                                                 
int SnapshotTipsy::nextFrame(const int * index_tab, const int nsel)
{
  int status=0;
  stv.clear();
  parseSelectTime();
  load_vel = go->vel_req;

  if (go->select_part=="all" || (go->select_part.find("gas")!=std::string::npos))
    take_gas = true;

  if (go->select_part=="all" || (go->select_part.find("halo")!=std::string::npos))
    take_halo = true;

  if (go->select_part=="all" || (go->select_part.find("stars")!=std::string::npos))
    take_stars = true;



  if (valid && checkRangeTime(tipsy_io->getTime())) {
    status=1;
    if (nsel > *part_data->nbody) {
      //pos
      if (part_data->pos) delete [] part_data->pos;
      part_data->pos = new float[nsel*3];
      //vel
      if (load_vel) {
        if (part_data->vel) delete [] part_data->vel;
        part_data->vel = new float[nsel*3];
      }
      //rho
      if (part_data->rho) delete part_data->rho;
      part_data->rho = new PhysicalData(PhysicalData::rho,nsel);
      for (int i=0; i<nsel; i++) part_data->rho->data[i]=-1.;
      //rneib
      if (part_data->rneib) delete part_data->rneib;
      part_data->rneib = new PhysicalData(PhysicalData::neib,nsel);
      for (int i=0; i<nsel; i++) part_data->rneib->data[i]=-1.;
      //temp
      if (part_data->temp) delete part_data->temp;
      part_data->temp = new PhysicalData(PhysicalData::temperature,nsel);
      for (int i=0; i<nsel; i++) part_data->temp->data[i]=-1.;
      //Ids
      part_data->id.clear();
      for (int i=0; i<nsel; i++) part_data->id.push_back(-1);
    }
    *part_data->nbody = nsel;
    
    std::cerr << "vector size ="<<part_data->id.size() <<"  nsel="<<nsel<<"\n";

    if (tipsy_io->loadData(&part_data->id,part_data->pos,part_data->vel,part_data->rho->data, part_data->rneib->data,part_data->temp->data,
                            index_tab,nsel,load_vel, take_gas, take_halo, take_stars)) {

      part_data->computeVelNorm();
      part_data->rho->computeMinMax();
      part_data->temp->computeMinMax();
      part_data->rneib->computeMinMax();
      if ( ! part_data->rho->isValid()) {
        if (part_data->rho) delete part_data->rho;
        part_data->rho=NULL;
        if (part_data->rneib) delete part_data->rneib;
        part_data->rneib=NULL;
        if (part_data->temp) delete part_data->temp;
        part_data->temp=NULL;
      }
      if (! part_data->timu ) part_data->timu = new float;
      *part_data->timu = tipsy_io->getTime();
      end_of_data=true; // only one frame from an gadget snapshot
    }
    else {
      end_of_data=true;
    }
  }
  return status;
}
// ============================================================================
// close()                                                                     
int SnapshotTipsy::close()
{
  tipsy_io->close();
  end_of_data=false;
  return 1;
}
// ============================================================================
// endendOfDataMessage()                                                       
QString SnapshotTipsy::endOfDataMessage()
{
  QString message=tr("Tipsy Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
}
// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotTipsy

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(tipsyplugin, glnemo::SnapshotTipsy);
#endif
