// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
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
#include "snapshotphigrape.h"
#include "zlib.h"
#include <cstdio>
#include <sstream>

namespace glnemo {
#define KB        512
#define N_MAX     128*KB
// ============================================================================
// constructor                                                                 
SnapshotPhiGrape::SnapshotPhiGrape():SnapshotInterface()
{ 
  valid=false;
  part_data = NULL;
  part_data = new ParticlesData();
  interface_type    ="PhiGRAPE";
  full_nbody = 0;
}

// ============================================================================
// destructor                                                                  
SnapshotPhiGrape::~SnapshotPhiGrape()
{
  std::cerr << "SnapshotPhiGrape::~SnapshotPhiGrape()\n";
  if (obj)       delete obj;
  if (part_data) delete part_data;
  
  if (pos)   free ((float *) pos);
  if (vel)   free ((float *) vel);

  if (valid) close();
}
// ============================================================================
// newObject                                                                   
// instantiate a new object and return a pointer on it                         
SnapshotInterface * SnapshotPhiGrape::newObject(const std::string _filename)
{
  filename = _filename;
  obj      = new SnapshotPhiGrape();
  obj->setFileName(_filename);
  return obj;
}
// ============================================================================
// isValidData()                                                               
// return true if it's a SnapshotPhiGrape snapshot. Standard input '-' is assumed
// to be a valid SnapshotPhiGrape snapshot.                                      
bool SnapshotPhiGrape::isValidData()
{
  valid=false;
  file = gzopen(filename.c_str(),"r");
  if (file) {
    try {
      valid = detectHeader();
    }
    catch (int e) {
      std::cerr << "SnapshotPhiGrape::detectHeader failed...\n";
    }
  }
 return valid; 
}
// ============================================================================
// detectHeader()                                                              
// return  true if it's a valid phiBar snapshot                                
bool SnapshotPhiGrape::detectHeader()
{
  bool ok=false;
  char buff[KB];
  std::stringstream ss;
  if (gzgets(file, buff, KB)==Z_NULL)          //  frame number
    throw 10;
  ss.str(buff);
  ss >> frame_number;
  
  if (gzgets(file, buff, KB)==Z_NULL)          //  particle number
    throw 10;
  ss.str(buff);
  ss >> full_nbody;
  
  if (gzgets(file, buff, KB)==Z_NULL)          //  snapshot time
    throw 10;
  ss.str(buff);
  ss >> *part_data->timu;
  
  int n1,n2;
  z_off_t fpos;
  fpos = gztell(file); // save the current position
  
  if (gzgets(file, buff, KB)==Z_NULL)          // particle record #0
    throw 10;
  ss.str(buff);
  ss >> n1;
  
  if (gzgets(file, buff, KB)==Z_NULL)          // particle record #1
    throw 10;
  ss.str(buff);
  ss >> n2;

  if (n2 == n1+1 && full_nbody>0 && n1 >=0 && n2 >0) { // BINGO !! , it's a valid phiGRAPE file
    ok = true;
    gzseek(file, fpos, SEEK_SET); // go back to the first record
  }
  return ok;
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * SnapshotPhiGrape::getSnapshotRange()
{
  crv.clear();
  if (valid) {
    ComponentRange * cr = new ComponentRange();
    cr->setData(0,full_nbody-1);
    cr->setType("all");
    crv.push_back(*cr);
    //ComponentRange::list(&crv);
    delete cr;
    if (first) {
      first       = false;
      crv_first   = crv;
      nbody_first = full_nbody;
      //time_first  = 0.0;
    }
  }
  return &crv;
}
// ============================================================================
// initLoading()                                                               
int SnapshotPhiGrape::initLoading(GlobalOptions * so)
{
  load_vel = so->vel_req;
  select_part="all";
  select_time=so->select_time;
  std::cerr << "SnapshotPhiGrape::initLoading select_time = " << select_time << "\n";
  return 1;
}
// ============================================================================
// nextFrame()                                                                 
int SnapshotPhiGrape::nextFrame(const int * index_tab, const int nsel)
{
  int status=0;
  if (valid) {
    status=1;
    if (nsel > *part_data->nbody) {
      if (part_data->pos) delete [] part_data->pos;
      part_data->pos = new float[nsel*3];
      if (load_vel) {
        if (part_data->vel) delete [] part_data->vel;
        part_data->vel = new float[nsel*3];
      }
    }
    assert(nsel<=full_nbody);
    *part_data->nbody=nsel;
    int cpt=0;
    char buff[KB];
    std::stringstream ss;
    bool first=true;
    for (int i=0; i<full_nbody; i++) {
  
      gzgets(file, buff, KB);   // read a line from the file
      int idx=index_tab[i];
      if (idx!=-1) { // it's a valid particle
        int nn; 
        float v[3],dummy;
        ss.str(buff);
        // particle index + mass
        ss >> nn >> dummy; 
        // x y z
        ss >> part_data->pos[cpt*3+0] >> part_data->pos[cpt*3+1] >> part_data->pos[cpt*3+2];
        // vx vy vz
        ss >> v[0] >> v[1] >> v[2];
        if (load_vel) 
          for (int j=0; j<3; j++)
            part_data->vel[cpt*3+j] = *(v+j);
        float rho1,hsml;
        
        ss >> rho1 >> hsml;    // try to read RHO and HSML
        
        if (first && !ss.eof()) { // first time and rho+hsml exist
          first = false;
          if (nsel > *part_data->nbody|| !part_data->rho) {
            if (part_data->rho) delete part_data->rho;
            part_data->rho = new PhysicalData(PhysicalData::rho,nsel);
          }
          if (nsel > *part_data->nbody || !part_data->rneib) {
            if (part_data->rneib) delete part_data->rneib;
            part_data->rneib = new PhysicalData(PhysicalData::neib,nsel);
          }
        }
        if (!ss.eof()) {
          part_data->rho->data[cpt]   = rho1;
          part_data->rneib->data[cpt] = hsml;
        }
        ss.clear(); // clear error state flag like eof()
        cpt++;
      }
      
    }
    if (part_data->rho) part_data->rho->computeMinMax();
  }
  end_of_data = true;
  return status;
}
// ============================================================================
// close()                                                                     
int SnapshotPhiGrape::close()
{
  int ret=0;
  if (valid) {
    gzclose(file);
    end_of_data=false;
    valid = false;
    ret=1; 
  }
  return ret;
}
// ============================================================================
// endendOfDataMessage()                                                       
QString SnapshotPhiGrape::endOfDataMessage()
{
  QString message=tr("PhiGrape Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
}

// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotNemo                            

Q_EXPORT_PLUGIN2(phigrapeplugin, glnemo::SnapshotPhiGrape);
