// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//   <        13388 Marseille cedex 13 France                                   
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

SnapshotPhiGrape::SnapshotPhiGrape():SnapshotInterface()
{ 
  valid=false;
  part_data = NULL;
  part_data = new ParticlesData();
  interface_type    ="PhiGRAPE";
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
    valid = detectHeader();
  }
 return valid; 
}
// ============================================================================
// detectHeader()                                                              
// returne true if it's a valid phiBar snapshot                                
bool SnapshotPhiGrape::detectHeader()
{
  bool ok=false;
  char buff[KB];
  std::stringstream ss;
  gzgets(file, buff, KB);          //  frame number
  ss.str(buff);
  ss >> frame_number;
  
  gzgets(file, buff, KB);          //  particle number
  ss.str(buff);
  ss >> full_nbody;
  
  gzgets(file, buff, KB);          //  snapshot time
  ss.str(buff);
  ss >> *part_data->timu;
  
  int n1,n2;
  z_off_t fpos;
  fpos = gztell(file); // save the current position
  
  gzgets(file, buff, KB);          // particle record #0
  ss.str(buff);
  ss >> n1;
  
  gzgets(file, buff, KB);          // particle record #1
  ss.str(buff);
  ss >> n2;

  if (n2 == n1+1) {
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
int SnapshotPhiGrape::initLoading(const bool _load_vel, const std::string _select_time)
{
  load_vel = _load_vel;
  select_part="all";
  select_time=_select_time;
  std::cerr << "SnapshotPhiGrape::initLoading select_time = " << select_time << "\n";
  return 1;
}
// ============================================================================
// nextFrame()                                                                 
int SnapshotPhiGrape::nextFrame(const int * index_tab, const int nsel)
{
  if (valid) {
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
      }
      cpt++;
    }
    
  }
  end_of_data = true;
}
// ============================================================================
// close()                                                                     
int SnapshotPhiGrape::close()
{
  if (valid) {
    gzclose(file);
    end_of_data=false;
    valid = false;
    return 1; 
  }
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
