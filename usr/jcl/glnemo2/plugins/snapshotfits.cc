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
#include "snapshotfits.h"
#include <cmath>
#include <limits>

#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
Q_PLUGIN_METADATA(IID "fr.glnemo2.fitsPlugin")
#endif

namespace glnemo {
using namespace CCfits;

SnapshotFits::SnapshotFits():SnapshotInterface()
{
  valid=false;
  part_data = NULL;
  part_data = new ParticlesData();
  image = NULL;
  pInfile = NULL;
  interface_type="Fits";
}

// ============================================================================
// Destructor
SnapshotFits::~SnapshotFits()
{
  if (obj)       delete obj;
  if (part_data) delete part_data;
  if (valid) close();
  if (valid) {
      delete pInfile;

  }
  //if (fits_io) delete fits_io;
}
// ============================================================================
// newObject
// instantiate a new object and return a pointer on it
SnapshotInterface * SnapshotFits::newObject(const std::string _filename, const int x)
{
  if (x) {;} // get rid of compiler warning
  filename = _filename;
  obj      = new SnapshotFits();
  std::cerr << "SnapshotFits:: " << _filename << "\n";
  obj->setFileName(_filename);

  return obj;
}
// ============================================================================
// isValidData()
// return true if it's a Fits snapshot.
bool SnapshotFits::isValidData()
{
    valid=false;
    try {
        pInfile = new FITS(filename,Read,false);
        try {
            image =  &pInfile->pHDU();
            std::cerr << "There is a pHDU\n";
            std::cerr << "#axis ="<<image->axes() << std::endl;
            if (image->axes() < 2) {
                throw CCfits::FITS::OperationNotSupported("",true);
            } else {
                //(pInfile->pHDU()).read(contents);
                ((PHDU *) image)->read(contents);
            }

        } catch (FitsException& e) { // it fails to read pHDU
            image  = &pInfile->extension(1);
            std::cerr << "There is an extension\n";
            std::cerr << "#axis ="<<image->axes() << std::endl;
            //(pInfile->extension(1)).read(contents);
            ((ExtHDU *) image)->read(contents);
        }
        valid = true;
    } catch (FitsException&) {
        std::cerr << " Fits Exception Thrown....\n";
        valid = false;
    }
    return valid;
}
// ============================================================================
// getSnapshotRange
ComponentRangeVector * SnapshotFits::getSnapshotRange()
{

  crv.clear();
  if (valid) {
      // read all user-specifed, coordinate, and checksum keys in the image
      image->readAllKeys();
      std::cerr << "#axis ="<<image->axes() << std::endl;
      std::cerr << "NBODY = " << contents.size() << "\n";
      // MUST BE SET from CLI or GUI
      float dmin = std::numeric_limits<float>::min();
      float dmax = std::numeric_limits<float>::max();
      if (go->phys_min_glob!=-1) {
          //dmin = go->phys_min_glob;
      }
      if (go->phys_max_glob!=-1) {
          //dmax = go->phys_max_glob;
      }
      int   zmin = std::numeric_limits<int>::min();
      int   zmax = std::numeric_limits<int>::max();

      long ax(image->axis(0));
      long ay(image->axis(1));
      long az;
      if (image->axes()>=3) {
          az=(image->axis(2));
      } else {
          az=0;
      }
      // Load data and fill up arrays
      int valid_nbody=0;
      for (unsigned long i=0; i<contents.size(); i++) {
          if (std::isfinite(contents[i]) && contents[i]>=dmin && contents[i]<=dmax) {
              long z_i=int(i/(ax*ay)); // current Z plane
              assert(z_i<=az);
              if (z_i>=zmin && z_i<=zmax) { // inside Z selection
                  valid_nbody++;
              }
          }
      }

      glnemo::ComponentRange cr;
      // all
      cr.setData(0,valid_nbody-1);
      cr.setType("all");
      crv.clear();
      crv.push_back(cr);
      ComponentRange::list(&crv);
      if (first) {
          first       = false;
          crv_first   = crv;
          time_first = 0.0;
          nbody_first = valid_nbody;
      }
  }
  return &crv;
}
// ============================================================================
// initLoading()
int SnapshotFits::initLoading(GlobalOptions * so)
{
  load_vel = so->vel_req;
  select_time = so->select_time;
  std::cerr << "select_time ="<<select_time<<"\n";
  go = so;
  return 1;
}
// ============================================================================
// nextFrame()
int SnapshotFits::nextFrame(const int * index_tab, const int nsel)
{
  int status=0;
  stv.clear();
  parseSelectTime();
  load_vel = go->vel_req;


  //if (valid && checkRangeTime(fits_io->getTime())) {
  if (valid && checkRangeTime(0.0)) {
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

      //Ids
      part_data->id.clear();
      for (int i=0; i<nsel; i++) part_data->id.push_back(-1);
    }
    *part_data->nbody = nsel;

    std::cerr << "vector size ="<<part_data->id.size() <<"  nsel="<<nsel<<"\n";


    long ax(image->axis(0));
    long ay(image->axis(1));
    long az;
    if (image->axes()>=3) {
        az=(image->axis(2));
    } else {
        az=0;
    }

    // MUST BE SET from CLI or GUI
    float dmin = std::numeric_limits<float>::min();
    float dmax = std::numeric_limits<float>::max();
    if (go->phys_min_glob!=-1) {
        //dmin = go->phys_min_glob;
    }
    if (go->phys_max_glob!=-1) {
        //dmax = go->phys_max_glob;
    }
    int   zmin = std::numeric_limits<int>::min();
    int   zmax = std::numeric_limits<int>::max();

    // Load data and fill up arrays
    long nan=0;
    int valid_nbody=0;
    int cpt=0;
    for (unsigned long i=0; i<contents.size(); i++) {

        if (std::isfinite(contents[i]) && contents[i]>=dmin && contents[i]<=dmax) {

            long z_i=int(i/(ax*ay)); // current Z plane
            assert(z_i<=az);
            if (z_i>=zmin && z_i<=zmax) { // inside Z selection
                if (index_tab[cpt]!=-1) {
                    long nxy=i-z_i*(ax*ay);// #pixels (x/y) of the latest Z plane
                    long y_i=int(nxy/ax);  // current Y coordinate
                    assert(y_i<ay);
                    long x_i=nxy-(y_i*ax); // current X coordinate
                    assert(z_i<=az);
                    part_data->pos[valid_nbody*3+0] = x_i*1.0;
                    part_data->pos[valid_nbody*3+1] = y_i*1.0;
                    part_data->pos[valid_nbody*3+2] = z_i*1.0;
                    part_data->rho->data[valid_nbody]     = contents[i];
                    part_data->rneib->data[valid_nbody]   = 1.0;
                    part_data->id.push_back(valid_nbody);
                    valid_nbody++;
                }
                cpt++;
                assert(valid_nbody<=nsel);
            }
            //std::cerr << x_i << " " << y_i << " " << z_i << " " << contents[i] << "\n";
        } else { // NAN or INF
            nan++;
        }
    }
    glnemo::ComponentRange cr;
    *part_data->nbody = valid_nbody;
    part_data->computeVelNorm();
    part_data->rho->computeMinMax();
    part_data->rneib->computeMinMax();

    if (! part_data->timu ) part_data->timu = new float;
    *part_data->timu = 0.0;
    end_of_data=true; // only one frame from an gadget snapshot
  }
  return status;
}
// ============================================================================
// close()
int SnapshotFits::close()
{
  //fits_io->close();
  if (valid) {

  }
  end_of_data=false;
  return 1;
}
// ============================================================================
// endendOfDataMessage()
QString SnapshotFits::endOfDataMessage()
{
  QString message=tr("Fits Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
}
// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotFits

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(fitsplugin, glnemo::SnapshotFits);
#endif

