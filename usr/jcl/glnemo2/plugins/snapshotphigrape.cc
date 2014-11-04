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
#include <QtGui> // Mandatory for plugins management
#include "snapshotphigrape.h"
#include "zlib.h"
#include <cstdio>
#include <sstream>

#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
Q_PLUGIN_METADATA(IID "fr.glnemo2.phigrapePlugin")
#endif

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
  size_buff = 2048;
  BUFF = NULL;
  sbuff="";
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

  close();
  if (BUFF) delete [] BUFF;

}
// ============================================================================
// newObject                                                                   
// instantiate a new object and return a pointer on it                         
SnapshotInterface * SnapshotPhiGrape::newObject(const std::string _filename, const int x)
{
  if (x) {;} // get rid of compiler warning
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
// readBuffer()
bool SnapshotPhiGrape::readBuffer()
{
  bool status = false;
  int bytes_read = gzread(file,&BUFF[0],size_buff-1); // read the buffer
  if (bytes_read == -1) { // reading error
    std::cerr << "SnapshotPhiGrape::readBuffer error on gzread\n";
    throw 11;
  }
  if (bytes_read==0) {
    std::cerr << "SnapshotPhiGrape::readBuffer END OF FILE\n";
    throw 12;
  }
  if (bytes_read>0) {
    BUFF[bytes_read]='\0';
    sbuff = sbuff + string(BUFF); // merge old + new buffer
    status=true;
  }
  return status;
}
// ============================================================================
// gzGetString();
bool SnapshotPhiGrape::gzGetLine()
{
  bool status=false;
  if (sbuff.size()>0) {     // buffer is not empty      
    status=getLine();   // get a line
  }
  if (!status) { // no line found or buffer is empty
    if ((readBuffer()      // then  append a new buffer
      && getLine())) { // and   get a new line
      status= true;
    }
    else {
      std::cerr << "SnapshotPhiGrape::gzGetLine can't find a new line \n";
      throw 13;
    }
  }
  return status;
}
#if 0
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
    bool first=true;
    for (int i=0; i<full_nbody; i++) {
  
      //gzgets(file, buff, KB);   // read a line from the file
      gzGetLine();
      int idx=index_tab[i];
      if (idx!=-1) { // it's a valid particle
        int nn; 
        float v[3],dummy;
        float rho1,hsml;
	int ret=0;
        //                 n  m  x  y  z  vx vy vz rho hsml
	ret = sscanf(line,"%d %f %f %f %f %f %f %f %f  %f",
		     &nn,&dummy,
		     &part_data->pos[cpt*3+0],
		     &part_data->pos[cpt*3+1],
		     &part_data->pos[cpt*3+2],
		     v,v+1,v+2,
		     &rho1,&hsml);
        if (load_vel) 
          for (int j=0; j<3; j++)
            part_data->vel[cpt*3+j] = *(v+j);

        if (first && ret>=10) { // first time and rho+hsml exist
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
        if (ret>=10) {
          part_data->rho->data[cpt]   = rho1;
          part_data->rneib->data[cpt] = hsml;
        }
        cpt++;
      }
      
    }
    if (part_data->rho) part_data->rho->computeMinMax();
  }
  end_of_data = true;
  return status;
}
#else
// ============================================================================
// nextFrame()                                                                 
int SnapshotPhiGrape::nextFrame(const int * index_tab, const int nsel)
{
  
  // note about strtof
  //
  // after  x=strtof(p, &endptr);
  // if (p==endptr) means that there was NOT data to read
  
  int status=0;
  if (valid) {
    BUFF = new char[size_buff];
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
    bool first=true;
    for (int i=0; i<full_nbody; i++) {
  
      gzGetLine();
      int idx=index_tab[i];
      if (idx!=-1) { // it's a valid particle
        long nn;
        float v[3],dummy;
        float rho1,hsml;
        char *endptr;
        char * p = line;
        nn     = strtol(p, &endptr, 10); // index
        p   = endptr;
        dummy  = strtof(p, &endptr); // mass
        p   = endptr;
        // pos
        for (int i=0;i<3;i++) {
            part_data->pos[cpt*3+i] = strtof(p, &endptr);
            p = endptr;
        }
        // vel
        for (int i=0;i<3;i++) {
            *(v+i) = strtof(p, &endptr);
            p = endptr;
        }
        if (load_vel)
          for (int j=0; j<3; j++)
            part_data->vel[cpt*3+j] = *(v+j);

        // try rho
        rho1  = strtof(p, &endptr);
        // rho ALLOC
        if (p!=endptr && first && *endptr!='\0') {
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
        // rho
        if (p!=endptr && *endptr!='\0') {
            //rho1  = strtof(p, &endptr);
            p   = endptr;
            part_data->rho->data[cpt]   = rho1;
            // hsml
            hsml  = strtof(p, &endptr);
            
            if (p != endptr) {                
                p   = endptr;
                part_data->rneib->data[cpt] = hsml;
            }
            else { // default if not exist
                part_data->rneib->data[cpt] = 1.;
            }
        }
        cpt++;
      }
      
    }
    if (part_data->rho) part_data->rho->computeMinMax();
  }
  end_of_data = true;
  return status;
}
#endif
// ============================================================================
// close()                                                                     
int SnapshotPhiGrape::close()
{
  int ret=0;
  if (file) gzclose(file);
  if (valid) {    
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

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(phigrapeplugin, glnemo::SnapshotPhiGrape);
#endif
