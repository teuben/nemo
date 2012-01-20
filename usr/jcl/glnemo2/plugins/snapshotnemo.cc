// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include <QtGui> // Mandatory for plugins management
#include <iostream>
#include "snapshotnemo.h"
#include "particlesdata.h"

extern "C" {
  int io_nemo(const char * , const char *, ...);
# include <stdinc.h>
# include <filestruct.h>
# include <nemo.h>
# include <snapshot/snapshot.h>
};
namespace glnemo {

// ============================================================================
// Constructor                                                                 
SnapshotNemo::SnapshotNemo():SnapshotInterface()
{
  valid=false;
  //part_data = new ParticlesData(ParticlesData::Malloc); // c memory alloc style
  part_data = new ParticlesData();
  interface_type="Nemo";
  first_stream=false;
}

// ============================================================================
// Destructor                                                                  
SnapshotNemo::~SnapshotNemo()
{
  std::cerr << "SnapshotNemo::~SnapshotNemo()\n";
  if (obj)       delete obj;
  if (part_data) delete part_data;
  //if (nbody) free ((int *)   nbody);
  if (pos)   free ((float *) pos);
  if (vel)   free ((float *) vel);
  //if (rho)   free ((float *) rho);
  //if (rneib) free ((float *) rneib);
  if (valid) close();
}
// ============================================================================
// newObject                                                                   
// instantiate a new object and return a pointer on it                         
SnapshotInterface * SnapshotNemo::newObject(const std::string _filename, const int x)
{
  if (x) {;} // get rid of compiler warning
  filename = _filename;
  obj      = new SnapshotNemo();
  obj->setFileName(_filename);
  return obj;
}
// ============================================================================
// isValidData()                                                               
// return true if it's a NEMO snapshot. Standard input '-' is assumed to be a  
// valid NEMO snapshot.                                                        
bool SnapshotNemo::isValidData()
{
  bool status;
  valid=true;

  if (filename == "-") {        // we assume here that "-"
    status=true;                // is standard input and a NEMO stream...
    npos=NULL; nvel=NULL; ntimu=NULL; nrho=NULL; nrneib=NULL;
    nnemobits = NULL; nnbody = NULL; nid = NULL;
    first_stream=true;
    select_part="all";
    select_time="all";
    load_vel=true;
    if (load_vel) { // velocities requested
      status_ionemo=io_nemo(filename.c_str(),"info,float,read,sp,n,pos,vel,dens,aux,k,t,st,b",
		     select_part.c_str(),&nnbody,&npos,&nvel,&nrho,&nrneib,&nid,
		     &ntimu, select_time.c_str(),&nnemobits);
    }
    else {          // velocities NOT requested
      status_ionemo=io_nemo(filename.c_str(),"float,read,sp,n,pos,dens,aux,k,t,st,b",
                     select_part.c_str(),&nnbody,&npos, &nrho,&nrneib,&nid,&ntimu,
		     select_time.c_str(),&nnemobits);
    }
    full_nbody = *nnbody;
    assert(npos!=NULL);
  } else {
    if (isFileExist()) {
      stream str=stropen(filename.c_str(),(char *) "r"); // open NEMO file for reading
      if ( ! str )  status = false; // failed to open
      if (qsf(str)) status = true;  // it's a structured binary file (NEMO)
      else          status = false; // it's not                            
      strclose(str);
      if (status)  {                // it's a NEMO snapshot
	int * ptr=NULL;      // get the full nbody
	ntimu=NULL;nnemobits = NULL; 
	int status1=io_nemo(filename.c_str(),"float,read,n,t,b",&ptr,&ntimu,&nnemobits);
        if (status1 > 0 || status1 == -1) {
          io_nemo(filename.c_str(),"close");
        }
	assert(ptr);
	full_nbody=*ptr;
	free((int *) ptr);
      }
    }
    else
      status=false; // file does not exist
  }
  valid=status;
  if (valid) {
    if ( ! ( *nnemobits & TimeBit)) { // no TimeBit
      time_first = 0.0;
    }
    else time_first = *ntimu;
  }
  return status;
}
// ============================================================================
// getSnapshotRange                                                            
ComponentRangeVector * SnapshotNemo::getSnapshotRange()
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
int SnapshotNemo::initLoading(GlobalOptions * so)
{
  load_vel = so->vel_req;
  select_part="all";
  select_time=so->select_time;
  std::cerr << "SnapshotNemo::initLoading select_time = " << select_time << "\n";
  return 1;
}
// ============================================================================
// nextFrame()                                                                 
int SnapshotNemo::nextFrame(const int * index_tab, const int nsel)
{
  int status;  // io_nemo status
  std::string force_select;
  if (keep_all) force_select = "all";        // we want to select all particles
  else          force_select = select_part;  // we want a substet of particles
 
#if 0
  if (load_vel) { // velocities requested
    status=io_nemo(filename.c_str(),"float,read,sp,n,pos,vel,t,st,b",
                   force_select.c_str(),&part_data->nbody,&part_data->pos,&part_data->vel,
                   &part_data->timu, select_time.c_str(),&part_data->nemobits);
  }
  else {          // velocities NOT requested
    status=io_nemo(filename.c_str(),"float,read,sp,n,pos,t,st,b",
                   force_select.c_str(),&part_data->nbody,&part_data->pos, &part_data->timu,
                   select_time.c_str(),&part_data->nemobits);
  }
#endif
  if ( ! first_stream) { // Normal NEMO file
    npos=NULL; nvel=NULL; ntimu=NULL; nrho=NULL; nrneib=NULL;
    nnemobits = NULL; ;nnbody = NULL; nid=NULL;
    std::cerr << "selected time = "<< select_time << "\n";
    if (load_vel) { // velocities requested
      status=io_nemo(filename.c_str(),"float,read,sp,n,pos,vel,dens,aux,k,t,st,b",
                   force_select.c_str(),&nnbody,&npos,&nvel,&nrho,&nrneib,&nid,
		   &ntimu, select_time.c_str(),&nnemobits);
    }
    else {          // velocities NOT requested
      status=io_nemo(filename.c_str(),"float,read,sp,n,pos,dens,aux,k,t,st,b",
                     force_select.c_str(),&nnbody,&npos, &nrho,&nrneib,&nid,&ntimu,
		     select_time.c_str(),&nnemobits);
    }
  }
  else { // nemo file from standard input "-"
    status = status_ionemo;
  }
  
 // std::cerr << "status next frame ="<< status <<" nsel ="<<nsel<<" part  nbody"<<*part_data->nbody<<" nbody = "<<*nnbody<<"\n";
 // assert(npos);
  if (first_stream || status != 0) {
    if (first_stream) first_stream=false;
    if (status ==  0) end_of_data=true;
    
    if ( ! ( *nnemobits & TimeBit)) { // no TimeBit
      if ( ! part_data->timu )
        part_data->timu = (float *) malloc(sizeof(float));
      std::cerr << "Forcing time to [0.0]\n";
      *(part_data->timu) = 0.0;
    } else {
      *part_data->timu=*ntimu;
    }
    if (status != -2) { // NEMO snapshot must have particles TAG
      // copy data from NEMO array to glnemo
      if (nsel > *part_data->nbody) {
	if (part_data->pos) delete [] part_data->pos;
	part_data->pos = new float[nsel*3];
	if (load_vel && nvel) {
	  if (part_data->vel) delete [] part_data->vel;
	  part_data->vel = new float[nsel*3];
	}
      }
      assert(nsel<=*nnbody);
      //*part_data->nbody=nsel; // !!!! 07 November 2010 !!!!!!
      int cpt=0;
      for (int i=0; i<*nnbody; i++) {
	int idx=index_tab[i];
	//std::cerr << "idx="<<idx<<"\n";
	if (idx!=-1) {
	  for (int j=0; j<3; j++) {
	    part_data->pos[cpt*3+j] = npos[idx*3+j];
	    if (load_vel && nvel) part_data->vel[cpt*3+j] = nvel[idx*3+j];
	  }
	  cpt++;
	}
      }
      assert(cpt==nsel);
      // garbage collecting
      if (npos)     free ((float *) npos);
      if (nvel)     free ((float *) nvel);

      // read density and rneib
      if (( *nnemobits & DensBit) && ( *nnemobits & AuxBit)) {
        cpt=0;
        if (nsel > *part_data->nbody || !part_data->rho) {
          if (part_data->rho) delete part_data->rho;
          //!!!part_data->rho = new float[nsel];
          part_data->rho = new PhysicalData(PhysicalData::rho,nsel);
          
        }
        if (nsel > *part_data->nbody || !part_data->rneib) {
          if (part_data->rneib) delete [] part_data->rneib;
          part_data->rneib = new PhysicalData(PhysicalData::neib,nsel);
        }
        for (int i=0; i<*nnbody; i++) {
          int idx=index_tab[i];
          if (idx!=-1) {
              part_data->rho->data[cpt]   = nrho[idx];
              cpt++;
          }
        }
        assert(cpt==nsel);
        // garbage collecting
        if (nrho)     free ((float *) nrho);
        cpt=0;
        for (int i=0; i<*nnbody; i++) {
          int idx=index_tab[i];
          if (idx!=-1) {
            part_data->rneib->data[cpt] = nrneib[idx];
            //std::cerr << "rneib="<< part_data->rneib[cpt] <<"\n";
            cpt++;
          }
        }
        assert(cpt==nsel);
        part_data->rho->computeMinMax();
        part_data->rneib->computeMinMax();
      }

      // read Ids
      if (*nnemobits & KeyBit) {        
        if (nsel > *part_data->nbody) {
          part_data->id.clear();
          for (int i=0; i<nsel; i++) part_data->id.push_back(-1);
        }
        cpt=0;
        for (int i=0; i<*nnbody; i++) {
          int idx=index_tab[i];
          if (idx!=-1) {
            part_data->id.at(cpt) = nid[idx];
            cpt++;
          }
        }
        assert(cpt==nsel);
        if (nid) free ((int*) nid);
      }      
      // garbage collecting
      if (nrneib)   free ((float *) nrneib   );
      if (ntimu)    free ((float *) ntimu    );
      if (nnbody)   free ((int   *) nnbody   );
      if (nnemobits)free ((int   *) nnemobits);
      //
      *part_data->nbody=nsel; // !!!! 07 November 2010 !!!!!!
      // compute velocity vector norm
      part_data->computeVelNorm();
    }
  } else {
    end_of_data=true;
  }
  return status;
}
// ============================================================================
// close()                                                                     
int SnapshotNemo::close()
{
  int status=0;
  if (valid) {
    status = io_nemo(filename.c_str(),"close");
    end_of_data = false;
    valid = false; // added 2009 June 19th 
  }
  return status;
}
// ============================================================================
// endendOfDataMessage()                                                       
QString SnapshotNemo::endOfDataMessage()
{
  QString message=tr("Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
  return message;
}
} // end of glnemo namespace

// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotNemo                            

Q_EXPORT_PLUGIN2(nemoplugin, glnemo::SnapshotNemo);

//
