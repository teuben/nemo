// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "cpart.h"
#include <QDir>

namespace ramses {
  using namespace std;
// ----------------------------------------------------------------------------
// READING constructor                                                                 
CPart::CPart(const std::string _indir, const int _select,const bool _v)
{
  valid=false;
  nbody     = 0;  
  ndm       = 0;
  ndm_box   = 0;
  nstar_box = 0;
  nselect   = 0;
  verbose=_v;
  select = _select;
  indir = _indir;
  infile="";
  
  // keep filename untill last /
  int found=indir.find_last_of("/");
  if (found != (int) string::npos && (int) indir.rfind("output_")<found) {
    indir.erase(found,indir.length()-found);
  }
  std::cerr << "indir =" << indir <<"\n";
  
  found=(int) indir.rfind("output_"); 
  if (found!=std::string::npos) {
    s_run_index= indir.substr(found+7,indir.length()-1); // output_ = 7 characters
    
    while ((found=s_run_index.find_last_of("/"))>0) { // remove trailing "/"
      s_run_index.erase(found,found);
    }
    std::cerr << "Run index = " << s_run_index << "\n";
    infile = indir + "/part_" + s_run_index + ".out00001";
    std::cerr << "infile =" << infile <<"\n";
  }
}

// ----------------------------------------------------------------------------
// Destructor                                                                 
CPart::~CPart()
{

}
// ----------------------------------------------------------------------------
//
bool CPart::isValid()
{    
  if (part.open(infile)) {
    valid=true;
    readHeader();
    part.close();
  }
  else
    valid=false;
  return valid;
}

// ============================================================================
// readHeader
int CPart::readHeader()
{
  part.readDataBlock((char *) &ncpu);
  part.readDataBlock((char *) &ndim);
  part.readDataBlock((char *) &npart);
  part.skipBlock();
  part.readDataBlock((char *) &nstar);  
  return 1;
}
// ============================================================================
// loadData
int CPart::loadData(bool take_halo, bool take_stars,float * pos, float * vel,const int *index,
                    const int nsel, const bool load_vel, const int namr_box)
{
  int cpt_dm=0;
  int cpt_star=0;
  bool count_only=false;
  QString str_status;
  
  if (nsel) {;} // remove compiler warning
  
  if (index==NULL) {
    ndm_box=0;
    nstar_box=0;
    count_only=true;  
    nselect=0;
  }
  else {
    //assert(nsel==nselect);
  }
  nbody=0;
  
  for (int i=0; i<ncpu; i++) {
    std::ostringstream osf;
    osf << std::fixed << std::setw(5) << std::setfill('0') <<i+1;
    std::string infile = indir + "/part_" + s_run_index + ".out" + osf.str();
    if (verbose) std::cerr << "reading file : " << infile << "\n";
    str_status = std::string("Loading file : " + infile).c_str();
    emit stringStatus(str_status);
    part.open(infile);
    readHeader();
    
    double * tmp[6];//=new double[npart];
    part.skipBlock(3);
    
    // read positions
    
    for (int j=0; j<ndim; j++) {
      tmp[j] = new double[npart]; // alloc
      part.readDataBlock((char *) tmp[j]);
    }
    
    
    // read velocities
    for (int j=0; j<ndim; j++) {
      tmp[3+j] = new double[npart]; // alloc
      if (load_vel) {
        part.readDataBlock((char *) tmp[3+j]);
      } else {
        part.skipBlock();
      }
    }
    
    // skip masses
    //tmp[6] = new double[npart]; // alloc
    //part.readDataBlock((char *) tmp[6]);
    part.skipBlock();
    
    double * age;
    if (nstar>0) { // there are stars
      part.skipBlock(); // skip identity
      part.skipBlock(); // skip level
      age = new double[npart];
      part.readDataBlock((char *) age);
      for (int k=0; k<npart; k++) {
        if ((age[k]==0.&& (select==0 || select==2))  || // it's DM    && (DM    sel || DM + Stars sel)
            (age[k]!=0.&& (select==1 || select==2))) {  // its' stars && (Stars sel || DM + Stars sel)
          if ((tmp[0][k]>=xmin && tmp[0][k]<=xmax) &&
              (tmp[1][k]>=ymin && tmp[1][k]<=ymax) &&
              (tmp[2][k]>=zmin && tmp[2][k]<=zmax)) {
            
            if (count_only) {
              if (age[k]==0) { // it's DM
                ndm_box++;
              } else {         // it's a star
                nstar_box++;
              }
              nselect++;
            }
            else {
              int idx=666;//index[nbody];
              if (idx!=-1) { // it's a valid particle
                if (take_halo && age[k]==0) { // DM selected and it's a  DM
                  int cpt = namr_box+cpt_dm;
                  assert(cpt<(nselect+namr_box));
                  for (int l=0;l<3;l++) {
                    pos[cpt*3+l]=tmp[l][k];
                    if (load_vel) {
                      vel[cpt*3+l]=tmp[3+l][k];
                    }
                  }
                  cpt_dm++;
                } 
                if (take_stars && age[k]!=0) { // STARS selected and it's a star
                  int cpt = namr_box+(take_halo?ndm_box:0)+cpt_star;
                  assert(cpt<(nselect+namr_box));
                  for (int l=0;l<3;l++) {
                    pos[cpt*3+l]=tmp[l][k];
                    if (load_vel) {
                      vel[cpt*3+l]=tmp[3+l][k];
                    }
                  }
                  cpt_star++;
                }
                
              }
              nbody++;
            }
          }
        }                
      }
      // garbage
      delete [] age;
    }
    else {  // there are no stars
      if (select==0 || select==2) { // DM sel|| DM + stars sel        
        for (int k=0; k<npart; k++) {
          if ((tmp[0][k]>=xmin && tmp[0][k]<=xmax) &&
              (tmp[1][k]>=ymin && tmp[1][k]<=ymax) &&
              (tmp[2][k]>=zmin && tmp[2][k]<=zmax)) {
            
            if (count_only) {
              if (1 /*age[k]==0*/) { // it's DM
                ndm_box++;
              } else {         // it's a star
                nstar_box++;
              }
              nselect++;
            }
            else {
              int idx=666;//index[nbody];
              if (idx!=-1) { // it's a valide particle
                if (take_halo/* && age[k]==0*/) { // DM selected and it's a  DM
                  int cpt = namr_box+cpt_dm;
                  assert(cpt<(nselect+namr_box));
                  for (int l=0;l<3;l++) {
                    pos[cpt*3+l]=tmp[l][k];
                    if (load_vel) {
                      vel[cpt*3+l]=tmp[3+l][k];
                    }
                  }
                  cpt_dm++;
                } 
                if (take_stars /*&& age[k]!=0*/) { // STARS selected and it's a star
                  int cpt = namr_box+(take_halo?ndm_box:0)+cpt_star;
                  assert(cpt<(nselect+namr_box));
                  for (int l=0;l<3;l++) {
                    pos[cpt*3+l]=tmp[l][k];
                    if (load_vel) {
                      vel[cpt*3+l]=tmp[3+l][k];
                    }
                  }
                  cpt_star++;
                }
                
              }
              nbody++;  
            }
          }
        }
      }
      // garbage
      //delete [] age;
    } 
    // garbage collecting
    for (int i=0; i<6; i++)
        delete [] tmp[i];    
     
    part.close(); // close current file  
  } // for ... 
  return nselect;
}
} // namespace ramses
