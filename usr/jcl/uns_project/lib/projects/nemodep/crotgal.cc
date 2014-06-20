// ============================================================================
// Copyright Jean-Charles LAMBERT - 2010                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#include <iostream>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <cstring>
#include "crotgal.h"
#include "csnaptools.h"
#include "uns.h"

using namespace uns_proj;
using namespace jclut;
using namespace std;
using namespace uns;
// ----------------------------------------------------------------------------
// contructor
CRotgal::CRotgal(uns::CunsIn * _uns)
{
  unsin = _uns;
  nbody = 0;
  density = NULL;
}
// ----------------------------------------------------------------------------
// Destructor
CRotgal::~CRotgal()
{
  clearVectors();
  pvec.clear();
  prtvec.clear();
  pvecselect.clear();
  delete density;
}

// ----------------------------------------------------------------------------
// loadData
bool CRotgal::loadData()
{
  bool ret=false;
  if (unsin->snapshot->nextFrame()) {
    clearVectors(); // flush vectors
    ret=true;
    bool ok;
    int nn;
    float * p;
    int   * pid;
    // get the input number of bodies according to the selection
    ok=unsin->snapshot->getData("nsel",&nbody);
    assert(ok==true);
    // get the simulation time
    ok=unsin->snapshot->getData("time",&time);
    std::cerr << "nbody=" << nbody << " time="<<time <<"\n";
        
    // get POS from input snapshot
    ok=unsin->snapshot->getData("pos" ,&nn,&p);
    assert(ok==true);
    pos.reserve(sizeof(float)*nbody*3);       // reserve vector space
    memcpy(&pos[0],p,sizeof(float)*nbody*3);  // cp memory to vectors
    
    // get VEL from input snapshot
    ok=unsin->snapshot->getData("vel" ,&nn,&p);
    if (ok) {
      vel.reserve(sizeof(float)*nbody*3);       // reserve vector space
      memcpy(&vel[0],p,sizeof(float)*nbody*3);  // cp memory to vectors
    }
    // get MASS from input snapshot
    ok=unsin->snapshot->getData("mass",&nn,&p);
    if (ok) {
      mass.reserve(sizeof(float)*nbody);       // reserve vector space
      memcpy(&mass[0],p,sizeof(float)*nbody);  // cp memory to vectors
    }
    // get IDs from input snapshot
    ok=unsin->snapshot->getData("id",&nn,&pid);
    if (ok) {
      id.reserve(sizeof(int)*nbody);       // reserve vector space
      memcpy(&id[0],pid,sizeof(int)*nbody);  // cp memory to vectors
    }
    assert(ok==true);
  }
  return ret;
}
// ----------------------------------------------------------------------------
// process
void CRotgal::process()
{
  if (!density) delete density;
  // Instantiate a density object  
  density = new CDensity(nbody,&pos[0],&mass[0]);
  density->compute(0,32,1,8); // estimate density
  
  // shift to COD
  double cod[6];
  float * vv=NULL;
  if (vel.size()>0) {
    vv = &vel[0];
  }
  CSnaptools::moveToCod<float>(nbody,&pos[0],vv,&mass[0],density->getRho(),cod,true);
  
  // put particles into a vector
  pvec.clear();
  pvec.reserve(nbody);
  for (int i=0;i<nbody;i++) {
    CPartVec p(this,i);
    pvec.push_back(p);
  }
  // Descending sort particles according to their densities
  sortRho();
  
}
// ----------------------------------------------------------------------------
// selectPart()
void CRotgal::selectPart()
{
  
  pvecselect.clear();
  //for (int i=0.2*nbody; i<0.25*nbody; i++) {
  //for (int i=0; i<0.05*nbody; i++) {
  for (int i=0.4*nbody; i<0.45*nbody; i++) {
    int ii= pvec.at(i).index;
    CPartVec p(this,ii);
    p.computeR2();
    pvecselect.push_back(p);
  }
  
  // Descending sort particles according to their Ids
  std::sort(pvecselect.begin(),pvecselect.end(),CPartVec::sortId);
}
// ----------------------------------------------------------------------------
// selectPart()
void CRotgal::saveSelectPart(std::vector <CPartVec> * ppvec)
{
  // clear RT vector
  prtvec.clear();
  
  // Descending sort particles according to their Ids
  std::sort(pvec.begin(),pvec.end(),CPartVec::sortId);
    
  int cpt=0;
  std::cerr << "selectPart : pvec.size "<< pvec.size() << " ppvec->size() = "<<ppvec->size()<<"\n";
  int start=0;
  int last=0;
  // loop on particles selected from the previous time
  for (int j=0; j<(int) ppvec->size(); j++){ 
    int jidx = ppvec->at(j).index;
    bool stop=false;
    start=last;
    // loop on particles selected from the current time
    for (int i=start; i<(int)pvec.size()&&!stop; i++) {
      int iidx = pvec.at(i).index;      
      
      //std::cerr <<  ppvec->at(j).rotgal->id[jidx] << " " << pvec.at(i).rotgal->id[iidx]<< "\n";
      if (ppvec->at(j).rotgal->id[jidx] == pvec.at(i).rotgal->id[iidx]) { // matches!!
        
        //std::cerr << j<<" FOUND!!!\n";
        stop=true; // we stop to speedup

        // compute the new position of the particles
        computeRadiusTheta(&ppvec->at(j),&pvec.at(i));
        
        cpt++;
        last=i; // we position the last pointer to speedup
      }
    }
  }  
}
// ----------------------------------------------------------------------------
// computeTheta
#define   PI         3.141592653589793238462643
void CRotgal::computeRadiusTheta(CPartVec * p1, CPartVec * p2)
{

  float r1=sqrt(p1->computeR2());  // old radius
  float r2=sqrt(p2->computeR2());  // new radius
  
  //p1
  float x1=p1->rotgal->pos[p1->index*3+0];
  float y1=p1->rotgal->pos[p1->index*3+1];
  float z1=p1->rotgal->pos[p1->index*3+2];
  //p2
  float x2=p2->rotgal->pos[p2->index*3+0];
  float y2=p2->rotgal->pos[p2->index*3+1];
  float z2=p2->rotgal->pos[p2->index*3+2];
  
  float c=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
  float theta = acos((r1*r1+r2*r2-c*c)/(2*r1*r2));
  
  float diff_radius = fabs(r2-r1)*100./r1;
  //std::cerr << "Diff radius="<<diff_radius<<" theta="<<theta<<" "<< theta*180/PI<<"\n";
  CPartRT p(diff_radius,theta);
  prtvec.push_back(p);
}
// ----------------------------------------------------------------------------
// computeRotation
void CRotgal::computeRotation()
{
  // sort by radius
  sort(prtvec.begin(),prtvec.end(),CPartRT::sortRadius);
  //sort(prtvec.begin(),prtvec.end(),CPartRT::sortTheta);
  
  for (std::vector<CPartRT>::iterator p=prtvec.begin(); p<prtvec.end(); p++) {
    std::cerr << "Diff radius="<<(*p).diff_radius<<" theta="<<(*p).theta<<" "<< (*p).theta*180/PI<<"\n";
    
  }
}

// ----------------------------------------------------------------------------
// selectPart()
void CRotgal::saveSelectPart(std::string out, std::vector <CPartVec> * ppvec)
{
  
  // Descending sort particles according to their Ids
  std::sort(pvec.begin(),pvec.end(),CPartVec::sortId);
    
  std::vector<float> pos1;
  pos1.reserve(ppvec->size()*3);
  //float * vel1 = new float[3*pvec.size()];
  std::vector<float> mass1;
  mass1.reserve(ppvec->size());
  std::vector<float> rho1;
  rho1.reserve(ppvec->size());
  std::vector<float> hsml1;
  hsml1.reserve(ppvec->size());
  int cpt=0;
  std::cerr << "selectPart : pvec.size "<< pvec.size() << " ppvec->size() = "<<ppvec->size()<<"\n";
  int start=0;
  int last=0;
  for (int j=0; j<(int) ppvec->size(); j++){      
    int jidx = ppvec->at(j).index;
    bool stop=false;
    start=last;
    //std::cerr << "last = "<< last <<"\n";
    for (int i=start; i<(int)pvec.size()&&!stop; i++) {
      int iidx = pvec.at(i).index;      
      
      //std::cerr <<  ppvec->at(j).rotgal->id[jidx] << " " << pvec.at(i).rotgal->id[iidx]<< "\n";
      if (ppvec->at(j).rotgal->id[jidx] == pvec.at(i).rotgal->id[iidx]) {
        
        //std::cerr << j<<" FOUND!!!\n";
        stop=true;
        int ii= iidx;
        pos1 [cpt*3+0]=pos[ii*3+0];
        pos1 [cpt*3+1]=pos[ii*3+1];
        pos1 [cpt*3+2]=pos[ii*3+2];
        mass1[cpt]    = mass[ii];
        rho1 [cpt]    = density->getRho()[ii];
        hsml1[cpt]    = density->getHsml()[ii];
        cpt++;
        last=i;
      }
    }
  }
  std::cerr << "cpt="<<cpt<<"\n";
  uns::CunsOut * unsout = new uns::CunsOut(out,"nemo",false); 
  unsout->snapshot->setData("time",time);
  unsout->snapshot->setData("pos" ,cpt,&pos1 [0],false);
  unsout->snapshot->setData("mass",cpt,&mass1[0],false);
  unsout->snapshot->setData("rho" ,cpt,&rho1 [0],false);
  unsout->snapshot->setData("hsml",cpt,&hsml1[0],false);
  unsout->snapshot->save();
}
// ----------------------------------------------------------------------------
//                           CPartVec methods
// ----------------------------------------------------------------------------
// sortRho
bool CPartVec::sortRho(const CPartVec& a, const CPartVec& b)
{
  return a.rotgal->getDensity()->getRho()[a.index] > b.rotgal->getDensity()->getRho()[b.index];
}
// ----------------------------------------------------------------------------
// sortId
bool CPartVec::sortId(const CPartVec& a, const CPartVec& b)
{
  return a.rotgal->id[a.index] < b.rotgal->id[b.index];
}
// ----------------------------------------------------------------------------
// computeR2
float CPartVec::computeR2() 
{
  float x=rotgal->pos[index*3+0];
  float y=rotgal->pos[index*3+1];
  float z=rotgal->pos[index*3+2];
  return x*x+y*y+z*z;
}
// ----------------------------------------------------------------------------
//                           CPartRT methods
// ----------------------------------------------------------------------------
// sortRadius
bool CPartRT::sortRadius(const CPartRT& a, const CPartRT& b)
{
  return a.diff_radius < b.diff_radius;
}
// ----------------------------------------------------------------------------
// sortTheta
bool CPartRT::sortTheta(const CPartRT& a, const CPartRT& b)
{
  return a.theta < b.theta;
}
// ----------------------------------------------------------------------------
