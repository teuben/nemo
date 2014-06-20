// ============================================================================
// Copyright Jean-Charles LAMBERT - 2011                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
#include <iostream>
#include <cmath>
#include <sstream>
#include "cbar.h"
#include "uns.h"
#include "csnaptools.h"

using namespace uns_proj;
using namespace std;
#define PI 3.141592653589793238462643
// ----------------------------------------------------------------------------
// contructor
CBar::CBar(const int _nbody, float * _pos, float * _vel, float * _mass, float *_rho, float *_hsml)
{
  // link arguments
  nbody = _nbody;
  pos = _pos;
  vel = _vel;
  mass= _mass;
  density = NULL;
  rho = _rho;
  hsml=_hsml;
  // sort by rho
  sortRho();
}
// ----------------------------------------------------------------------------
// contructor
CBar::~CBar()
{
  if (density) delete density;
}
// ----------------------------------------------------------------------------
// computeAngle();
// Based on anglebar_impr.F
// with new method using density instead of mass
float CBar::computeAngle(const float dmin, const float dmax, const bool mvcod)
{  
  // shift to COD
  double cod[6]={0.,0.,0.,0.,0.,0.};
  if (mvcod)
    CSnaptools::moveToCod<float>(nbody,&pos[0],&vel[0],&mass[0],rho,cod,false);
  
  float alpha=0.0;
  float beta =0.0;
  float minrho=log(rho[vec_rho.at(0).index]);
  float maxrho=log(rho[vec_rho.at(nbody-1).index]);
  float binf=minrho+dmin*(maxrho-minrho);
  float bsup=minrho+dmax*(maxrho-minrho);
  std::cerr << "binf/bsup :"<< binf<<"/"<<bsup<<"\n";
  int cpt=0;
  //for (int i=dmin*nbody; i<dmax*nbody; i++) {
  for (int i=0; i<nbody;i++) {
    int ii= vec_rho.at(i).index;
    float rhoii = log(rho[ii]);
    if ((rhoii)>=binf && (rhoii)<=bsup ) {
      cpt++;
      float x = pos[ii*3+0]-cod[0];
      float y = pos[ii*3+1]-cod[1];
      float x2 = x*x;
      float y2 = y*y;
      float r2 = x2+y2;
      float cos2th = (x2 - y2) / r2;
      float sin2th = 2. * x * y / r2;
#if 1
      alpha += rho[ii]*sin2th;
      beta  += rho[ii]*cos2th;    
#else
      alpha += mass[ii]*sin2th;
      beta  += mass[ii]*cos2th;
#endif
    }
  }
  std::cerr << "Found ["<<cpt<<"] particles into the range.\n";
  assert(cpt>0);
  float bar_angle = 0.5 * atan2(alpha,beta);
  return bar_angle;  
}
// ----------------------------------------------------------------------------
// computeAngle();
// Based on anglebar_impr.F
// with new method using density instead of mass
// select max #particles per density shells
// guess rotation angle from selected particles
float CBar::computeAngle(const bool mvcod)
{  
  float minrho=log(rho[vec_rho.at(0).index]);
  float maxrho=log(rho[vec_rho.at(nbody-1).index]);
  
  // reset array
  for (int i=0; i<100; i++) {
    data_histo[i] = 0;
  }
  // compute #particles per density shell
  for (int i=0; i<nbody;i++) {
    int ii= vec_rho.at(i).index;
    int index=((log(rho[ii])-minrho)*99.)/(maxrho-minrho);
    //std::cerr << "CBar::computeAngle index="<<index<<"\n";
    assert(index<100);
    data_histo[index]++;
  }
  // find max shell
  int maxshell=data_histo[0];
  int ishell=0;
  for (int i=1;i<100;i++) {
    if (data_histo[i]>maxshell) {
      maxshell = data_histo[i];
      ishell=i;
    }    
  }
  
  float dmax=std::max(ishell+5,ishell);
  float dmin=std::max(0.,(ishell-20.));
  std::cerr << "CBar::computeAngle dmin="<<dmin<<"/ dmax="<<dmax<<"\n";
  return computeAngle(dmin/100.,dmax/100.,mvcod);
}
// ----------------------------------------------------------------------------
// rotateOnX
void CBar::rotateOnX(const float bar_angle)
{
  float theta = -bar_angle;
  rotate(theta);
}
// ----------------------------------------------------------------------------
// rotateOnY
void CBar::rotateOnY(const float bar_angle)
{
  float theta =  0.5*M_PI-bar_angle;
  rotate(theta);
}
// ----------------------------------------------------------------------------
// rotate
void CBar::rotate(const float theta)
{
  for (int i=0; i<nbody; i++) {
    float rx = pos[i*3+0] * cos(theta) - pos[i*3+1] * sin(theta);
    float ry = pos[i*3+0] * sin(theta) + pos[i*3+1] * cos(theta);
    pos[i*3+0] = rx;
    pos[i*3+1] = ry;
  }
  if (vel) { // velocities exist
    for (int i=0; i<nbody; i++) {
      float rx = vel[i*3+0] * cos(theta) - vel[i*3+1] * sin(theta);
      float ry = vel[i*3+0] * sin(theta) + vel[i*3+1] * cos(theta);
      vel[i*3+0] = rx;
      vel[i*3+1] = ry;   
    }
  }
}
// ----------------------------------------------------------------------------
//  save
void CBar::save(std::string out, const float timu, const bool mvcod)
{
  // shift to COD
  double cod[6]={0.,0.,0.,0.,0.,0.};
  if (mvcod)
    CSnaptools::moveToCod<float>(nbody,&pos[0],&vel[0],&mass[0],rho,cod,true);
  // save data
  uns::CunsOut * unsout = new uns::CunsOut(out,"nemo",false);   
  unsout->snapshot->setData("time",timu);
  unsout->snapshot->setData("mass",nbody,mass,false);
  unsout->snapshot->setData("pos" ,nbody,pos ,false);
  unsout->snapshot->setData("vel" ,nbody,vel ,false);
  unsout->snapshot->setData("rho" ,nbody,rho ,false);
  unsout->snapshot->setData("hsml",nbody,hsml,false);
  unsout->snapshot->save();
  delete unsout;
}

// ----------------------------------------------------------------------------
//  saveAllRho
void CBar::saveAllRho(std::string out)
{
  int cpt=0;
  for (int i=0;i<99;i++) {
    int start=nbody*(float)(i  )/100.;
    int end  =nbody*(float)(i+1)/100.;
    end = std::min(end,nbody);
    int n=end-start;
    if (n>0) {
      float * p = new float[n*3]; // pos
      float * r = new float[n];   // rho
      float * h = new float[n];   // hsml
      int ii=0;
      for (int i=start; i<start+n; i++) {
      //for (int i=0; i<end; i++) {
        int indx = vec_rho.at(i).index;
        p[ii*3+0] = pos[indx*3+0];
        p[ii*3+1] = pos[indx*3+1];
        p[ii*3+2] = pos[indx*3+2];
        r[ii]     = rho[indx];
        h[ii]     = hsml[indx];
        ii++;
      }
      assert(ii==n);
      stringstream ss;
      ss << out << "." << setw(5) << setfill('0') << cpt++;
      uns::CunsOut * unsout = new uns::CunsOut(ss.str(),"nemo",false);   
      unsout->snapshot->setData("pos" ,n,p,false);
      unsout->snapshot->setData("rho" ,n,r,false);
      unsout->snapshot->setData("hsml",n,h,false);
      unsout->snapshot->save();
      delete unsout;
      delete [] p;
      delete [] r;
      delete [] h;
    }
  }
}
// ----------------------------------------------------------------------------
// sortRho();
void CBar::sortRho()
{
  if (!rho) {
    std::cerr << "Density NULL during instantiation, we gonna compute density!!\n";
    // Instantiate a density object  
    density = new CDensity(nbody,&pos[0],&mass[0]);
    density->compute(0,32,1,8); // estimate density
    rho = density->getRho();
    hsml= density->getHsml();
  }
  
  // put particles into a vector
  vec_rho.clear();
  vec_rho.reserve(nbody);
  for (int i=0;i<nbody;i++) {
    CVecRho p(this,i);
    vec_rho.push_back(p);
  }
  
  // Descending sort particles according to their densities
  std::sort(vec_rho.begin(),vec_rho.end(),CVecRho::sortRho);
}
// ----------------------------------------------------------------------------
// VecRho
bool CVecRho::sortRho(const CVecRho& a, const CVecRho& b) 
{
  return a.bar->getRho()[a.index] <  b.bar->getRho()[b.index];      
}

