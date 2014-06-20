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
#include <cpgplot.h>
#include <iostream>
#include <cmath>
#include "cfitsellipse.h"
#include "csnaptools.h"
#include "uns.h"
#include "cgaussian.h"
using namespace uns_proj;
using namespace jclut;
using namespace std;

// ----------------------------------------------------------------------------
// contructor
CFitsEllipse::CFitsEllipse(const int _xaxis, const int _yaxis, const int nm, const float _tmax)
{
  tmax = _tmax;
  nmesh   = nm;
  xaxis       = _xaxis;
  yaxis       = _yaxis;
  
  // create grid
  grid = new float[nmesh*nmesh];
}
// ----------------------------------------------------------------------------
// destructor
CFitsEllipse::~CFitsEllipse()
{
  delete [] grid;
}
// ----------------------------------------------------------------------------
// buildGrid()
void CFitsEllipse::buildGrid(const int nbody, const float * xi, const float * val)
{
  for (int i=0; i<nmesh*nmesh; i++) 
    grid[i] = 0.0;
  // find min value
  float min = CSnaptools::minArray<float>(nbody,val);
  std::cerr << "buildGrid minimum="<<min<<"\n";
  int sm=1;
  for (int i=0; i<nbody; i++) {
    int x=(xi[3*i+xaxis ]/tmax+1)*nmesh/2-.5;
    int y=(xi[3*i+yaxis ]/tmax+1)*nmesh/2-.5;
    if (x>sm && x<nmesh-sm && y>=sm && y<nmesh-sm) {
      // smoothing
      for (int j=-sm;j<=sm;j++) {
        for (int k=-sm;k<=sm;k++) {
          grid[(y+j)*nmesh+(x+k)] += val[i]/min;          
        }
      }
    }
  }
  // compute the log
  for (int i=0; i<nmesh*nmesh; i++) {
    grid[i] = log(min+grid[i]);
  }
#if 0
  CGaussian<float> * gaussian = new CGaussian<float>(20,5.);
  for (int i=0;i<nmesh;i++) {
    for (int j=0;j<nmesh;j++) {        
        gaussian->applyOnArrayXY(grid,nmesh,nmesh,j,i,1.0);
      }
  }
#endif
}
// ----------------------------------------------------------------------------
// intensity
float CFitsEllipse::intensity(const float x, const float y)
{  
  float tmin=-tmax;
  
  // to transpose in alignment with the table
  int   ix=(x-tmin)*nmesh/(tmax-tmin)-1;
  int   iy=(y-tmin)*nmesh/(tmax-tmin)-1;
  float xi=tmin+(ix+1)*(tmax-tmin)/nmesh;
  float yi=tmin+(iy+1)*(tmax-tmin)/nmesh;
  // if out of limits, return zero */
  if(ix<1||ix>=nmesh-1||iy<1||iy-1>=nmesh)return 0.0;
  
  float p=nmesh*(x-xi)/(tmax-tmin);
  float q=nmesh*(y-yi)/(tmax-tmin);
  
  float f00=grid[ ix   +nmesh*    iy];
  float f10=grid[(ix+1)+nmesh*    iy];
  float f01=grid[ ix   +nmesh*(iy+1)];
  float f11=grid[(ix+1)+nmesh*(iy+1)];
  
  float rimage=(1-p)*(1-q)*f00+
               p*(1-q)*f10+
               q*(1-p)*f01+
               p*q*f11;
  return rimage;
}

#define CLEV 20
// ----------------------------------------------------------------------------
// displayGrid()
void CFitsEllipse::displayGrid()
{
  //
  std::string outdev="/xw";
  cpgopen(outdev.c_str());
  
  float tr[6];
  float lev1[CLEV];
  
  tr[0] = -tmax ; tr[1] = 2 * (tmax)/nmesh  ; tr[2]=0             ;
  tr[3] = -tmax ; tr[4] = 0                 ; tr[5]=2*(tmax)/nmesh;
  
  float maxlev = 0.95 * intensity(0,0.1);
  float minlev= 0.25 * (intensity(tmax/2, 0.) + 
                        intensity(tmax/2, 0.) + 
                        intensity(0., tmax/2) + 
                        intensity(0., tmax/2));
  std::cerr << "minlev = "<<minlev<< " maxlevel="<<maxlev<<"\n";
  for (int i=0;i<CLEV;i++) 
    lev1[i]=minlev+(maxlev-minlev)*i/CLEV;
  
  //cpgsvp (0.05, 0.40, 0.1, 0.9);
  cpgsvp (0.01, 0.01, 0.01, 0.99);
  cpgwnad (-tmax/2, tmax/2, -tmax/2, tmax/2);
  cpgsls(1);
  cpgcont (grid, nmesh, nmesh, 1, nmesh, 1, nmesh, lev1, CLEV, tr);
  cpgbox ("BCTN", 0.0, 0, "BCTN", 0.0, 0);  
  cpglab ("x","y","");
  
#if 0
  cpgsch(0.8);
  cpgmtxt("T",5.0,0.5,0.5,simname.c_str());
  cpgmtxt("T",7.0,0.5,0.5,select_c);
  sprintf(label,"Bulge radius : %f (%f)",bulgerad,Rb);
  cpgmtxt("T",3.0,0.5,0.5,label);
  cpgsch(1);
#endif
  cpgask(1);
  cpgend();
}
// ----------------------------------------------------------------------------
// saveGrid()
void CFitsEllipse::saveGrid(std::string outname)
{

  uns::CunsOut * unsout = new uns::CunsOut(outname,"nemo",false); 
  
  float * pos  = new float[nmesh*nmesh*3];
  float * hsml = new float[nmesh*nmesh  ];
  
  for (int i=0;i<nmesh;i++) {
    for (int j=0;j<nmesh;j++) {        
        pos[i*nmesh*3+j*3 +0] = j;
        pos[i*nmesh*3+j*3 +1] = i;
        pos[i*nmesh*3+j*3 +2] = 0;
        hsml[i*nmesh+j] = 1.0;
      }
  }
  
  unsout->snapshot->setData("pos"  ,nmesh*nmesh,pos,false);
  unsout->snapshot->setData("rho"  ,nmesh*nmesh,grid ,false);
  unsout->snapshot->setData("hsml" ,nmesh*nmesh,hsml ,false);
  unsout->snapshot->save();
  
  delete [] hsml;
  delete [] pos;
}
