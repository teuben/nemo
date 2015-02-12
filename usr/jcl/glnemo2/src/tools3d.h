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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#ifndef GLNEMOTOOLS3D_H
#define GLNEMOTOOLS3D_H
#include "glselection.h"
#include "particlesobject.h"
#include "vec3d.h"
#include <limits>

namespace glnemo {

#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mModel[col*4+row]

class Tools3D{
public:
    Tools3D(double * _m, double * _p);

    ~Tools3D();
    inline Vec3D projPoint(const float x,const float y,const float z) {
      float w=1.;
      // do the product Mmodel X point = mxyzw
      float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3)*w;
      float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3)*w;
      float Mmz= MM(2,0)*x + MM(2,1)*y + MM(2,2)*z;
      float mz = Mmz + MM(2,3)*w;
      float mw = MM(3,0)*x + MM(3,1)*y + MM(3,2)*z + MM(3,3)*w;
      // do the product Mproj X mxyzw  = pxyzw
      float Ppx= MP(0,0)*mx + MP(0,1)*my + MP(0,3)*mw;
      float px = Ppx + MP(0,2)*mz;
      float Ppy= MP(1,0)*mx + MP(1,1)*my + MP(1,3)*mw;
      float py = Ppy + MP(1,2)*mz;
      float Ppz= MP(2,0)*mx + MP(2,1)*my + MP(2,3)*mw;
      float pz = Ppz + MP(2,2)*mz;
      float Ppw= MP(3,0)*mx + MP(3,1)*my + MP(3,3)*mw;
      float pw = Ppw + MP(3,2)*mz;
      // normalyze
      px /= pw;
      py /= pw;
      pz /= pw;
      Vec3D v3D(px,py,pz);
      return v3D;
    };
    static void bestZoomFromList(double * mProj,double * mModel,
                                const int * viewport, const std::vector <int> * list,
                                const ParticlesData * part_data, GlobalOptions * store_options);
    static void bestZoomFromObject(double * mProj,double * mModel,
                                  const int * viewport, const ParticlesObjectVector * pov,
                                  const ParticlesData * part_data, GlobalOptions * store_options);
private:
  float mModel[16], mProj[16];
};

}

#endif
