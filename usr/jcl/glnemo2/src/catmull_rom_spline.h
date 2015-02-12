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

#ifndef CATMULL_ROM_SPLINE_H
#define CATMULL_ROM_SPLINE_H

#include "vec3d.h"
#include <vector>

namespace glnemo {
class CRSpline
{
public:

  // Constructors and destructor
  CRSpline();
  CRSpline(const CRSpline&);
  ~CRSpline();

  // Operations
  void AddSplinePoint(const Vec3D& v);
  Vec3D GetInterpolatedSplinePoint(float t);   // t = 0...1; 0=vp[0] ... 1=vp[max]
  int GetNumPoints();
  Vec3D& GetNthPoint(int n);
  
  // Static method for computing the Catmull-Rom parametric equation
  // given a time (t) and a vector quadruple (p1,p2,p3,p4).
  static Vec3D Eq(float t, const Vec3D& p1, const Vec3D& p2, const Vec3D& p3, const Vec3D& p4);

  // Clear ctrl points
  void clearCPoints() { vp.clear();}
private:
  std::vector<Vec3D> vp;
  float delta_t;
};
}
#endif
