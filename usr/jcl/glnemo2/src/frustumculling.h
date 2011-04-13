// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
@author Jean-Charles Lambert
*/
#ifndef FRUSTUMCULLING_H
#define FRUSTUMCULLING_H

namespace glnemo {
class FrustumCulling
{
public:
    FrustumCulling();

    ~FrustumCulling();
   void getFC(double *, double *);
   bool isPointInside(float, float, float);
private:
   float   proj[16];
   float   modl[16];
   float   clip[16];
   float   frustum[6][4];

};
} // end of namespace
#endif
