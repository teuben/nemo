// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef FRUSTUMCULLING_H
#define FRUSTUMCULLING_H

/**
@author Jean-Charles Lambert
*/
class FrustumCulling
{
public:
    FrustumCulling();

    ~FrustumCulling();
   void getFC();
   bool isPointInside(float, float, float);
private:
   float   proj[16];
   float   modl[16];
   float   clip[16];
   float   frustum[6][4];

};

#endif
