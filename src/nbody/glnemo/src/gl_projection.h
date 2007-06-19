// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef GL_PROJECTION_H
#define GL_PROJECTION_H

class GLProjection{
public:
    GLProjection();
    ~GLProjection();
    void update(int ratio);    
private:
    void setPerspective();
    void setOrthographic();
};
#endif
// ============================================================================
