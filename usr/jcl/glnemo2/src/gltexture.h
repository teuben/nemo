// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
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
#ifndef GLNEMOGLTEXTURE_H
#define GLNEMOGLTEXTURE_H
#include <QString>
#include <QGLWidget>
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
namespace glnemo {

class GLTexture;

typedef std::vector <GLTexture> GLTextureVector;

class GLTexture{
public:
    GLTexture();
    ~GLTexture();
    bool load(QString _texture_name, QString _path, bool embeded=true);
    static  int loadTextureVector(GLTextureVector & gtv);
    int glBindTexture();
    float U() { return u_max;};
    float V() { return v_max;};
    static QString TEXTURE[][2];
    QString getName() { return texture_name; }
    void deleteTexture() { if (texture) glDeleteTextures(1,&texture); };
private:
  float u_max, v_max;
  QString texture_name, path;
  GLuint texture; // storage for the texture
  bool status;
  float evalHermite(float pA, float pB, float vA, float vB, float u);
  unsigned char* createGaussianMap(int N);
  void createGaussian(int resolution);
};

}

#endif
