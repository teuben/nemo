// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
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
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOGLOBJECTPARTICLES_H
#define GLNEMOGLOBJECTPARTICLES_H

#include "globject.h"
#include "gltexture.h"
#include <QObject>

//#include "particlesdata.h"
namespace glnemo {

class GLObject;
class GLObjectParticles;
class ParticlesData;
class ParticlesObject;
class GlobalOptions;
class GLTexture;
class GLObjectIndexTab;

typedef std::vector <GLObjectIndexTab> GLObjectIndexTabVector;

class GLObjectIndexTab {
  public:
    GLObjectIndexTab() {};
    const GLObjectIndexTab & operator=(const GLObjectIndexTab& m) {
      index   = m.index;
      value   = m.value;
      i_point = m.i_point;
      return *this;
    }
    static bool compareLow(const  GLObjectIndexTab &a ,const GLObjectIndexTab &b) {
      return a.value < b.value;
    };
    static bool compareHigh(const  GLObjectIndexTab &a ,const GLObjectIndexTab &b) {
      return a.value > b.value;
    };

    int index, i_point;
    float value;
};

typedef std::vector <GLObjectParticles> GLObjectParticlesVector;

class GLObjectParticles : public GLObject {
  
  public:
    GLObjectParticles(GLTextureVector *);
    GLObjectParticles(const ParticlesData   *,
                      ParticlesObject *,
                      const GlobalOptions   *,
		      GLTextureVector *);
    ~GLObjectParticles();
    void update(const ParticlesData   *,
           ParticlesObject *,
           const GlobalOptions   *);
    void updateVel();
    void updateVbo();
    void updateColorVbo();
    const ParticlesData * getPartData() const { return part_data; };
    const ParticlesObject * getPartObj() const { return (const_cast <ParticlesObject *>(po)); };
    void buildDisplayList();
    void buildVelDisplayList();
    void buildOrbitsDisplayList();
    void buildVboPos();
    void buildVboColor();
    void buildVboColorTempGasSorted();
    void buildVboColorGasGasSorted();
    void buildVboSize();
    void buildVboSize2();
    void display(const double * mModel, int);
    void setTexture(QString);
    void setTexture(const int);
    void setTexture();
    void toto();
    void checkVboAllocation(const int sizebuf);

  private:
    // Data
    const ParticlesData * part_data;
    ParticlesObject * po;
    const GlobalOptions * go;
    GLuint vel_dp_list, orb_dp_list;
    GLTexture * texture;
    GLTextureVector * gtv;
    // method
    void displaySprites(const double *mModel);
    void displayVboSprites(int,const bool);
    void displayVboPoints();
    void sortByDepth();
    void sortByDensity();
    // vbo
    GLuint vbo_pos, vbo_color , vbo_size, vbo_index, vbo_index2;
    int nvert_pos;
    // Rho
    GLObjectIndexTabVector rho,zdepth;
    GLuint * indexes_sorted, nind_sorted;
    //
    void initShader();
    static int compareZ(const void * a, const void * b);
};

}

#endif
