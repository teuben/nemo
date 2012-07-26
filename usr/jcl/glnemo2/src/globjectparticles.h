// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOGLOBJECTPARTICLES_H
#define GLNEMOGLOBJECTPARTICLES_H
#include "cshader.h"
#include "globject.h"
#include "gltexture.h"
#include <QObject>
#include <iostream>
#include <vector>

//#include "particlesdata.h"
namespace glnemo {

class GLObject;
class GLObjectParticles;
class ParticlesData;
class PhysicalData;
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
		      GLTextureVector *, CShader *);
    ~GLObjectParticles();
    void update(const ParticlesData   *,
                ParticlesObject *,
                const GlobalOptions   *,
                const bool update_obj=true);
    void updateVel();
    void updateVbo();
    void updateColorVbo();
    void updateBoundaryPhys();
    const ParticlesData * getPartData() const { return part_data; }
    const ParticlesObject * getPartObj() const { return (const_cast <ParticlesObject *>(po)); }
    void buildDisplayList();
    void buildVelDisplayList();
    void buildOrbitsDisplayList();
    void buildVboPos();
    void buildVboHsml();
    void buildVboPhysData();
    void display(const double * mModel, int);
    void setTexture(QString);
    void setTexture(const int);
    void setTexture();
    void checkVboAllocation(const int sizebuf);
    void updateColormap();
    
  private:
    // shader
    CShader * shader;
    
    // manage min/max index for the physical quantity selected
    int min_index, max_index;
    //static const int nhisto;// #entries in index_histo
    std::vector <int> index_histo;// index_histo[100]; // store first part's index in the percentage
    // Data
    const ParticlesData * part_data;
    ParticlesObject * po;
    const GlobalOptions * go;
    GLuint vel_dp_list, orb_dp_list;
    GLTexture * texture;
    GLTextureVector * gtv;
    const PhysicalData * phys_select;
    int phys_select_id;
    // local color map
    std::vector <float> cmap; 
 
    // method
    void displaySprites(const double *mModel);
    void displayVboShader(const int,const bool use_point=false);
    void sortByDepth();
    void sortByDensity();
    void selectParticles();
    void buildIndexHisto();    
    void sendShaderColor(const int, const bool use_point);

    // vbo
    GLuint vbo_pos, vbo_color , vbo_size, vbo_index, vbo_index2, vbo_data;
    int nvert_pos;
    // Rho
    GLObjectIndexTabVector vindex_sel,phys_itv,rho_itv;
    GLuint * indexes_sorted, nind_sorted;
    //
    bool hasPhysic; // Does object has physic value
    void initShader();
    static int compareZ(const void * a, const void * b);
    void checkGlError(std::string s) {
      int err = glGetError();
      if (err)
        std::cerr << "!!!! OpenGL error["<<s<<"] error = "<<err<<"\n";
    }
};

}

#endif
