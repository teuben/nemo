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
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
 */
#ifndef GLNEMOGLOBJECTPARTICLES_H
#define GLNEMOGLOBJECTPARTICLES_H
#include "cshader.h"
#include "globject.h"
#include "gltexture.h"
#include <QObject>
#include <iostream>
#include <vector>
#include "particlesdata.h"

//#include "particlesdata.h"
namespace glnemo {

  class PhysObject;

  typedef std::vector <PhysObject> PhysObjectVector;
  // ----------------------------------------------------------------------------
  // class ObjectPhys to store per Physical Object properties
  class PhysObject {
    public:
      PhysObject(const enum PhysicalData::PHYS _type,
                 const float * _data,const int _npart, const int * _index_tab):type(_type){}
      ~PhysObject() {}
      bool hasPhyisc() {
        return valid;
      }
      // real Physical value
      void setMinPhys(const float _v) { min_phys = _v; }
      void setMaxPhys(const float _v) { max_phys = _v; }
      float getMinPhys()  const { return min_phys; }
      float getMaxPhys()  const { return max_phys; }
      // Percentage Physical Value
      void setMinPercenPhys(const int _v) { min_percen_phys = _v; }
      void setMaxPercenPhys(const int _v) { max_percen_phys = _v; }
      int getMinPercenPhys()  const { return min_percen_phys; }
      int getMaxPercenPhys()  const { return max_percen_phys; }

      std::vector<int> * getSortedIndexes() { return &sorted_indexes; }

    private:
      const enum PhysicalData::PHYS type;
      int npart;
      float * data; // pointer to physical data
      int * index_tab(); // index of object particles
      std::vector<int> sorted_indexes; // store indes sorted

      // data for object dialog box
      int min_percen_phys, max_percen_phys;
      float min_phys, max_phys;
      // local color map
      std::vector <float> cmap; // store color map
      bool valid; // true if has physic

  };

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
              GLTextureVector *, CShader *, CShader *);
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
    void buildVboVelFactor();
    void buildVboHsml();
    void buildVboPhysData();
    //void buildVboVel();
    void display(const double * mModel, int);
    void setTexture(QString);
    void setTexture(const int);
    void setTexture();
    void checkVboAllocation(const int sizebuf);
    void updateColormap();
    
  private:
    // shader
    CShader * shader, * vel_shader;
    
    // manage min/max index for the physical quantity selected
    int min_index, max_index;
    //static const int nhisto;// #entries in index_histo
    std::vector <int> index_histo;// index_histo[100]; // store first part's index in the percentage
    // Data
    const ParticlesData * part_data;
    ParticlesObject * po; // address of the selected object
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
    void displayVboVelShader330(const int,const bool use_point=false);
    void displayVboVelShader130(const int,const bool use_point=false);
    void sortByDepth();
    void sortByDensity();
    void selectParticles();
    void buildIndexHisto();    
    void sendShaderColor(const int, const bool use_point);

    // vbo
    GLuint vbo_pos, vbo_color , vbo_size, vbo_index, vbo_index2, vbo_data, vbo_vel, vbo_vel_X2;
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
