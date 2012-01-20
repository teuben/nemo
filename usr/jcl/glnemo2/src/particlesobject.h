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
#ifndef GLNEMOPARTICLESOBJECT_H
#define GLNEMOPARTICLESOBJECT_H
#include <vector>
#include <iostream>
#include <QColor>
#include "orbits.h"

namespace glnemo {
class ParticlesObject;

typedef std::vector <ParticlesObject> ParticlesObjectVector; 
typedef std::vector <OrbitsList> OrbitsVector;

// ----------------------------------------------------------------------------
// class ParticlesObject describes the proprerties of the object:              
// #particles, type of particles, color, indexes.                              
class ParticlesObject{
  public:
    enum ObjFrom { Range, File, Select, ntype };
    //ParticlesObject();
    ParticlesObject(const ObjFrom _of=Range,
                    const std::string _name="none");
    ParticlesObject(const int, const int, const int, const int _step=1,
                    const ObjFrom _of=Range,
                    const std::string _name="none");

    ParticlesObject(const ParticlesObject&);
    const ParticlesObject& operator=(const ParticlesObject&);
    void copyProperties(const ParticlesObject&);
    ~ParticlesObject();
    static bool compareFirst(const ParticlesObject &a ,const ParticlesObject &b) {
         return a.first < b.first;
    };
    static bool comparePos(const ParticlesObject &a ,const ParticlesObject &b) {
         return a.pos < b.pos;
    };
    static void copyVVkeepProperties(ParticlesObjectVector&,ParticlesObjectVector&, const int nbody);
    static void backupVVProperties(ParticlesObjectVector& src,ParticlesObjectVector& dest, const int nsel);
    static void clearOrbitsVectorPOV(ParticlesObjectVector&);   
    static void initOrbitsVectorPOV(ParticlesObjectVector&);   
    void buildIndexList(const int, const int, const int, const int _step=1);
    void buildIndexList();
    void buildIndexList(std::vector<int> &);
    int resizeRange(const int, const int, int&);
    static int nobj; // object's index           
    int npart;       // #particles in the object
    int first;       // index of the first particle
    int last;        // index of the last particle
    int step;        // incremental step between particles.
    int * index_tab; // particles's indexes
    //
    // Orbits stuff
    //
    OrbitsVector ov;
    void addOrbits(const ParticlesData * p_data);
    //
    // Object managing method
    //
    bool isVisible()   const  { return visible;}
    void setVisible(const bool _v)  { visible = _v;  }
    const QColor getColor()  const  { return color;  }
    void setColor(const QColor _v)  { color = _v;    }
    // particles
    bool isPartEnable() const { return part;  }
    void setPart(const bool _v)     { part = _v; }
    void setPartSize(const float _v){ part_size=_v;  }
    void setPartAlpha(const int _v) { part_alpha=_v; }
    float getPartSize()       { return part_size;}
    int getPartAlpha()        { return part_alpha;}
    // gas
    bool isGazEnable() const  { return gaz;    }
    void setGaz(const bool _v)      { gaz = _v; }
    void setGazSize(const float _v) { gaz_size=_v;  }
    void setGazAlpha(const int _v)  { gaz_alpha=_v;  }
    float getGazSize()        { return gaz_size;}
    int getGazAlpha()        { return gaz_alpha; }
    bool isGazRotate()         { return gaz_rotate; }
    float getGazSizeMax()    { return gaz_size_max;}
    void setGazSizeMax(const float _v) { gaz_size_max = _v;}
    void setGazRotate( const bool _b) { gaz_rotate = _b; }
    void setTextureIndex(const int _tex) { texture_index = _tex;}
    int getTextureIndex()     { return texture_index; }
    void setGazGlsl(const bool _v)  { gaz_glsl=_v;}
    bool  isGazGlsl()         { return gaz_glsl;}
    
    // velocity
    bool isVelEnable() const  { return vel;    }
    void setVel(const bool _v)      { vel = _v; }
    void setVelSize(const float _v)   { vel_size=_v;   }
    void setVelAlpha(const int _v)  { vel_alpha=_v;  }
    float getVelSize()         { return vel_size;}
    int getVelAlpha()        { return vel_alpha; }
    void setVelFactor(const float _v) { vel_factor =_v;  }
    float getVelFactor()        {return vel_factor;}
    void setVelSizeMax(const float _v){ vel_size_max =_v; }
    float getVelSizeMax()       { return vel_size_max;}
    // orbits
    bool isOrbitsEnable() const  { return orbits; }
    bool isOrbitsRecording() const  { return o_record; }
    void setOrbits(const bool _o)      { orbits=_o;}
    void setOrbitsRecord(const bool _o){ o_record=_o;}
    void setOrbitsMax(const int _m)    { orbits_max = _m;}
    void setOrbitsHistory(const int _h)    { orbits_history = _h;}
    void setOrbitsAnimate(const bool _a) { orbits_animate = _a;}
    int getOrbitsMax()           { return orbits_max; }
    int getOrbitsHistory()       { return orbits_history;}
    bool getOrbitsAnimate()      { return orbits_animate;}
    ObjFrom selectFrom()        { return obj_from;}
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
    bool rhoSorted() const { return rho_sorted;}
    void setRhoSorted(const  bool _v) {
      rho_sorted=_v;
    }
    void setPhysic(const bool _v) { has_physic = _v;}
    bool  hasPhysic() const { return has_physic;}
  private:

    void copyDataObject(const ParticlesObject&, const bool garbage=false);
    QColor color;    // object color                    
    bool visible;    // TRUE if object is visible       
    // part
    bool part;       // TRUE to display particles
    float part_size;
    int   part_alpha;
    // gas
    bool gaz;        // TRUE to display gaz effect
    float gaz_size;
    int gaz_alpha;
    int texture_index;
    bool gaz_rotate;
    float gaz_size_max;
    bool gaz_glsl;
    // velocity
    bool vel;        // TRUE to display velocity vectors
    float vel_size, vel_size_max;
    int vel_alpha;
    float vel_factor;
    // orbits
    bool orbits;
    bool o_record;
    int orbits_max;
    int orbits_history;
    bool orbits_animate;
    // physicial value
    bool has_physic;
    float min_phys;
    float max_phys;
    int min_percen_phys;
    int max_percen_phys;
    bool rho_sorted;
    void setColor();
    std::string obj_name;  // or file name
    ObjFrom obj_from;
    void init(const ObjFrom,const std::string);
    bool freed;
    int pos;
    static long long int cpt;
};

}

#endif
