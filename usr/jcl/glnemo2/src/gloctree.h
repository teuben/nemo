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
// getOctant(x,y,z) ??                                                         
//                                                                             
//      b_____________a                                                        
//     /|            /|                                                        
//    / |           / |                                                        
//   f--|----------e  |                                                        
//   |  c__________|__d                                                        
//   | /           | /                                                         
//   |/            |/                                                          
//   g-------------h                                                           
//                                                                             
// a = 111 = 7                                                                 
// b = 011 = 3                                                                 
// c = 001 = 1                                                                 
// d = 101 = 5                                                                 
// e = 110 = 6                                                                 
// f = 010 = 2                                                                 
// g = 000 = 0                                                                 
// h = 100 = 4                                                                 
// ============================================================================
#ifndef GLNEMOGLOCTREE_H
#define GLNEMOGLOCTREE_H
/**
        @author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
 */
#include <globject.h>
#include "particlesobject.h"

namespace glnemo {
class GlobalOptions;
class ParticlesData;

class Node {
  public:
    Node(int index_=-1, int obj_=-1, Node * _levelup=NULL);
    ~Node();
    int npart;
    int type;
    int index;
    int obj;
    Node * levelup;
    Node * node[8];
};

class GLOctree : public GLObject {
  Q_OBJECT
public:
    GLOctree(GlobalOptions * _options);
    ~GLOctree();
    void update(ParticlesData   * ,
                ParticlesObjectVector *);
    void update();
    void display();
  public slots:
    void buildDisplayList();
    void computeDensity();
    
  private:
    Node * root;
    int level_max;     // level tree
    int level_requested;
    int size_max;      // max square size in power of 2
    bool new_data;     // true if new
    int nbody_keeped;  // #bodies keep after tree walk
    float rmid[3];
    ParticlesData * p_data;
    GlobalOptions * store_options;
    ParticlesObjectVector  * pov;
    bool first;
    void computeSizeMax();
    int build();
    void expandTree(float );
    void insertParticle(int, Node **, float *, int,  int, Node *);
    int hackTreeDL(Node *, float *, int );
    int hackTreeDensity(Node *, const int);
    int hackTreeCountPart(Node *, const int);
    void init();

};

}

#endif
