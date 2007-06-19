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
#ifndef PARTICLES_DATA_H
#define PARTICLES_DATA_H

/**
@author Jean-Charles Lambert
*/
class ParticlesData{
public:
    ParticlesData();
    const ParticlesData& operator=(const ParticlesData& m);
    
    int allocVar();
    int allocTree();
    int  * nbody, i_max[3] ;
    float * pos,
          * vel, * vel_norm,
          coo_max[3], 
          * timu;
    int * tree_depth;   // array to store tree_depth for each particles
    int tree_size_max;
    int * nemobits;      
    ~ParticlesData();
     char * mallocate(char *, int, bool force=false);
     void computeVelNorm();
};

#endif
