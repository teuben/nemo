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
#ifndef FRAME_DATA_H
#define FRAME_DATA_H

#include <qdatetime.h>
#include <vector>

#include "global_options.h"

class FrameData;
typedef std::vector <FrameData> FrameDataVector;

class FrameData{
public:
    FrameData();

    ~FrameData();
    float time;        // real frame time                   
    int elapsed;       // elapsed time for the current frame
    int cumul_elapsed; // cumul elapsed from the beginning  
    GlobalOptions store_options;
};

#endif
