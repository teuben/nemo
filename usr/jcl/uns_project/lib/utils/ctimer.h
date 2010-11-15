// ============================================================================
// Copyright Jean-Charles LAMBERT - 2009                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pole de l'Etoile, site de Chateau-Gombert                         
//           38, rue Frederic Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================

/*
  @author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
*/

#ifndef CTIMER_H
#define CTIMER_H

#include <ctime>
#include <sys/time.h>
#include <cstdio>
#include <unistd.h>

namespace jclut {

  class CTimer {
  public:
    
    CTimer() {
      restart();
    }
    void restart() {
      restartElapsed();
      restartCpu();
    }

    void restartElapsed() {
      struct timeval start;
      gettimeofday(&start, NULL);
      elapsed_start = start.tv_sec+(start.tv_usec/1000000.0);
    }

    void restartCpu() {
      _start_time = std::clock(); 
    }
    
    double elapsed() {
      struct timeval end;
      gettimeofday(&end, NULL);
      double elapsed_end = end.tv_sec+(end.tv_usec/1000000.0);
      return ((double) (elapsed_end-elapsed_start));
    }

    double cpu() {
      return  double(std::clock() - _start_time) / CLOCKS_PER_SEC;
    }

  private:
    double elapsed_start;
    std::clock_t _start_time;
 
  };

}

#endif
//
  
