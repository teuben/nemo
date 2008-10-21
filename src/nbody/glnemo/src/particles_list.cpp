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
#include <iostream>
#include <assert.h>
#include <qmessagebox.h>
#include "particles_list.h"
#define LOCAL_DEBUG 0
#include "print_debug.h"

const std::string ParticlesList::header = "#p_index_list";
// ============================================================================
// constructor                                                                 
ParticlesList::ParticlesList():VirtualParticlesSelect()
{
  v_type    = 2;         // record VirtualParticlesSetect type
  step_part = 1;
  npart     = 0;         // raz particles number              
  index_tab = NULL;
  ni_index  = 0;
}

// ============================================================================
// destructor                                                                  
ParticlesList::~ParticlesList()
{
  //index_list.clear(); // raz index list                    
  if (index_tab) {
    delete [] index_tab;
  }
}
// ============================================================================
// ParticlesList::defaultIndexTab()                                            
// fill the index tab array with the default number of particles               
inline int ParticlesList::defaultIndexTab()
{
  if (index_tab) {
    delete [] index_tab;
  }
  index_tab = new int[npart];
  ni_index = 0;
  for (int i=0; i<npart; i+=step_part) {
    index_tab[ni_index++] = (int) index_list[i];
  }
  return 1;
}
// ============================================================================
// ParticlesList::parseSelectedString                                          
int ParticlesList::parseSelectedString(char * select_string, const int _nbody, 
				      ParticlesSelectVector * psv)
{
  if (_nbody); // remove compiler warning
  if (psv);   // remove compiler warning
  loadFile(select_string,psv);
  return 1;
}
// ============================================================================
// ParticlesList::loadFile()                                                   
int ParticlesList::loadFile(const char * select_file,ParticlesSelectVector * psv)
{
  list_file = select_file; // list file
  std::ifstream fi(list_file.c_str());

  try {  
    // try to open file
    if (! fi.is_open()) {
      
      //exit(1);
      throw(-1);
    }
    else {
      // count nparticles in the range vector
      int range_npart=npartSelected(psv,1);
      PRINT_D std::cerr << "ParticlesList::loadFile range_npart="
                        << range_npart <<"\n";
      //read header
      std::string _header;
      fi >> _header;
      if ( _header != header ) {
        throw(-2);
      }
      // read file    
      while (! fi.eof()) {
        int index;
        fi >> index;  // read index from file           
        if ( ! fi.eof() ) {
          if (index < range_npart) {
            npart++;
            index_list.push_back(index); // fill up vector
          }
        }
      } 
    }
  } catch (int n) {
      std::string message;
      switch (n) {
      case -1: 
        message = "ParticlesList::loadFile, "
                  "Failed to open List file[" +
                  list_file+"]";
        PRINT_D std::cerr << message<<"\n";         
        break;
      case -2:
        message = "ParticlesList::loadFile, "
                  " unknown type file ["
                  + list_file+"]";
        PRINT_D std::cerr << message << "\n";            
          break;
      default:
          assert(1);
     } //switch    
     error_message = message;
     throw(n); // throw back the exception
  }
  defaultIndexTab();
  return npart;
}
// ============================================================================
// ParticlesList::getIndex()                                                  
// return index of the particle                                                
int ParticlesList::getIndex(int index)
{
  //assert(index_list[index]<nbody);
  return (index_list[index]);
}
// ============================================================================

