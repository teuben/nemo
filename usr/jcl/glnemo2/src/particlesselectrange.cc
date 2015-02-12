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
#include "particlesselectrange.h"
#include <nemo.h>
namespace glnemo {
// ============================================================================
// constructor                                                                 
ParticlesSelectRange::ParticlesSelectRange(const int _nbody,              // nbody in the snapshot  
                                           const std::string _select,     // particles to select
                                           const bool _keep_all,          // keep all particles ?
                                           ParticlesObjectVector * _pov)  // vector to store objects
{
  nbody    = _nbody;
  select   = _select;
  keep_all = _keep_all;
  pov      = _pov;
  parseSelect();
}

// ============================================================================
// destructor                                                                  
ParticlesSelectRange::~ParticlesSelectRange()
{
}
// ============================================================================
// parseSelect                                                                 
void ParticlesSelectRange::parseSelect()
{
  std::string current_s,next_s;
  next_s = select;
  while ((current_s=parseString(next_s)) != "") {
    std::cerr << "current string["<<current_s<<"]...\n";
    storeObject(current_s);
  }
}
// ============================================================================
// parseString                                                                 
// return the string at the position after the next 'coma' otherwise ""        
std::string ParticlesSelectRange::parseString(std::string & next_string)
{
  std::string return_string;
  std::string::size_type coma=next_string.find(",",0);  // try to find "'"
  if (coma != std::string::npos) { // found coma
    return_string = next_string.substr(0,coma);
    next_string   = next_string.substr(coma+1,next_string.length());
  } else {                         // not found
    return_string = next_string;
    next_string = "";
  }
  return return_string;
}
// ============================================================================
// storeObject
void ParticlesSelectRange::storeObject(const std::string current)
{
  int * indexes = new int[nbody]; // store the indexes
  int npart,first,last;
  // find out particles in the current string
  if ( current != "all" ) {
    // gonna use nemoinpi to figure out how many particles there
    // are in current string
    npart = nemoinpi((char *) current.c_str(),indexes,nbody);
    if (npart <=0 ) {
      std::cerr << "WARNING : misformated selected string <" << current
                << "," << select << ">\n";
      throw(-40); // misformated selected string
    }
#if 0
    int nb=0; // #body selected
    for (int i=0;i<npart; i++)
        if (indexes[i] < nbody) nb++;
    npart = nb;
#endif
  }
  else { // current==all
    for (int i=0;i<nbody;i++)
      indexes[i]=i;
    npart=nbody;
  }
  // store particles in vector
  if (keep_all) {     // we keep all the particles
    first = indexes[0];
    last  = indexes[npart-1];
  }
  else {            // we keep the selected particles only
    if (pov->size() == 0) {
      first = 0;
      last  = npart-1;
    }
    else {
      first = (*pov)[pov->size()-1].last + 1;
      last  = first + npart-1;
    }
  }
  ParticlesObject * po = new ParticlesObject(); // new object
  po->buildIndexList( npart,first,last);        // object's particles indexes
  pov->push_back(*po);                          // store in vector
  delete po;
  delete [] indexes; // useless anymore
}
}
