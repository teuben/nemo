// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
//                                                                             
// VirtualData Class implementation                                            
//                                                                             
// ============================================================================
#include <math.h> // happy gcc 3.4.x
#include "virtual_data.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"
  
// ============================================================================
//
int VirtualData::loadPos(ParticlesRangeVector * prv)
{
  if (prv) ; // to remove compiler warning
  return 0;
}
// ============================================================================
//
int VirtualData::getNbody()
{
  return 0;
}
// ============================================================================
//
float * VirtualData::getPos()
{
  return NULL;
}
// ============================================================================
//
float  VirtualData::getTime()
{
  return 0.;
}
// ============================================================================
//
bool VirtualData::isValidData()
{
  return FALSE;
}
// ============================================================================
//
int * VirtualData::getCooIndexMax()
{
  return NULL;
}
// ============================================================================
//
float * VirtualData::getCooMax()
{
  return NULL;
}
// ============================================================================
//
void VirtualData::uploadGlData(ParticlesRangeVector * prv)
{
  if (prv) ; // to remove compiler warning
}
// ============================================================================
//
QString VirtualData::endOfDataMessage()
{
  return "No message in VirtualData base class";
}
// ============================================================================
//
const char * VirtualData::getDataName()
{
  return NULL;
}
// ============================================================================
//
const char * VirtualData::getDataType()
{
  return NULL;
}
// ============================================================================
//
int VirtualData::fillParticleRange(ParticlesRangeVector * prv,
                                   const int nbody,const char * sel2)
{
  ParticlesRange * pr;
  const char * s = sel2;

  while (s) {
    pr = new ParticlesRange();
    s=pr->parseString(s,nbody,prv);
    if (s) 
      PRINT_D cerr << " >>>> s sring = ["<< s << "]\n";
    prv->push_back(*pr);
    PRINT_D cerr << "In globwin, prv->size() = " << prv->size() << "\n";	  
    delete pr;
  }
  for (int i=0; i< (int) prv->size(); i++) {
    PRINT_D cerr << " - - - - - - - - - - - \n";
    PRINT_D cerr << i << "\n";
    PRINT_D (*prv)[i].printRange();
  }
  return prv->size();
}
// ----------------------------------------------------------------------------
//  compute extremum coordinates
void VirtualData::computeCooMax()
{
  coo_max[0]= fabs(pos[0]);
  i_max[0]  = 0;
  coo_max[1]= fabs(pos[1]);
  i_max[1]  = 0;
  coo_max[2]= fabs(pos[2]);
  i_max[2]  = 0;
  
  for (int i=0;i<*nbody;i++) {
    if (fabs(pos[i*3  ]) > coo_max[0]) {
      coo_max[0] = fabs(pos[i*3  ]);
      i_max[0]   = i;
    }
    if (fabs(pos[i*3+1]) > coo_max[1]) {
      coo_max[1] = fabs(pos[i*3+1]);
      i_max[1]   = i;
    }
    if (fabs(pos[i*3+2]) > coo_max[2]) {
      coo_max[2] = fabs(pos[i*3+2]);
      i_max[2]   = i;
    }
  }
  PRINT_D cerr << "Max coordinates \n";
  PRINT_D cerr << coo_max[0] << " " << coo_max[1] << " " << coo_max[2] << "\n";
  PRINT_D cerr << i_max[0] << " " << i_max[1] << " " << i_max[2] << "\n";
}
// 
//
