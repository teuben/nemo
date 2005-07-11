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
//                                                                             
// VirtualData Class implementation                                            
//                                                                             
// ============================================================================
#include <math.h> // happy gcc 3.4.x
#include "virtual_data.h"
#include "particles_range.h"

#define LOCAL_DEBUG 0
#include "print_debug.h"
  
// ============================================================================
// VirtualData::loadPos()                                                      
int VirtualData::loadPos(ParticlesSelectVector * psv)
{
  if (psv) ; // to remove compiler warning
  return 0;
}
// ============================================================================
// VirtualData::getNbody()                                                     
int VirtualData::getNbody()
{
  return 0;
}
// ============================================================================
// VirtualData::getPos()                                                       
float * VirtualData::getPos()
{
  return NULL;
}
// ============================================================================
// VirtualData::getTime()                                                      
float VirtualData::getTime()
{
  return 0.;
}
// ============================================================================
// VirtualData::isValidData()                                                  
bool VirtualData::isValidData()
{
  return FALSE;
}
// ============================================================================
// VirtualData::getCooIndexMax()                                               
int * VirtualData::getCooIndexMax()
{
  return NULL;
}
// ============================================================================
// VirtualData::getCooMax()                                                    
float * VirtualData::getCooMax()
{
  return NULL;
}
// ============================================================================
// VirtualData::uploadGlData()                                                 
void VirtualData::uploadGlData(ParticlesSelectVector * psv)
{
  if (psv) ; // to remove compiler warning
}
// ============================================================================
// VirtualData::endOfDataMessage()                                             
QString VirtualData::endOfDataMessage()
{
  return "No message in VirtualData base class";
}
// ============================================================================
// VirtualData::getDataName()                                                  
const char * VirtualData::getDataName()
{
  return NULL;
}
// ============================================================================
// VirtualData::getDataType()                                                  
const char * VirtualData::getDataType()
{
  return NULL;
}
// ============================================================================
// VirtualData::fillParticleRange()                                            
// fill up particles range vector                                              
int VirtualData::fillParticleRange(ParticlesSelectVector * psv,
                                   const int nbody,const char * sel2)
{
  VirtualParticlesSelect * vps;
  const char * s = sel2;

  while (s) {
    vps = new ParticlesRange();
    vps->setColor();
    s=vps->parseString(s,nbody,psv);
    if (s) {
      PRINT_D cerr << " >>>> s sring = ["<< s << "]\n";
    }
    ParticlesSelect * ps = new ParticlesSelect();
    ps->vps = vps;
    psv->push_back(*ps);
    PRINT_D cerr << "In VirtualData::fillParticleRange, psv->size() = " 
                 << psv->size() << "\n";	
     //PRINT_D vps->printRange();              
    delete ps;
  }
  for (int i=0; i< (int) psv->size(); i++) {
    PRINT_D cerr << " - - - - - - - - - - - \n";
    PRINT_D cerr << i << "\n";
    PRINT_D (*psv)[i].vps->npart;
    PRINT_D (*psv)[i].vps->printRange();
  }
  return psv->size();
}
// ============================================================================
// VirtualData::computeCooMax()                                                
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
// ============================================================================
