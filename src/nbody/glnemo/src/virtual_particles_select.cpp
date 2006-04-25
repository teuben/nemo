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
#include <iostream>
#include <assert.h>

#define LOCAL_DEBUG 0
#include "print_debug.h"

#include "virtual_particles_select.h"
#include "particles_list.h"
int VirtualParticlesSelect::nb_select=0;

// ============================================================================
// Constructor                                                                 
VirtualParticlesSelect::VirtualParticlesSelect()
{
  nb_select++;               
  //setColor();
  PRINT_D std::cerr << "VirtualParticlesSelect Constructor nb_select : " 
                    << nb_select << "\n";
  npart       = -1;
  first_part  = -1;
  last_part   = -1;
  step_part   = -1;
  is_visible  = TRUE;
  v_type      = 0;
  index_tab   = NULL;
  ni_index    = 0;
}
// ============================================================================
// Copy Constructor                                                            
VirtualParticlesSelect::VirtualParticlesSelect(const VirtualParticlesSelect&m) {
  npart            = m.npart;
  first_part       = m.first_part;
  last_part        = m.last_part;
  col              = m.col;
  is_visible       = m.is_visible;
  step_part        = m.step_part;
  v_type           = m.v_type;
  list_file        = m.list_file;
  index_list       = m.index_list;
  ni_index         = m.ni_index;
  index_tab        = new int[ni_index];
  memcpy(index_tab,m.index_tab,sizeof(int)*ni_index);
};
// ============================================================================
// Copy Constructor                                                            
const VirtualParticlesSelect::VirtualParticlesSelect& VirtualParticlesSelect::operator=(const VirtualParticlesSelect&m) {
  npart            = m.npart;
  first_part       = m.first_part;
  last_part        = m.last_part;
  col              = m.col;
  is_visible       = m.is_visible;
  step_part        = m.step_part;
  v_type           = m.v_type;
  list_file        = m.list_file;
  index_list       = m.index_list;
  ni_index         = m.ni_index;
  index_tab        = new int[ni_index];
  memcpy(index_tab,m.index_tab,sizeof(int)*ni_index);  
  return *this;
};
// ============================================================================
// Destructor                                                                  
VirtualParticlesSelect::~VirtualParticlesSelect()
{
}
// ============================================================================
// VirtualParticlesSelect::setColor()                                          
void VirtualParticlesSelect::setColor()
{
  //QColor modulo_col[5] = { Qt::white, Qt::green, Qt::yellow, Qt::red, Qt::blue };
  
  int modulo_col[5][3] = { 
                            { 255, 75, 39 },
                            { 214,214, 52 },
                            { 114,214, 32 },
                            {  58, 61,214 },
                            { 214, 47,197 }
                          };
  col = QColor(modulo_col[(nb_select-1)%5][0],
               modulo_col[(nb_select-1)%5][1],
               modulo_col[(nb_select-1)%5][2]);
}
// ============================================================================
// VirtualParticlesSelect::defaultIndexTab()                                   
// fill index_tab with default value.                                          
int VirtualParticlesSelect::defaultIndexTab()
{
  std::cerr << "[VirtualParticlesSelect::defaultIndexTab()], Should not be here "
            << "\naborted....\n";
  std::exit(1); 
}
// ============================================================================
// VirtualParticlesSelect::getIndex()                                          
// return particle index at the requested index.                               
int VirtualParticlesSelect::getIndex(int index)
{
  if (index);
  std::cerr << "[VirtualParticlesSelect::getIndex()], Should not be here "
            << "\naborted....\n";
  std::exit(1);
}
// ============================================================================
// ParticlesRange::parseString()                                               
// parse 'select_string' according to the 'nemoinpi' rules.                    
// return the string at the position after the next 'coma' otherwise NULL      
char * VirtualParticlesSelect::parseString(const char * select_string, const int nbody, 
				   ParticlesSelectVector * psv)
{
  char * status;
  PRINT_D std::cerr << "In parseString...["<< select_string << "]\n";

  char * c = strchr(select_string,',');
  int sup;
  if ( c) {
    status = c+1;
    sup = c-select_string;
  } else {
    status = NULL;
    sup = strlen(select_string)+1;
  }
  char tmp[100];
  strncpy(tmp,select_string,sup);
  tmp[sup] = '\0';
  parseSelectedString(tmp,nbody,psv);
  return status;
}
// ============================================================================
// VirtualParticlesSelect::storeParticles()                                    
// fill up particles range vector                                              
template <class T> int VirtualParticlesSelect::storeParticles(ParticlesSelectVector * psv,
                                   const int nbody,const char * sel2)
{
  VirtualParticlesSelect * vps;
  const char * s = sel2;

  while (s) {
    vps = new T();
    vps->setColor();
    s=vps->parseString(s,nbody,psv);
    if (s) {
      PRINT_D std::cerr << " >>>> s sring = ["<< s << "]\n";
    }
    ParticlesSelect * ps = new ParticlesSelect();
    ps->vps = vps;
    psv->push_back(*ps);
    PRINT_D std::cerr << "In VirtualParticlesSelect::storeParticles, psv->size() = " 
                 << psv->size() << "\n";	
     //PRINT_D vps->printRange();              
    delete ps;
  }
  for (int i=0; i< (int) psv->size(); i++) {
    PRINT_D std::cerr << " - - - - - - - - - - - \n";
    PRINT_D std::cerr << i << "\n";
    PRINT_D std::cerr << (*psv)[i].vps->npart << "\n";
    PRINT_D (*psv)[i].vps->printRange();
  }
  return psv->size();
}
// ============================================================================
// VirtualParticlesSelect::storeParticlesList()                                
// fill up particles range vector                                              
int VirtualParticlesSelect::storeParticlesList(ParticlesSelectVector * psv,
                                   const int nbody,const char * sel2)
{
  VirtualParticlesSelect * vps;
  const char * s = sel2;

  while (s) {
    vps = new ParticlesList();
    vps->setColor();
    try { // try to parse
      s=vps->parseString(s,nbody,psv);
      if (s) {
        PRINT_D std::cerr << " >>>> s sring = ["<< s << "]\n";
      }
      ParticlesSelect * ps = new ParticlesSelect();
      ps->vps = vps;
      psv->push_back(*ps);
      PRINT_D std::cerr << "In VirtualParticlesSelect::storeParticles,"
                        <<   " psv->size() = " 
                        << psv->size() << "\n";
      PRINT_D std::cerr << "In VirtualParticlesSelect::storeParticlesList  :"
                        << vps->index_list.size() <<"\n";             	
      //PRINT_D vps->printRange();              
      delete ps;
    }
    catch (int n) {
        switch (n) {
        case -1: 
            break;
        case -2:
            break;	
        default:
            assert(1);
        } //switch
        error_message = vps->error_message;
        delete vps;
        throw(n); // send back errot
    } //catch
  }
  for (int i=0; i< (int) psv->size(); i++) {
    PRINT_D std::cerr << " - - - - - - - - - - - \n";
    PRINT_D std::cerr << i << "\n";
    PRINT_D std::cerr << (*psv)[i].vps->npart << "\n";
    PRINT_D (*psv)[i].vps->printRange();
  }
  return psv->size();
}
// ============================================================================
// VirtualParticlesSelect::printRange()                                        
// print particles range, usefull for debugging                                
void VirtualParticlesSelect::printRange()
{
  PRINT_D std::cerr << "Npart       = " << npart      << "\n";
  PRINT_D std::cerr << "First_part  = " << first_part << "\n";
  PRINT_D std::cerr << "Last_part   = " << last_part  << "\n";
  PRINT_D std::cerr << "Step_part   = " << step_part  << "\n";
  PRINT_D std::cerr << "v_type      = " << v_type     << "\n";
  PRINT_D std::cerr << "Is visible? = " << is_visible << "\n";
}
// ============================================================================
// VirtualParticlesSelect::npartSelected()                                     
// return the #particles according to v_type                                   
int VirtualParticlesSelect::npartSelected(ParticlesSelectVector * psv,int v_type)
{ 
  int npart=0;
  for (int i=0; i< (int) psv->size(); i++) {
    if ((*psv)[i].vps->v_type == v_type) {
      npart += (*psv)[i].vps->npart;
    }
  }
  return npart;
}
// ============================================================================
// VirtualParticlesSelect::resetIndexTab()                                     
// reset ni_index value                                                        
inline int VirtualParticlesSelect::resetIndexTab()
{
  ni_index = 0;
}
// ============================================================================
// VirtualParticlesSelect::addIndexTab()                                       
// add a particle in the index_tab array                                       
inline int VirtualParticlesSelect::addIndexTab(int index)
{
  index_tab[ni_index++] = index;
  assert(ni_index<=npart);
}
// ============================================================================
