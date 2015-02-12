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
#include "gloctree.h"
#include "globaloptions.h"
#include "particlesdata.h"
#include "particlesobject.h"
#include <math.h>
#include <assert.h>
#define MAX(A,B) ((A)>(B)?(A):(B))
#define PRINT_D if (0)
namespace glnemo {

// ============================================================================
// Constructor                                                                 
Node::Node(int index_, int obj_, Node * _levelup)
{
  index = index_;
  if (index<0) {
    type = 0;  // node
    obj  = -1;
  }
  else {
    type=1;    // leaf
    obj=obj_;
  }
  levelup = _levelup;
  npart   = 0;
  for (int i=0; i<8; i++) {
    node[i]=NULL;
  }
}
// ============================================================================
// Destructor                                                                  
Node::~Node()
{
  //std::cerr << p << " ";
  for (int i=0; i<8; i++) {
    delete node[i];
  }
}
// ============================================================================
// Constructor                                                                 
GLOctree::GLOctree(GlobalOptions * _options):GLObject()
{
  // copy global options
  store_options = _options;
  root=NULL;
  pov=NULL;
  dplist_index = glGenLists( 1 );    // get a new display list index
  setColor(Qt::green);
  setActivate(store_options->octree_display);
  init();
  first=false;
}
// ============================================================================
// Destructor                                                                  
GLOctree::~GLOctree()
{
  if (root) {
    delete root;
  }
}
// ============================================================================
// Display                                                                     
void GLOctree::display()
{
#if 1  
  glEnable (GL_LINE_SMOOTH);
  
  //
  //glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
#endif
  glLineWidth (0.7);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  GLObject::display();
  glDisable(GL_BLEND);
}
// ============================================================================
// Constructor                                                                 
void GLOctree::init()
{
   // init
  level_max = 0;
  size_max  = 1;
  rmid[0]=rmid[1]=rmid[2]=0.0;
  if (root) {
    delete root;
  }
  root=NULL;
}
// ============================================================================
// GLOctree::update()                                                          
// update tree with new Particles object                                       
void GLOctree::update(ParticlesData * _p_data ,
                      ParticlesObjectVector  * pov_)
{
  first=true;
  // copy pointers
  p_data = _p_data;
  pov   = pov_;
  new_data = true;

  // build Display list
  //!!buildDisplayList();
  build();
  computeDensity();
}
// ============================================================================
// GLOctree::update()                                                          
// update tree according store_options conditions                              
void GLOctree::update()
{
#if 0
  //if (!first) return;
  setActivate(store_options->octree_display);
  // restore default object particles index
  for (int obj=0; obj< (int ) pov->size(); obj++ ) {
    if ((*pov)[obj].vps->is_visible) {
      (*pov)[obj].vps->defaultIndexTab();
    }
  }
  if ( pov && ! store_options->octree_enable) {
    setActivate(false); // do not display tree
  } else {
    buildDisplayList();
  }
#endif
}
// ============================================================================
// GLOctree::build()                                                           
// build octree                                                                
int GLOctree::build()
{
  if (pov && store_options->octree_enable) {
    new_data=false;
    init();
    computeSizeMax();
    p_data->tree_size_max = size_max;
    
    PRINT_D std::cerr << "Start build\n";
    
    //for (int obj=0; obj< (int ) pov->size(); obj++ ) {
    for (ParticlesObjectVector::iterator po=pov->begin(); po != pov->end(); po++) {
      if (po->isVisible()) {
        //for (int i=0; i < vps->npart; i+=vps->step_part) {
        for (int i=0; i < po->npart; i+=po->step) {
          int index=po->index_tab[i];
          assert(index<*p_data->nbody);
          insertParticle(index,&root, rmid, 1, po-pov->begin(),NULL);
        }
      }
    }
    hackTreeCountPart(root,1);
    PRINT_D std::cerr << "stop build Level max= " << level_max << "\n";
    return 1;
  }
  return 0;
}
// ============================================================================
// hackTreeCountPart()                                                           
int GLOctree::hackTreeCountPart(Node * tree, const int level)
{
  if ( !tree  ) { // free node
    return 0;
  }
  else
    if (tree->type == 1) { // a leaf
      return 1;
    }
    else { // a node
      for (int i=0; i<8; i++) {
        tree->npart+=hackTreeCountPart(tree->node[i],level+1);
      }
      return tree->npart;
    }
}

// ============================================================================
// GLOctree::insertParticle()                                                  
// insert a new particle in the tree                                           
void GLOctree::insertParticle(int index, Node ** node, float * rmid, int level, int obj, Node * _levelup)
{
  int bit;
  int octant=0;
  float new_rmid[3];
  //assert(level < 30);
  if (level > 30) {
    //std::cerr << "Skip particle with index [" << index << "]\n";
    return;
  }
  level_max = MAX(level,level_max);
  if ( ! (*node) ) {
    (*node) = new Node();
  }
  (*node)->levelup = _levelup; // keep track of the parent

  // look for position in the cube (octant)
  for (int i=0; i<3; i++) {
    float comp=p_data->pos[index*3+i];
    if ((comp-rmid[i]) < 0) {
      bit=0;
      new_rmid[i]=((float ) (-size_max)/(float) (1<<level))+rmid[i];
    } else {
      bit=1;
      new_rmid[i]=((float ) (size_max)/(float ) (1<<level))+rmid[i];
    }
    octant+=(bit << (2-i));
  }
  // check if octant is free
  if (((*node)->node[octant])==NULL) { // yes it's free !
    Node * new_node=new Node(index,obj); // insert a particle (leaf)
    (*node)->node[octant] = new_node;
    (*node)->node[octant]->levelup = *node;
  }
  else { // octant is not free
    if ((*node)->node[octant]->type ==0) { // it's a node
      insertParticle(index,&((*node)->node[octant]),new_rmid,level+1,obj,*node);
    }
    else { // it's a leaf
      // we must insert particles's leaf itself in a new octant
      int save_index=(*node)->node[octant]->index;
      //std::cerr << "save index="<< save_index <<"\n";
      (*node)->node[octant]->type=0; // leaf become node
      int old_obj = (*node)->node[octant]->obj; // must use old obj !!
      insertParticle(save_index,&((*node)->node[octant]),new_rmid,level+1,old_obj,*node);
      // we must now insert our particle
      insertParticle(index,&((*node)->node[octant]),new_rmid,level+1,obj,*node);
    }
  }
}
// ============================================================================
// GLOctree::computeDensity()
void GLOctree::computeDensity()
{ 
  // display list
  if (pov && store_options->octree_enable) {
    if (!p_data->rho) {
      p_data->rho = new PhysicalData(PhysicalData::rho,*p_data->nbody);
      for (int i=0; i<*p_data->nbody; i++) {
        p_data->rho->data[i] = -1.;
      }
    }
    if (!p_data->rneib) {
      p_data->rneib = new PhysicalData(PhysicalData::neib,*p_data->nbody);
    }
    hackTreeDensity(root,1);
    PRINT_D std::cerr << "Stop build display List #bodies kept ["<<nbody_keeped<<"]\n";
    std::cerr << "Density estimation completed...\n";
    p_data->rho->computeMinMax();
  }
}
// ============================================================================
// hackTreeDensity()                                                           
int GLOctree::hackTreeDensity(Node * tree, const int level)
{
  if ( !tree  ) { // free node
    return 0;
  }
  else
    if (tree->type == 1) { // a leaf
    float density;
    if (tree->levelup->npart > 4) {
      float cell_size = (float ) size_max/(1<<(level));
      density = tree->levelup->npart/(cell_size * cell_size * cell_size);

    }
    else {
      density = 0.000000001;
    }
      //int index=(*pov)[tree->obj].index_tab[tree->index];
      p_data->rho->data[tree->index] = density;
      int level2;
      if (level-3>0) level2=level-3;
      else if (level-2>0) level2=level-2;
        else if (level-1>0) level2=level-1;
          else level2 = level;
      p_data->rneib->data[tree->index] = (float ) size_max/(1<<(level2));
      return 1;
   }
   else { // a node
     //std::cerr << "level="<<level<<" npart="<<tree->npart<<"\n";
     for (int i=0; i<8; i++) {
       hackTreeDensity(tree->node[i],level+1);
     }
     return 1;
   }
}
// ============================================================================
// GLOctree::buildDisplayList()
void GLOctree::buildDisplayList()
{
  // display list
  if (pov && store_options->octree_enable) {
    for (ParticlesObjectVector::iterator po=pov->begin(); po != pov->end(); po++) {
      if (po->isVisible()) {
        //vps->resetIndexTab();
        //aa[obj]=0;
      }
    }
    nbody_keeped=0;
    PRINT_D std::cerr << "Start build display List\n";
    glNewList( dplist_index, GL_COMPILE );
    hackTreeDL(root,rmid,1);
    glEndList();
    PRINT_D std::cerr << "Stop build display List #bodies kept ["<<nbody_keeped<<"]\n";
    //std::cerr << "aa0 = " << aa[0] << "   aa1 = " << aa[1] << "\n";
  }
}
// ============================================================================
// GLOctree::hackTreeDL()
int GLOctree::hackTreeDL(Node * tree, float * rmid, int level)
{
  float vv[8][3]; // 8 vertex to store
  //std::cerr << "hi\n";
  float new_rmid[3];
  if ( !tree  ) { // free node
    return 0;
  }
  else if (tree->type == 1) { // a leaf
    if (level >= store_options->octree_level) { // selected particles
      nbody_keeped++; // one more particles to save
      //(*pov)[tree->obj].vps->addIndexTab(tree->index);
    }
    return 1;
  }
  else { // a node
    
    
    // Build
    float square_minus1=(float ) size_max/(1<<(level-1));
    float square=(float ) size_max/(1<<(level));
    int fac[2] = { -1,1 };
    for (int i=0;i<2;i++) {
      for (int j=0;j<2;j++) {
        for (int k=0;k<2;k++) {
          int octant=(i<<2)+(j<<1)+k;
          if (store_options->octree_display) {
            // new middle coordinates
            new_rmid[0]=fac[i]*square+rmid[0];
            new_rmid[1]=fac[j]*square+rmid[1];
            new_rmid[2]=fac[k]*square+rmid[2];
            // vertex's octant coordinates
            vv[octant][0]=fac[i]*square_minus1+rmid[0];
            vv[octant][1]=fac[j]*square_minus1+rmid[1];
            vv[octant][2]=fac[k]*square_minus1+rmid[2];
            //
          }
          assert(octant<8);
          hackTreeDL(tree->node[octant],new_rmid,level+1);
        }
      }
    }
    // plot the 6 squares  North - South - Est - West - Top - Bottom
    if (level >= store_options->octree_level-1) {
      if (store_options->octree_display) {
        int square_index[6][4] = { {7,3,1,5}, {6,2,0,4},
          { 3,2,0,1}, {7,6,4,5}, {7,3,2,6}, {1,5,4,0}  };
          for (int i=0; i < 6; i++) {
            glBegin(GL_LINE_LOOP);
            for (int j=0; j<4; j++) {
            //std::cerr << " \n";
              glVertex3f(vv[square_index[i][j]][0],
                         vv[square_index[i][j]][1],
                         vv[square_index[i][j]][2]);
            }
            glEnd();
          }
      }
    }
    return 1;
  }
}
// ============================================================================
// GLOctree::computeSizeMax()                                                  
//  compute extremum coordinates                                               
void GLOctree::computeSizeMax()
{
  bool visible=false;
  float coo_max[3];
  int i_max[3];
  for (ParticlesObjectVector::iterator po=pov->begin(); po != pov->end(); po++) {
    if (po->isVisible()) {
      visible = true;
      coo_max[0]= fabs(p_data->pos[po->index_tab[0]]);
      i_max[0]  = 0;
      coo_max[1]= fabs(p_data->pos[po->index_tab[0]]);
      i_max[1]  = 0;
      coo_max[2]= fabs(p_data->pos[po->index_tab[0]]);
      i_max[2]  = 0;
    }
  }
  for (ParticlesObjectVector::iterator po=pov->begin(); po != pov->end(); po++) {
    if (po->isVisible()) {
      for (int i=0; i < po->npart; i+=po->step) {
        int index=po->index_tab[i];
        if (fabs(p_data->pos[index*3  ]) > coo_max[0]) {
          coo_max[0] = fabs(p_data->pos[index*3  ]);
          i_max[0]   = index;
        }
        if (fabs(p_data->pos[index*3+1]) > coo_max[1]) {
          coo_max[1] = fabs(p_data->pos[index*3+1]);
          i_max[1]   = index;
        }
        if (fabs(p_data->pos[index*3+2]) > coo_max[2]) {
          coo_max[2] = fabs(p_data->pos[index*3+2]);
          i_max[2]   = index;
        }
      }
    }
  }
  if (visible) {
    PRINT_D std::cerr << "Max coordinates \n";
    PRINT_D std::cerr << coo_max[0] << " " << coo_max[1] << " " << coo_max[2] << "\n";
    PRINT_D std::cerr << i_max[0] << " " << i_max[1] << " " << i_max[2] << "\n";
    float max=MAX(MAX(coo_max[0],coo_max[1]),coo_max[2]);
    expandTree(max);
    PRINT_D std::cerr << "Size max = " << size_max << "\n";
  }
}
// ============================================================================
// GLOctree::expandTree
// compute the max square size 
void GLOctree::expandTree(float new_size)
{
  while (size_max < new_size) {
    size_max = size_max << 1;
  }
}
// ============================================================================

}
