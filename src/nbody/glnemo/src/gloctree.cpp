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
//
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
#include "gloctree.h"
#include <math.h>
#define LOCAL_DEBUG 0
#include "print_debug.h"
#include <assert.h>
#define MAX(A,B) ((A)>(B)?(A):(B))
// ============================================================================
// Constructor                                                                 
Node::Node(int index_, int obj_)
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
  psv=NULL;
  dplist_index = glGenLists( 1 );    // get a new display list index
  setColor(green);
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
		      ParticlesSelectVector  * psv_)
{
  first=true;
  // copy pointers
  p_data = _p_data;
  psv   = psv_;
  new_data = true;

  // build Display list
  buildDisplayList();
}
// ============================================================================
// GLOctree::update()                                                          
// update tree according store_options conditions                              
void GLOctree::update()
{
  //if (!first) return;
  setActivate(store_options->octree_display);
  // restore default object particles index
  for (int obj=0; obj< (int ) psv->size(); obj++ ) { 
    if ((*psv)[obj].vps->is_visible) {
      (*psv)[obj].vps->defaultIndexTab();
    }
  }
  if ( psv && ! store_options->octree_enable) {
    setActivate(false); // do not display tree
  } else {
    buildDisplayList();
  }
}
// ============================================================================
// GLOctree::build()                                                           
// build octree                                                                
int GLOctree::build()
{
  if (psv && store_options->octree_enable) {
    new_data=false;
    init();
    computeSizeMax();
    p_data->tree_size_max = size_max;
    
    PRINT_D std::cerr << "Start build\n";
    
    for (int obj=0; obj< (int ) psv->size(); obj++ ) { 
      VirtualParticlesSelect * vps = (*psv)[obj].vps;
      if (vps->is_visible) {
        //for (int i=0; i < vps->npart; i+=vps->step_part) {
        for (int i=0; i < vps->ni_index; i++) {
          int index=vps->index_tab[i];
	  assert(index<*p_data->nbody);
          //std::cerr << "Nbody = " << nbody << " Obj = " << obj << " I = " << i << " Index = " << index << "\n";
          insertParticle(index,&root, rmid, 1, obj );
        }
      } 
    }
    
    PRINT_D std::cerr << "stop build Level max= " << level_max << "\n";
    return 1;
  }   
  return 0;
}
// ============================================================================
// GLOctree::insertParticle()                                                  
// insert a new particle in the tree                                           
void GLOctree::insertParticle(int index, Node ** node, float * rmid, int level, int obj)
{
  int bit;
  int octant=0;
  float new_rmid[3];
  //assert(level < 30);
  if (level > 30) {
    PRINT_D ;
    //std::cerr << "Skip particle with index [" << index << "]\n";
    return;
  }
  level_max = MAX(level,level_max);
  if ( ! (*node) ) {
    (*node) = new Node();
  }

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
    p_data->tree_depth[index] = level; // Save particle's level
  }
  else { // octant is not free
    if ((*node)->node[octant]->type ==0) { // it's a node
      insertParticle(index,&((*node)->node[octant]),new_rmid,level+1,obj);
    }
    else { // it's a leaf
      // we must insert particles's leaf itself in a new octant
      int save_index=(*node)->node[octant]->index;
      //std::cerr << "save index="<< save_index <<"\n";
      (*node)->node[octant]->type=0; // leaf become node
      int old_obj = (*node)->node[octant]->obj; // must use old obj !!
      insertParticle(save_index,&((*node)->node[octant]),new_rmid,level+1,old_obj); 
      // we must now insert our particle
      insertParticle(index,&((*node)->node[octant]),new_rmid,level+1,obj);
    }
  }
}
//int aa[10];
// ============================================================================
// GLOctree::buildDisplayList()
void GLOctree::buildDisplayList()
{
  build();   // build octree
  
  // display list
  if (psv && store_options->octree_enable) {
    for (int obj=0; obj< (int ) psv->size(); obj++ ) { 
      VirtualParticlesSelect * vps = (*psv)[obj].vps;
      if (vps->is_visible) {
        vps->resetIndexTab();
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
      (*psv)[tree->obj].vps->addIndexTab(tree->index);
#if 0
      aa[tree->obj]++;
      std::cerr << nbody_keeped << " : level =" << level << " obj= "<< tree->obj << " ----- " << tree->index <<"\n";
#endif
    }
    return 1;
  }
  else { // a node
    
    
    // Build
    float square_minus1=(float ) size_max/(1<<level-1);
    float square=(float ) size_max/(1<<level);
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
  for (int obj=0; obj< (int ) psv->size(); obj++ ) {
    VirtualParticlesSelect * vps = (*psv)[obj].vps;
    if (vps->is_visible) {
      visible = true;
      coo_max[0]= fabs(p_data->pos[vps->index_tab[0]]);
      i_max[0]  = 0;
      coo_max[1]= fabs(p_data->pos[vps->index_tab[0]]);
      i_max[1]  = 0;
      coo_max[2]= fabs(p_data->pos[vps->index_tab[0]]);
      i_max[2]  = 0;
    }
  }
  for (int obj=0; obj< (int ) psv->size(); obj++ ) { 
    VirtualParticlesSelect * vps = (*psv)[obj].vps;
    if (vps->is_visible) {
      //for (int i=0; i < vps->npart; i+=vps->step_part) {
      //  int index=vps->getIndex(i);
      for (int i=0; i < vps->ni_index; i++) {
        int index=vps->index_tab[i];
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
// GLParticlesObject::displayPolygons()                                        
// displayPolygons to create gas like particles effect.                        
// use billboarding technique to display polygons :                            
// - transforms coordinates points according model view matrix                 
// - draw quad (2 triangles) around new coordinates and facing camera          
void GLOctree::displayPolygons(const double * _mModel,const double * mProj, const GLuint _texture,const float _u_max,const float _v_max)
{
  // get variables
  mModel   = _mModel;
  texture  = _texture;
  u_max    = _u_max;
  v_max    = _v_max;

  textureX2 = texture*2;
  // parse tree if activated
  if (psv && store_options->octree_enable) {
    for (int obj=0; obj< (int ) psv->size(); obj++ ) {
      VirtualParticlesSelect * vps = (*psv)[obj].vps;
      if (vps->is_visible) {
        vps->resetIndexTab();
      }
    }
    nbody_keeped=0;
    glEnable( GL_TEXTURE_2D );
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, texture); // Select texture
    // Get Frustum
    frustum.getFC(mModel,mProj);

    hackTreePolygons(root,rmid,1,true);

    glDisable( GL_TEXTURE_2D );
  }
}
// ============================================================================
// GLOctree::hackTreePolygons
int GLOctree::hackTreePolygons(Node * tree, float * rmid, int level, bool check_tree)
{
#define MM(row,col)  mModel[col*4+row]
  
float vv[8][3]; // 8 vertex to store
  //std::cerr << "hi\n";
  float new_rmid[8][3];
  if ( !tree  ) { // free node
    return 0;
  }
  else if (tree->type == 1) { // a leaf
    if (level >= store_options->octree_level) { // selected particles
      (*psv)[tree->obj].vps->addIndexTab(tree->index);
      //computePolygons( tree );
    }
    return 1;
  }
  else { // a node
    
    
    // Build
    float square_minus1=(float ) size_max/(1<<level-1);
    float square=(float ) size_max/(1<<level);
    int fac[2] = { -1,1 };
    for (int i=0;i<2;i++) {
      for (int j=0;j<2;j++) {
        for (int k=0;k<2;k++) {
          int octant=(i<<2)+(j<<1)+k;
          
          // new middle coordinates
          new_rmid[octant][0]=fac[i]*square+rmid[0];
          new_rmid[octant][1]=fac[j]*square+rmid[1];
          new_rmid[octant][2]=fac[k]*square+rmid[2];
          if (check_tree) {
	    // vertex's octant coordinates
            vv[octant][0]=fac[i]*square_minus1+rmid[0];//+fac[i]*textureX2;
            vv[octant][1]=fac[j]*square_minus1+rmid[1];//+fac[j]*textureX2;
            vv[octant][2]=fac[k]*square_minus1+rmid[2];//+fac[k]*textureX2;
            //
          }
          assert(octant<8);
          //hackTreePolygons(tree->node[octant],new_rmid,level+1,true);
        } //k
      } //j
    } //i
    bool inside=true;  // true if octant inside FC
    // must check_tree if not all the vertex of the octant
    // are not in the frustum                             
    if (check_tree) {
      
      int np_inside=0;
      // check out octant in FC?
      for (int i=0; i<8; i++) {
        float x=vv[i][0]+store_options->xtrans;
        float y=vv[i][1]+store_options->ytrans;
        float z=vv[i][2]+store_options->ztrans;
#if 0
        // compute point coordinates according to model via matrix
        float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
        float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
        float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
#endif
        if (frustum.isPointInside(x,y,z)) { //!!mx,my,mz)) {
          np_inside++;
        }
      }
      // check if center of octant is inside
      if (np_inside == 0) {
        float x=rmid[0]+store_options->xtrans;
        float y=rmid[1]+store_options->ytrans;
        float z=rmid[2]+store_options->ztrans;
#if 0
        // compute point coordinates according to model via matrix
        float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
        float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
        float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
#endif
        if (frustum.isPointInside(x,y,z)) { //!!mx,my,mz)) {
          np_inside++;
        }
      }
      if (np_inside==0) {   // the whole octant 
        inside = false;     // is out of the FC 
      }
      else
        if (np_inside==8) {   // the whole octant is in the FC    
          check_tree = false; // no need to keep compute the FC   
        }
    }
    if (inside) {                // octant inside FC     
      for (int i=0; i<8; i++) {  // check out sub octants
        hackTreePolygons(tree->node[i],new_rmid[i],level+1,check_tree);
      }
    }
  }
}
// ============================================================================
// GLOctree::computePolygons( Node * )
int GLOctree::computePolygons(Node * tree)
{

  //glEnable(GL_DEPTH_TEST);	
  setColor((*psv)[tree->obj].vps->col);    // set the color
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),store_options->texture_alpha_color);

  float uv[4][2] = { {0.0,   1.0-v_max}, {0.0,   1.0},{u_max, 1.0},
    {u_max, 1.0-v_max}
  };
  static int modulo=0;
  float rot=0.0;
  // Get Frustum
  //frustum.getFC();
  int visible=0;

  int index=tree->index;
  float
      x=p_data->pos[index*3  ]+store_options->xtrans,
      y=p_data->pos[index*3+1]+store_options->ytrans,
      z=p_data->pos[index*3+2]+store_options->ztrans;
    //w=1.0;

  // compute point coordinates according to model via matrix
  float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
  float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
  float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
  modulo += (modulo%4);
  rot++;
  if (1) { // frustum.isPointInside(mx,my,mz)) {
    float new_texture_size=store_options->texture_size;
    visible++;
    glPushMatrix();
    glTranslatef(mx,my,mz);         // move to the transform particles
    glRotatef(rot,0.0,0.0,1.0);   // rotate triangles around z axis
    glBegin(GL_TRIANGLES);

    // 1st triangle
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    modulo = (modulo+1)%4;
    glVertex3f(-new_texture_size , new_texture_size  ,0. );
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    modulo = (modulo+1)%4;
    glVertex3f(-new_texture_size , -new_texture_size  ,0. );
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    glVertex3f(new_texture_size , -new_texture_size  ,0. );
    // second triangle
    modulo = (modulo+2)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(-new_texture_size , new_texture_size  ,0. );
    modulo = (modulo+2)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(+new_texture_size , -new_texture_size  ,0. );
    modulo = (modulo+1)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(new_texture_size , new_texture_size  ,0. );
    glEnd();
    glPopMatrix();
  }
  return 1;

  //std::cerr << "visible = " << visible << "\n";
}
// ============================================================================
