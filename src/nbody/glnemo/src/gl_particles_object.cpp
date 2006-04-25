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
// GLParticlesObject class implementation                                      
//                                                                             
// Manage OpenGL Particles Object                                              
// ============================================================================

#include <assert.h>
#include <iostream>
#include "gl_particles_object.h"
#include <errno.h>
#include <qimage.h>
#include <qstring.h>
#include <math.h>
#if GL_EXT_ENABLE
#include <GL/glext.h>
#define GLX_GLXEXT_LEGACY
#include <GL/glx.h>
#endif
#define LOCAL_DEBUG 0
#include "print_debug.h"

#if GL_EXT_ENABLE
PFNGLPOINTPARAMETERFARBPROC  glPointParameterfARB  = NULL;
PFNGLPOINTPARAMETERFVARBPROC glPointParameterfvARB = NULL;
#endif
using namespace std;
// ============================================================================
// Constructor                                                                 
GLParticlesObject::GLParticlesObject(const int * _nbody, const float * _pos,
				     VirtualParticlesSelect * _vps
				    ):GLObject()
{

  // get parameters
  pos   = _pos;
  vps   = _vps;
  nbody = _nbody;
  texture_size = 0.52;               // default texture size
  texture_alpha_color=255;
  particles_alpha = 255;
  dplist_index = glGenLists( 1 );    // get a new display list index
  PRINT_D perror("on display list");
  PRINT_D cerr << "gl get error= [" << glGetError() << "]\n";
  PRINT_D cerr << "GLParticlesObject My dplist_index = " << dplist_index << "\n";
  buildDisplayList(nbody, pos, vps); // build display list
  setColor(vps->col);                // set the color
  is_activated=vps->is_visible;      // Object is visible?
  computeCooMax();                   // compute extrem coordinates
  
#if GL_EXT_ENABLE 
  // sprites stuff
  glPointParameterfARB  = (PFNGLPOINTPARAMETERFEXTPROC) glXGetProcAddressARB((const GLubyte *) "glPointParameterfARB");
  glPointParameterfvARB = (PFNGLPOINTPARAMETERFVARBPROC) glXGetProcAddressARB((const GLubyte *) "glPointParameterfvARB");
  if( !glPointParameterfARB || !glPointParameterfvARB ) {
    std::cerr << "Error on if( !glPointParameterfARB || !glPointParameterfvARB )\n";
    std::exit(1);
  } 
#endif  

}
// ============================================================================
// Destructor                                                                  
GLParticlesObject::~GLParticlesObject()
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // Really weird problem :
  // if glDeleteLists is called, it's not possible to generate a new
  // list with the command glGenLists, which return all the time 0
  // so I can not delete the object....
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  glDeleteLists( dplist_index, 1 );
  perror("glDeleteLists ~GLParticlesObject()");
}
// ============================================================================
// GLParticlesObject::updateObject()                                           
// update the current object with a new particles range                        
int GLParticlesObject::updateObject(const int * _nbody, const float * _pos,
				    VirtualParticlesSelect * _vps
				    )
{
  // get parameters
  pos   = _pos;
  vps   = _vps;
  nbody = _nbody;
  
  buildDisplayList(nbody, pos, vps); // build display list
  setColor(vps->col);                // set the color
  is_activated=vps->is_visible;      // Object is visible?
  computeCooMax();                   // compute extrem coordinates
  return 1;
}

// ============================================================================
// GLParticlesObject::buildDisplayList()                                       
// build particles object display list                                         
void GLParticlesObject::rebuildDisplayList()
{
  buildDisplayList(nbody, pos, vps);
}
// ============================================================================
// GLParticlesObject::buildDisplayList()                                       
// build particles object display list                                         
void GLParticlesObject::buildDisplayList(const int            * nbody, 
					 const float          * pos, 
					 VirtualParticlesSelect * _vps)
{
  if (nbody);   // remove compiler warning 
  // display list
  glNewList( dplist_index, GL_COMPILE );
  glBegin(GL_POINTS);
  
  // draw all the selected points 
  //for (int i=0; i< _vps->npart; i+=vps->step_part) {
  //  int index=_vps->getIndex(i);
    for (int i=0; i < _vps->ni_index; i++) {
      int index=_vps->index_tab[i];    
      float 
      x=pos[index*3  ],
      y=pos[index*3+1],
      z=pos[index*3+2];
      //std::cerr << x << " " << y << " " << z << " " << index << "\n";
    // One point
    glVertex3f(x , y  ,z );
  }
  glEnd();
  glEndList();
}
#if GL_EXT_ENABLE
// ============================================================================
// GLParticlesObject::displaySprites()                                         
// displaySprites instead of polygones. Sprites does not need billboarding and 
// should give a faster rendering but it seems not to be the case.... Morover  
// sprites have a limited size which give a small cloud of gaz.                
void GLParticlesObject::displaySprites(GLuint texture)
{
  if (!is_activated)
    return;
  static bool first=true;
  if (first) {
    first=false;
    //loadImage();
  }
  glEnable( GL_TEXTURE_2D );
  //glDisable( GL_TEXTURE_2D );
  glDisable(GL_DEPTH_TEST);
  //glEnable(GL_DEPTH_TEST);
  setColor(vps->col);                // set the color
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),texture_alpha_color);
  
  glBindTexture(GL_TEXTURE_2D, texture);				// Select Our Texture  
  // This is how will our point sprite's size will be modified by 
  // distance from the viewer
  float quadratic[] =  { 1.0f, 0.0f, 0.01f };
  glPointParameterfvARB( GL_POINT_DISTANCE_ATTENUATION_ARB, quadratic );

  // Query for the max point size supported by the hardware
  float maxSize = 0.0f;
  glGetFloatv( GL_POINT_SIZE_MAX_ARB, &maxSize );

  // Clamp size to 100.0f or the sprites could get a little too big on some  
  // of the newer graphic cards. My ATI card at home supports a max point 
  // size of 1024.0f!
  if( maxSize > 100.0f )
      maxSize = 100.0f;
  glPointSize( maxSize );

  // The alpha of a point is calculated to allow the fading of points 
  // instead of shrinking them past a defined threshold size. The threshold 
  // is defined by GL_POINT_FADE_THRESHOLD_SIZE_ARB and is not clamped to 
  // the minimum and maximum point sizes.
  glPointParameterfARB( GL_POINT_FADE_THRESHOLD_SIZE_ARB, 60.0f );
  glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, 1.0f );
  glPointParameterfARB( GL_POINT_SIZE_MAX_ARB, maxSize );

  // Specify point sprite texture coordinate replacement mode for each 
  // texture unit
  glTexEnvf( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );

  //
  // Render point sprites...
  //
  glEnable( GL_POINT_SPRITE_ARB );
  display();   
  glDisable( GL_POINT_SPRITE_ARB );
  glDisable( GL_TEXTURE_2D );

}
#endif
// ============================================================================
// GLParticlesObject::displayPolygons()                                        
// displayPolygons to create gaz like particles effect.                        
// use billboarding technique to display polygons :                            
// - transforms coordinates points according model view matrix                 
// - draw quad (2 triangles) around new coordinates and facing camera          
void GLParticlesObject::displayPolygons(const double * mModel,GLuint texture,float u_max, float v_max)
{

#define MM(row,col)  mModel[col*4+row]  
  if (!is_activated)
    return;
  static bool first=true;
  //float sizeq=.52;
  if (first) {
    first=false;
    //loadImage();
  }
#if 0  
  // specify that back facing polygons are to be culled
   glCullFace( GL_BACK );
   // enable culling
   glEnable( GL_CULL_FACE );  
   glEnable( GL_POLYGON_SMOOTH );
#endif    

  glEnable( GL_TEXTURE_2D );
  glDisable(GL_DEPTH_TEST);
  //glEnable(GL_DEPTH_TEST);	
  setColor(vps->col);    // set the color
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),texture_alpha_color);
  
  glBindTexture(GL_TEXTURE_2D, texture); // Select Our Texture
  float uv[4][2] = { {0.0,   1.0-v_max}, {0.0,   1.0},{u_max, 1.0},
                     {u_max, 1.0-v_max}
                   };
  int modulo;
  float rot=0.0;
  // Get Frustum
  //frustum.getFC();
  int visible=0;
  // draw all the selected points           
  // method to shuffle triangles orientation
  // rotate uv map and rotate triangles     
  //for (int i=0; i < vps->npart; i+=vps->step_part) {
  //  int index=vps->getIndex(i);
  for (int i=0; i < vps->ni_index; i++) {
      int index=vps->index_tab[i]; 
      float 
      x=pos[index*3  ],
      y=pos[index*3+1],
      z=pos[index*3+2];
      //w=1.0;
      
    // compute point coordinates according to model via matrix  
    float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
    float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    modulo = i%4;
    rot++;
    if (1) { // frustum.isPointInside(mx,my,mz)) {
	visible++;
	glPushMatrix();
	glTranslatef(mx,my,mz);         // move to the transform particles
	glRotatef(rot,0.0,0.0,1.0);   // rotate triangles around z axis
	glBegin(GL_TRIANGLES);
	
	// 1st triangle
	glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
	modulo = (modulo+1)%4;
	glVertex3f(-texture_size , texture_size  ,0. );
	glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
	modulo = (modulo+1)%4;
	glVertex3f(-texture_size , -texture_size  ,0. );
	glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
	glVertex3f(texture_size , -texture_size  ,0. );
	// second triangle
	modulo = (modulo+2)%4;
	glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
	glVertex3f(-texture_size , texture_size  ,0. );
	modulo = (modulo+2)%4;
	glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
	glVertex3f(+texture_size , -texture_size  ,0. );
	modulo = (modulo+1)%4;
	glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
	glVertex3f(texture_size , texture_size  ,0. );     
	glEnd();
	glPopMatrix();    
    }
  }
   
  glDisable( GL_TEXTURE_2D );
  //std::cerr << "visible = " << visible << "\n";
}
// ============================================================================
// GLParticlesObject::computeCooMax()                                          
//  compute extremum coordinates                                               
void GLParticlesObject::computeCooMax()
{
  coo_max[0]= fabs(pos[vps->index_tab[0]]);
  i_max[0]  = 0;
  coo_max[1]= fabs(pos[vps->index_tab[0]]);
  i_max[1]  = 0;
  coo_max[2]= fabs(pos[vps->index_tab[0]]);
  i_max[2]  = 0;
  
  //for (int i=0; i < vps->npart; i+=vps->step_part) {
  //  int index=vps->getIndex(i);
  for (int i=0; i < vps->ni_index; i++) {
    int index=vps->index_tab[i];
    if (fabs(pos[index*3  ]) > coo_max[0]) {
      coo_max[0] = fabs(pos[index*3  ]);
      i_max[0]   = index;
    }
    if (fabs(pos[index*3+1]) > coo_max[1]) {
      coo_max[1] = fabs(pos[index*3+1]);
      i_max[1]   = index;
    }
    if (fabs(pos[index*3+2]) > coo_max[2]) {
      coo_max[2] = fabs(pos[index*3+2]);
      i_max[2]   = index;
    }
  }
  PRINT_D cerr << "Max coordinates \n";
  PRINT_D cerr << coo_max[0] << " " << coo_max[1] << " " << coo_max[2] << "\n";
  PRINT_D cerr << i_max[0] << " " << i_max[1] << " " << i_max[2] << "\n";
}
// ============================================================================

