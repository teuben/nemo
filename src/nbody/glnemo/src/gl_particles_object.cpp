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
GLParticlesObject::GLParticlesObject(const ParticlesData * _p_data,
				     VirtualParticlesSelect * _vps,
                                     float _vel_resize_factor
				    ):GLObject()
{

  // get parameters
  p_data = _p_data;
  vps    = _vps;
  vel_resize_factor = _vel_resize_factor;
  
  texture_size = 0.52;               // default texture size
  texture_alpha_color=255;
  particles_alpha = 255;
  vel_resize_factor = 1;
  vel_dp_list  = glGenLists( 1 );    // get a new display list index
                                     // for the velocity vectors    
  dplist_index = glGenLists( 1 );    // get a new display list index
  PRINT_D perror("on display list");
  PRINT_D cerr << "gl get error= [" << glGetError() << "]\n";
  PRINT_D cerr << "GLParticlesObject My dplist_index = " << dplist_index << "\n";
  buildDisplayList(p_data, vps); // build display list
  buildVelDisplayList(p_data, vps);  // build vel vector display list
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
int GLParticlesObject::updateObject(const ParticlesData * _p_data,
				    VirtualParticlesSelect * _vps,
                                    float _vel_resize_factor
				    )
{
  // get parameters
  p_data            = _p_data;
  vps               = _vps;
  vel_resize_factor = _vel_resize_factor;
  
  buildDisplayList(p_data, vps);     // build particles display list 
  buildVelDisplayList(p_data, vps);  // build vel vector display list
  setColor(vps->col);                // set the color
  is_activated=vps->is_visible;      // Object is visible?
  computeCooMax();                   // compute extrem coordinates
  return 1;
}


// ============================================================================
// GLParticlesObject::displayVelVector()                                       
// display velocity vector                                                     
void GLParticlesObject::displayVelVector()
{
  if (vps->is_visible && p_data->vel) {
    
    display(vel_dp_list);
  }
}
// ============================================================================
// GLParticlesObject::buildVelDisplayList()                                    
// build velocity vector display list                                          
void GLParticlesObject::buildVelDisplayList(const ParticlesData * p_data, 
					 VirtualParticlesSelect * _vps)
{
  if (p_data->vel && _vps->is_visible) { // there are velocity vector to display
    // display list
    glNewList( vel_dp_list, GL_COMPILE );
    //glPushMatrix();
    //glDisable(GL_BLEND);
    glBegin(GL_LINES);
    
    for (int i=0; i < _vps->ni_index; i++) {
      int index=_vps->index_tab[i];    
      float 
      x=p_data->pos[index*3  ],
      y=p_data->pos[index*3+1],
      z=p_data->pos[index*3+2];
      glVertex3f(x , y  ,z );     // draw starting point
      
      float 
      x1=p_data->vel[index*3  ] * vel_resize_factor,
      y1=p_data->vel[index*3+1] * vel_resize_factor,
      z1=p_data->vel[index*3+2] * vel_resize_factor;
      glVertex3f(x+x1 , y+y1  ,z+z1 );  // draw ending point  
    }
    glEnd();
    //glPopMatrix();
    glEndList();
    
  }
}
// ============================================================================
// GLParticlesObject::rebuildDisplayList()                                     
// build particles object display list                                         
void GLParticlesObject::rebuildDisplayList(float _vel_resize_factor)
{
  buildDisplayList(p_data, vps);
  rebuildVelDisplayList(_vel_resize_factor);
}
// ============================================================================
// GLParticlesObject::rebuildVelDisplayList()                                  
// build velocity vectot object display list                                   
void GLParticlesObject::rebuildVelDisplayList(float _vel_resize_factor)
{
  if ( _vel_resize_factor > 0) {
    vel_resize_factor = _vel_resize_factor;
  }
  buildVelDisplayList(p_data, vps);
}
// ============================================================================
// GLParticlesObject::buildDisplayList()                                       
// build particles object display list                                         
void GLParticlesObject::buildDisplayList(const ParticlesData * p_data, 
					 VirtualParticlesSelect * _vps)
{
  // display list
  glNewList( dplist_index, GL_COMPILE );
  glBegin(GL_POINTS);
  
  // draw all the selected points 
  //for (int i=0; i< _vps->npart; i+=vps->step_part) {
  //  int index=_vps->getIndex(i);
    for (int i=0; i < _vps->ni_index; i++) {
      int index=_vps->index_tab[i];    
      float 
      x=p_data->pos[index*3  ],
      y=p_data->pos[index*3+1],
      z=p_data->pos[index*3+2];
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
// displayPolygons to create gas like particles effect.                        
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
  
  glBindTexture(GL_TEXTURE_2D, texture); // Select texture
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
      x=p_data->pos[index*3  ],
      y=p_data->pos[index*3+1],
      z=p_data->pos[index*3+2];
      //w=1.0;
      
    // compute point coordinates according to model via matrix
    float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
    float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    
    //modulo = i%4;
    modulo=0;
    //rot++;
    //rot=0.; // 20% without rotations
    rot=index;
    if (1) { // frustum.isPointInside(mx,my,mz)) {
      float new_texture_size=(float) (p_data->tree_size_max) / (2<<(p_data->tree_depth[index]));
      if (texture_size<new_texture_size) {
        new_texture_size = texture_size;
      }
      else {
        new_texture_size *= p_data->tree_depth[index];
      }
      new_texture_size=texture_size;
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
  }
   
  glDisable( GL_TEXTURE_2D );
  //std::cerr << "visible = " << visible << "\n";
}
// ============================================================================
// GLParticlesObject::computeCooMax()                                          
//  compute extremum coordinates                                               
void GLParticlesObject::computeCooMax()
{
  coo_max[0]= fabs(p_data->pos[vps->index_tab[0]]);
  i_max[0]  = 0;
  coo_max[1]= fabs(p_data->pos[vps->index_tab[0]]);
  i_max[1]  = 0;
  coo_max[2]= fabs(p_data->pos[vps->index_tab[0]]);
  i_max[2]  = 0;
  
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
  PRINT_D cerr << "Max coordinates \n";
  PRINT_D cerr << coo_max[0] << " " << coo_max[1] << " " << coo_max[2] << "\n";
  PRINT_D cerr << i_max[0] << " " << i_max[1] << " " << i_max[2] << "\n";
}
// ============================================================================

