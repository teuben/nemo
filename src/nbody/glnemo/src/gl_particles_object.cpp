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

#define LOCAL_DEBUG 0
#include "print_debug.h"

using namespace std;
// ============================================================================
// Constructor
GLParticlesObject::GLParticlesObject(const int * _nbody, const float * _pos,
				     const ParticlesRange * _prv
				    ):GLObject()
{

  // get parameters
  pos   = _pos;
  prv   = _prv;
  nbody = _nbody;
  
  dplist_index = glGenLists( 1 );    // get a new display list index
  PRINT_D perror("on display list");
  PRINT_D cerr << "gl get error= [" << glGetError() << "]\n";
  PRINT_D cerr << "GLParticlesObject My dplist_index = " << dplist_index << "\n";
  buildDisplayList(nbody, pos, prv); // build display list
  setColor(prv->col);                // set the color
  is_activated=prv->is_visible;      // Object is visible?
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
// update the current object with a new particle range
int GLParticlesObject::updateObject(const int * _nbody, const float * _pos,
				    const ParticlesRange * _prv
				    )
{
  // get parameters
  pos   = _pos;
  prv   = _prv;
  nbody = _nbody;

  buildDisplayList(nbody, pos, prv); // build display list
  setColor(prv->col);                // set the color
  is_activated=prv->is_visible;      // Object is visible?
  return 1;
}
// ============================================================================
// 
void GLParticlesObject::buildDisplayList(const int            * nbody, 
					 const float          * pos, 
					 const ParticlesRange * prv)
{
  
  if (nbody);   // remove compiler warning
  

  // display list
  glNewList( dplist_index, GL_COMPILE );

  glBegin(GL_POINTS);
  
  // draw all the selected points 
  for (int i=prv->first_part; i<=prv->last_part; i+=prv->step_part) {
    float 
      x=pos[i*3  ],
      y=pos[i*3+1],
      z=pos[i*3+2];
    // One point
    glVertex3f(x , y  ,z );
  }
  glEnd();
  glEndList();
  
}
// ============================================================================
// displayPolygons                                                             
// use billboarding technique to display polygons :                            
// - transforms coordinates points according model view matrix                 
// - draw quad (2 triangles) around new coordinates and facing camera          
void GLParticlesObject::displayPolygons(const double * mModel,GLuint texture,float u_max, float v_max)
{

#define MM(row,col)  mModel[col*4+row]  
  if (!is_activated)
    return;
  static bool first=true;
  float sizeq=.52;
  if (first) {
    first=false;
    //loadImage();
  }
  // display list
  //glNewList( dplist_index, GL_COMPILE );
  
  glEnable( GL_TEXTURE_2D );
  //glDisable( GL_TEXTURE_2D );
  glDisable(GL_DEPTH_TEST);
  //glEnable(GL_DEPTH_TEST);
  setColor(prv->col);                // set the color
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),255);
  
  glBindTexture(GL_TEXTURE_2D, texture);				// Select Our Texture
  //glBegin(GL_QUADS);
  glBegin(GL_TRIANGLES);
  // draw all the selected points
  for (int i=prv->first_part; i<=prv->last_part; i+=prv->step_part) {
    float 
      x=pos[i*3  ],
      y=pos[i*3+1],
      z=pos[i*3+2];
      //w=1.0;
      
    // compute point coordinates according to model vie matrix  
    float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
    float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    //sizeq=sizeq*mz;
    //float nsizeq=z*sizeq;
    //glBegin(GL_QUADS);
    
    //glBegin(GL_TRIANGLES);
    
    glTexCoord2f(0.0,   1.0-v_max);
    glVertex3f(mx-sizeq , my+sizeq  ,mz );
    glTexCoord2f(0.0,   1.0);
    glVertex3f(mx-sizeq , my-sizeq  ,mz );
    glTexCoord2f(u_max, 1.0);
    glVertex3f(mx+sizeq , my-sizeq  ,mz );
#if 0   
    glTexCoord2f(u_max, 1.0-v_max);
    glVertex3f(mx+sizeq , my+sizeq  ,mz );
#else
    glTexCoord2f(0.0,   1.0-v_max);
    glVertex3f(mx-sizeq , my+sizeq  ,mz );
    glTexCoord2f(u_max, 1.0);
    glVertex3f(mx+sizeq , my-sizeq  ,mz );
    glTexCoord2f(u_max, 1.0-v_max);
    glVertex3f(mx+sizeq , my+sizeq  ,mz );       
#endif    


  }
  glEnd();
  glDisable( GL_TEXTURE_2D );
  //glEndList(); // end of displaylist
}

//
