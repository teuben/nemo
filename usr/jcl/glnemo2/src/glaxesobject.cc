// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#include "glaxesobject.h"

namespace glnemo {

using namespace std;
// ============================================================================
// Constructor 
GLAxesObject::GLAxesObject()
{
  dplist_index = glGenLists( 1 );
  quadric = gluNewQuadric();
  buildDisplayList();
}
// ============================================================================
// Destructor                                                                  
// Delete display list 
GLAxesObject::~GLAxesObject()
{
  gluDeleteQuadric(quadric);
  glDeleteLists( dplist_index, 1 );
}
// ============================================================================
// display
void GLAxesObject::display(const double * mScreen,const double * mScene, const int width, const int height, 
                           const int loc, const float psize, const bool perspective)
{
  int size=psize*width;
  
  int pwidth,pheight;
  switch (loc) {
  case 0: // bottom right
    pwidth = width-size;
    pheight= 0;
    break;
  case 1: // center
    pwidth = width/2-size/2;
    pheight= height/2-size/2;
    break;
  }

  glPushMatrix ();
  
  // set projection  
  //setProjection( width-size, 0, size, size);
  //setProjection( width/2-size/2, width/2, size, size);
  setProjection( pwidth, pheight, size, size,perspective);
  
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity (); // reset OGL rotations
  glTranslatef (0, 0 , -3);
  
  // apply screen rotation on the whole system
  glMultMatrixd (mScreen);  
  glMultMatrixd (mScene);  
  
  
  //glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  //glDepthMask(GL_FALSE);                               // Lock the Depth Mask so we cant edit it
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);                   // Set the type of blending we want
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  
  GLObject::display(); // call display list
  
  //glDepthMask(GL_TRUE);                                // Unlock the Depth Mask so we can edit it again
  glDisable(GL_BLEND);
  
  glPopMatrix ();
  
}

// ============================================================================
// buildDisplayList()                                            
// Build Display List                                                          
void GLAxesObject::buildDisplayList2()
{
  float ORG[3] = {0,0,0};
  
  float XP[3] = {1,0,0},  YP[3] = {0,1,0},
  ZP[3] = {0,0,1};
  
  // display list
  glNewList( dplist_index, GL_COMPILE );
  
 
  glLineWidth (1.2);
  
  glBegin (GL_LINES);
  glColor3f (1,0,0); // X axis is red.
  glVertex3fv (ORG);
  glVertex3fv (XP );
  glColor3f (0,1,0); // Y axis is green.
  glVertex3fv (ORG);
  glVertex3fv (YP );
  glColor3f (0,0,1); // z axis is blue.
  glVertex3fv (ORG);
  glVertex3fv (ZP );
  glEnd();
  
#if 0  
  glLineWidth (0.2);
  glColor3f (0.5,0.5,0.5); // 
  //quadratic = gluNewQuadric();
  gluQuadricNormals(quadratic, GLU_SMOOTH); 
//  gluQuadricTexture(quadratic, GL_TRUE);
  gluQuadricDrawStyle(quadratic,GLU_LINE);
 
  gluSphere(quadratic,1.f,12,12); 
  //gluDeleteQuadric(quadratic);
#endif
  glEndList();
}

// ============================================================================
// buildDisplayList()                                            
// Build Display List                                                          
void GLAxesObject::buildDisplayList()
{
  float length=1.0;
  float radius=length*0.05;
  //GLfloat color[4];
  // display list
  glNewList( dplist_index, GL_COMPILE );
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  gluQuadricNormals(quadric, GLU_SMOOTH); 
  gluQuadricDrawStyle(quadric, GLU_FILL); //this makes it solid
  // x axis  
  glPushMatrix();
  glColor3f (1,0,0); // x axis is red.
  glRotatef(90.0, 0.0, 1.0, 0.0);
  float color[4];
  color[0] = 0.7f;  color[1] = 0.7f;  color[2] = 1.0f;  color[3] = 1.0f;
  color[0] = 1.f;  color[1] = 0.f;  color[2] = 0.0f;  color[3] = 1.0f;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  buildArrow(length,radius,12);
  glPopMatrix();

  // y axis
  glPushMatrix();
  glColor3f (0,1,0); // y axis is red.
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  color[0] = 1.0f;  color[1] = 0.7f;  color[2] = 0.7f;  color[3] = 1.0f;
  color[0] = 0.0f;  color[1] = 1.f;  color[2] = 0.f;  color[3] = 1.0f;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  buildArrow(length,radius,12);
  glPopMatrix();

  // z axis
  glColor3f (0,0,1); // z axis is blue
  color[0] = 0.7f;  color[1] = 1.0f;  color[2] = 0.7f;  color[3] = 1.0f;
  color[0] = 0.0f;  color[1] = 0.0f;  color[2] = 1.f;  color[3] = 1.0f;
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
  buildArrow(length,radius,12);
  
#if 0  
  glLineWidth (0.2);
  glColor3f (0.5,0.5,0.5); // 
  //quadratic = gluNewQuadric();
  gluQuadricNormals(quadric, GLU_SMOOTH); 
//  gluQuadricTexture(quadratic, GL_TRUE);
  gluQuadricDrawStyle(quadric,GLU_LINE);
 
  gluSphere(quadric,length,12,12); 
  //gluDeleteQuadric(quadratic);
#endif
  glDisable(GL_LIGHTING);
  glEndList();
}

// ============================================================================
// buildArrow()
// Build axes arrow
void GLAxesObject::buildArrow(const float length, const float radius, const int nbSubdivisions)
{
  const float head =  2.5*(radius / length) + 0.1;
  const float coneRadiusCoef = 4.0 - 5.0 * head;

  gluCylinder(quadric, radius, radius, length * (1.0 - head/coneRadiusCoef), nbSubdivisions, 1);
  glTranslatef(0.0, 0.0, length * (1.0 - head));
  gluCylinder(quadric, coneRadiusCoef * radius, 0.0, head * length, nbSubdivisions, 1);
  glTranslatef(0.0, 0.0, -length * (1.0 - head));
}
}
