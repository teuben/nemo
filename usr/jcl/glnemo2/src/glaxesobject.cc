// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
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
void GLAxesObject::display(const double * mScreen,const double * mScene, const int width)
{
  int size=100;
  glPushMatrix ();
  // set projection  
  setProjection( width-size, 0, size, size);
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity (); // reset OGL rotations
  glTranslatef (0, 0 , -3);
  
  // apply screen rotation on the whole system
  glMultMatrixd (mScreen);  
  glMultMatrixd (mScene);  
  glEnable(GL_BLEND);
  //glEnable( GL_DEPTH_TEST );
  GLObject::display(); // call display list
  
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
  float radius=length*0.03;
  //GLfloat color[4];
  // display list
  glNewList( dplist_index, GL_COMPILE );
  gluQuadricNormals(quadric, GLU_SMOOTH); 
  // x axis  
  glPushMatrix();
  glColor3f (1,0,0); // x axis is red.
  glRotatef(90.0, 0.0, 1.0, 0.0);
  buildArrow(length,radius,12);
  glPopMatrix();

  // y axis
  glPushMatrix();
  glColor3f (0,1,0); // y axis is red.
  glRotatef(-90.0, 1.0, 0.0, 0.0);
  buildArrow(length,radius,12);
  glPopMatrix();

  // z axis
  glColor3f (0,0,1); // z axis is blue
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
