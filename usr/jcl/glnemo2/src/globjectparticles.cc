// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2008                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           P�le de l'Etoile, site de Ch�teau-Gombert                         
//           38, rue Fr�d�ric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include <GL/glew.h>
#include "globjectparticles.h"
#include "particlesobject.h"
#include "particlesdata.h"
#include "globaloptions.h"
#include "gltexture.h"
#include "glwindow.h"
#include <assert.h>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QTime>
#include <vector>
#include <sstream>
#include <algorithm>

#define BENCH 1
#define GLDRAWARRAYS 1
namespace glnemo {
float DENS_MIN=-1.;;
float DENS_MAX=0.000006;
float TMIN=0.;
float TMAX=1.;
int index_min, index_max;
// ============================================================================
// constructor                                                                 
GLObjectParticles::GLObjectParticles(GLTextureVector * _gtv ):GLObject()
{
  dplist_index = glGenLists( 1 );    // get a new display list index
  texture = NULL;                    // no texture yet              
  gtv = _gtv;
  if (GLWindow::GLSL_support) {
    glGenBuffersARB(1,&vbo_pos);
//    glGenBuffersARB(1,&vbo_color);
    glGenBuffersARB(1,&vbo_size);
    glGenBuffersARB(1,&vbo_index);
    glGenBuffersARB(1,&vbo_index2);
  }
  indexes_sorted = NULL;
  nind_sorted = 0;
}
// ============================================================================
// constructor                                                                 
GLObjectParticles::GLObjectParticles(const ParticlesData   * _part_data,
                                     ParticlesObject * _po,
                                     const GlobalOptions   * _go,
				     GLTextureVector * _gtv):GLObject()
{
  dplist_index = glGenLists( 1 );    // get a new display list index
  vel_dp_list  = glGenLists( 1 );    // get a new display vel list
  orb_dp_list  = glGenLists( 1 );    // get a new display orb list
  if (GLWindow::GLSL_support) {
    glGenBuffersARB(1,&vbo_pos);     // get Vertex Buffer Object
//    glGenBuffersARB(1,&vbo_color);   // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_size);    // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_index);   // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_index2);
  }
  indexes_sorted = NULL;
  nind_sorted = 0;
  texture = NULL;                    // no texture yet
  gtv = _gtv;
  //setTexture(":/images/textures/smoke10.png");
  assert(gtv->size()>0);
  setTexture(0);
  update(_part_data,_po,_go);
}

// ============================================================================
// desstructor                                                                 
GLObjectParticles::~GLObjectParticles()
{
}
// ============================================================================
// update                                                                      
void GLObjectParticles::display(const double * mModel, int win_height)
{
  if (po->isVisible()) {
    // display particles
    if (po->isPartEnable()) {
      glEnable(GL_BLEND);
      glEnable(GL_POINT_SMOOTH);
      glPointSize((float) po->getPartSize());
      GLObject::updateAlphaSlot(po->getPartAlpha());
      GLObject::setColor(po->getColor());
      if (GLWindow::GLSL_support) displayVboPoints();
      else GLObject::display();
      glDisable(GL_BLEND);
    }
    // display velocities
    if (po->isVelEnable() && part_data->vel) {
      glEnable(GL_BLEND);
      GLObject::updateAlphaSlot(po->getVelAlpha());
      GLObject::setColor(po->getColor());
      GLObject::display(vel_dp_list);
      glDisable(GL_BLEND);
    }
    // display sprites
    if (po->isGazEnable() && texture) {
      if (GLWindow::GLSL_support && po->isGazGlsl()) {
        //glEnable(GL_CULL_FACE);
        //glCullFace(GL_FRONT);
        displayVboSprites(win_height,true);
        //glCullFace(GL_BACK);
        //displayVboSprites(win_height,false);
      }
      else displaySprites(mModel); 
    }
  }
  if (po->isOrbitsEnable()) {
    glEnable (GL_LINE_SMOOTH);
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
    glLineWidth (1.0);
    setColor(Qt::yellow);
    GLObject::display(orb_dp_list);
    glDisable(GL_BLEND);
  }
}
// ============================================================================
// update                                                                      
void GLObjectParticles::update( const ParticlesData   * _part_data,
                                ParticlesObject * _po,
                                const GlobalOptions   * _go)
{
  // update variables
  part_data = _part_data;
  po        = _po;
  go        = _go;
  // color
  mycolor   = po->getColor();
  
  if (DENS_MIN == -1 && part_data->rho) {
    DENS_MIN = part_data->getMinRho();
    DENS_MAX = part_data->getMaxRho();
    TMIN     = part_data->getMinTemp();
    TMAX     = part_data->getMaxTemp();
  }
  //setColor();
  if (!GLWindow::GLSL_support) buildDisplayList();
  buildVelDisplayList();
  if (po->isOrbitsRecording()) {
    po->addOrbits(part_data);
    buildOrbitsDisplayList();
  } else {
    glNewList( orb_dp_list, GL_COMPILE );
    glEndList();
  }
  if (GLWindow::GLSL_support) {

    buildVboPos();
    buildVboColor();
    buildVboSize2();
#if ! GLDRAWARRAYS
    sortByDensity();
#endif
  }
  std::cerr << "DENS min="<< DENS_MIN<<" DENS max="<< DENS_MAX << "\n";
}
// ============================================================================
// updateVbo                                                                   
void GLObjectParticles::updateVbo()
{
  if (GLWindow::GLSL_support && po->isGazEnable()) {

    if (go->render_mode == 1 || go->render_mode == 2) {
      if (go->dens_local) {
        DENS_MAX = po->getMaxDensity();
        DENS_MIN = po->getMinDensity();
        TMIN     = po->getMinTemperature();
        TMAX     = po->getMaxTemperature();
      } else {
        DENS_MAX = go->dens_max_glob;
        DENS_MIN = go->dens_min_glob;
      }
      // We must re ordering indexes
      if (part_data->rho) {
#if ! GLDRAWARRAYS
        rho.clear();
#endif
        for (int i=0; i < po->npart; i+=po->step) {         
          if (part_data->rho) {
#if ! GLDRAWARRAYS
            int index=po->index_tab[i];
            GLObjectIndexTab myrho;
            myrho.index   = index;
            myrho.value   = part_data->rho[index];
            myrho.i_point = i;
            rho.push_back(myrho);
#endif
          }
        }
      }
      
      buildVboColor();
#if ! GLDRAWARRAYS
      sortByDensity();	
#else // GLDRAWARRAYS
      buildVboSize2();
#endif

    }
  }
}
// ============================================================================
// updateVbo                                                                   
void GLObjectParticles::updateColorVbo()
{
  if (1||go->render_mode == 1 || go->render_mode == 2) {
    buildVboColor();
  }
}
// ============================================================================
// update                                                                      
void GLObjectParticles::updateVel()
{
  buildVelDisplayList();
}
// ============================================================================
// buildDisplayList                                                            
void GLObjectParticles::buildDisplayList()
{
  QTime tbench;
  tbench.restart();
  // display list
  glNewList( dplist_index, GL_COMPILE );
  glBegin(GL_POINTS);
  
  // draw all the selected points 
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    float
      x=part_data->pos[index*3  ],
      y=part_data->pos[index*3+1],
      z=part_data->pos[index*3+2];
    // One point
    glVertex3f(x , y  ,z );
  }
  glEnd();
  glEndList();
  if (BENCH) qDebug("Time elapsed to build Pos Display list: %f s", tbench.elapsed()/1000.);
}
// ============================================================================
// buildDisplayList                                                            
void GLObjectParticles::buildVelDisplayList()
{
  if (part_data->vel) {
    QTime tbench;
    tbench.restart();
    // display list
    glNewList( vel_dp_list, GL_COMPILE );
    glBegin(GL_LINES);
    
    // draw all the selected points
    const float vfactor = po->getVelSize() / part_data->getMaxVelNorm();
    for (int i=0; i < po->npart; i+=po->step) {
      int index=po->index_tab[i];
      float
        x=part_data->pos[index*3  ],
        y=part_data->pos[index*3+1],
        z=part_data->pos[index*3+2];
      // Draw starting point
      glVertex3f(x , y  ,z );
      float
        x1=part_data->vel[index*3  ] * vfactor,
        y1=part_data->vel[index*3+1] * vfactor,
        z1=part_data->vel[index*3+2] * vfactor;
        glVertex3f(x+x1 , y+y1  ,z+z1 );  // draw ending point
    }
    glEnd();
    glEndList();
    if (BENCH) qDebug("Time elapsed to build Vel Display list: %f s", tbench.elapsed()/1000.);
  }
}
// ============================================================================
// buildOrbitsDisplayList                                                            
void GLObjectParticles::buildOrbitsDisplayList()
{
  OrbitsVector oo = po->ov;
  glNewList( orb_dp_list, GL_COMPILE );
  
  // loop on orbits_max
  for (OrbitsVector::iterator oiv =oo.begin();oiv!=oo.end() ; oiv++) {
    glBegin(GL_LINE_STRIP);
    // loop on orbit_history
    for (OrbitsList::iterator oil=(*oiv).begin(); oil != (*oiv).end(); oil++){
      glVertex3f(oil->x(),oil->y(),oil->z());
    }
    glEnd();
  }
  
  glEndList();
}
// ============================================================================
// ============================================================================
// compare 2 elements                                                          
int GLObjectParticles::compareZ( const void * a, const void * b )
{
  float (*pa)[3], (*pb)[3];
  pa = ( float(*)[3] ) a; 
  pb = ( float(*)[3] ) b;
  return (int)(*pb[2] - *pa[2]);
}
// ============================================================================
// buildVboPos                                                                 
// Build Vector Buffer Object for positions array                              
void GLObjectParticles::buildVboPos()
{
  QTime tbench;
  tbench.restart();
  nvert_pos=0;
  
  GLfloat* vertices = new GLfloat[((po->npart/po->step)+1)*3]; // create vertex array
  
  rho.clear();    // clear rho density vector
  zdepth.clear(); // clear zdepth vector     
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    if (part_data->rho) {
      GLObjectIndexTab myrho;
      myrho.index   = index;
      myrho.value   = part_data->rho[index];
      myrho.i_point = i;
      rho.push_back(myrho);
    }
#if 0 // used if z depth test activated
    GLObjectIndexTab myz;
    myz.index = index;
    myz.i_point = i;
    zdepth.push_back(myz);
#endif
  }
  // sort by density
#if GLDRAWARRAYS
  sort(rho.begin(),rho.end(),GLObjectIndexTab::compareLow);
#endif
  // fill vertices array sorted by density
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    vertices[nvert_pos*3+0] = part_data->pos[index*3  ];
    vertices[nvert_pos*3+1] = part_data->pos[index*3+1];
    vertices[nvert_pos*3+2] = part_data->pos[index*3+2];
    nvert_pos++;
  }
  //qsort(vertices,nvert_pos,sizeof(float)*3,compareZ);
  // bind VBO in order to use
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  //std::cerr << "vbo_pos = " << vbo_pos << "\n";
  assert( nvert_pos <= (po->npart/po->step)+1);
  // upload data to VBO
#if 0
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 3 * sizeof(float), vertices, GL_STATIC_DRAW_ARB);
#else
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * (3 + 4 )* sizeof(float), 0, GL_STATIC_DRAW_ARB);
  glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, nvert_pos * 3 * sizeof(float), vertices);
#endif
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  
  delete [] vertices;
  if (BENCH) qDebug("Time elapsed to build VBO pos: %f s", tbench.elapsed()/1000.);
}
// ============================================================================
// buildVboColor                                                               
// Build Vector Buffer Object for colors array                                 
void GLObjectParticles::buildVboColor()
{
  buildVboColorGasGasSorted();
  //buildVboColorTempGasSorted();
}
// ============================================================================
// buildVboColor                                                               
// Build Vector Buffer Object for colors array                                 
void GLObjectParticles::buildVboColorTempGasSorted()
{
  QTime tbench,vbobench;
  tbench.restart();
  nvert_pos=0;
  GLfloat * colors = new GLfloat[((po->npart/po->step)+1)*4]; // create colors array

  std::cerr << "buildVboColor Densmin="<<DENS_MIN<<"\n";
  std::cerr << "buildVboColor Densmax="<<DENS_MAX<<"\n";
  // density
  float LOGDMIN = log(DENS_MIN);
  float LOGDMAX = log(DENS_MAX);
  float DIFMM   = LOGDMAX - LOGDMIN;
  float INVDIFMM= 1./DIFMM;
  //float LOGMPMINRHO = log(part_data->getMinRho());
  //float LOGMPMAXRHO = log(part_data->getMaxRho());
  // temperature
  
  // Compute TMIN and TMAX according to the selected density
  bool first=true;
  
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    if (part_data->rho && part_data->temp) {
      float partrhoindex=part_data->rho[index];
      if (partrhoindex >= DENS_MIN && partrhoindex <= DENS_MAX) {
        if (first) {
          TMIN = TMAX = part_data->temp[index];
          first=false;
        }
        TMIN = std::min(TMIN,part_data->temp[index]);
        TMAX = std::max(TMIN,part_data->temp[index]);
      }
    }
  }
  std::cerr << "TMIN = "<< TMIN <<" TMAX = " << TMAX << "\n";
  float LOGTMIN = log(TMIN);
  float LOGTMAX = log(TMAX);
  float TDIFMM   = LOGTMAX - LOGTMIN;
  float INVTDIFMM= 1./TDIFMM;
  float LOGMPMINTEMP = log(part_data->getMinTemp());
  float LOGMPMAXTEMP = log(part_data->getMaxTemp());
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    float log_temp, log_rho;
    log_temp=0.0;
    if (part_data->rho && part_data->temp) {
      float partrhoindex=part_data->rho[index];
      float logpri=log(partrhoindex);
      float parttempindex=part_data->temp[index];
      float logpti=log(parttempindex);
      
#if 0
      log_temp=((log(part_data->rho[index])-log(part_data->getMinRho()))/
          (log(part_data->getMaxRho())-log(part_data->getMinRho())));
#endif
      if (partrhoindex >= DENS_MIN && partrhoindex <= DENS_MAX) {

        // log of the density
        log_rho=((logpri-LOGDMIN)*INVDIFMM);
        // log of temperature
        log_temp=((logpti-LOGTMIN)*INVTDIFMM);
          
        if (go->reverse_cmap) { // reversed colormap
          log_temp=((LOGTMAX-logpti)*INVTDIFMM);
        }
        else {                  // normal colormap

          if (! go->constant_cmap ) {
              // good density colors according to the selection
            log_temp=((logpti- LOGMPMINTEMP)/
                (LOGMPMAXTEMP-LOGMPMINTEMP));
          } else {
              // keep color map constant between 2 frames
            log_temp=((logpti-LOGTMIN)/
                (TDIFMM));
          }

        }
      }
      else {
        log_temp = 0.;
      }
    }
    GLObject::setColor(po->getColor());
    float powcolor=go->powercolor;
    float powalpha=go->poweralpha;
    float col=mycolor.redF();
    //int cindex=(int) (log_temp*(R.size()-1));
    int cindex=(int) (log_temp*(go->R->size()-1));
    bool ok=true;
    if (part_data->rho && part_data->rho[index] == -1) ok=false;
    if (part_data->rho && ok) {
      colors[nvert_pos*4+0] =  pow((*go->R)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+0] = mycolor.redF();
    }
    col=mycolor.greenF();
    if (part_data->rho && ok) {
      colors[nvert_pos*4+1] =  pow((*go->G)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+1] = mycolor.greenF();
    }
    col=mycolor.blueF();
    if (part_data->rho && ok) {
      colors[nvert_pos*4+2] =  pow((*go->B)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+2] = mycolor.blueF();
    }

    if (part_data->rho) {
      if (part_data->rho[index] != -1) { // Normalize rho for alpha color
#if 1
	//colors[nvert_pos*4+3] = pow(alpha * log_rho,2.);
        colors[nvert_pos*4+3] = pow(log_rho,powalpha);//* log_rho;
#else
        // Normalize rho for alpha color
        float alpha=po->getGazAlpha()/255.;
        colors[nvert_pos*4+3] = alpha+(1-alpha)*(part_data->rho[index]-part_data->getMinRho())/
            (part_data->getMaxRho()-part_data->getMinRho());
#endif
        if (colors[nvert_pos*4+3] > 1.0 || colors[nvert_pos*4+3] <0.0) {
          std::cerr << "Erreur alpha = "<<colors[nvert_pos*4+3]<<"\n";
        }
      } 
      else {
        colors[nvert_pos*4+3] = 1.0;
      }
    } 
    else {
      colors[nvert_pos*4+3] = 1.0;
    }
    nvert_pos++;
  }
  vbobench.restart();
  // bind VBO in order to use
  //glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_color);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  //std::cerr << "vbo_pos = " << vbo_pos << "\n";
  assert( nvert_pos <= (po->npart/po->step)+1);
  // upload data to VBO
#if 0
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 4 * sizeof(float), colors, GL_STATIC_DRAW_ARB);
#else
  glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 3 * sizeof(float), nvert_pos * 4 * sizeof(float), colors);
#endif
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

  if (BENCH) qDebug("WALL Time elapsed to build VBO color: %f s", tbench.elapsed()/1000.);
  if (BENCH) qDebug("TRANSFERT Time elapsed to build VBO color: %f sec %f MB/s", vbobench.elapsed()/1000.,
      (nvert_pos * 4 * sizeof(float)/1024./1024.)/(vbobench.elapsed()/1000.));
  
  delete [] colors;

}
// ============================================================================
// buildVboColor                                                               
// Build Vector Buffer Object for colors array                                 
void GLObjectParticles::buildVboColorGasGasSorted()
{
  QTime tbench,vbobench;
  tbench.restart();
  nvert_pos=0;
  GLfloat * colors = new GLfloat[((po->npart/po->step)+1)*4]; // create colors array
#if 0
  QColor mycolor(255,255,255);
  QColor leftcolor= po->getColor();
  int VV=256;
  float R0 = (leftcolor.red()-mycolor.red())/(VV-1.);
  float G0 = (leftcolor.green()-mycolor.green())/(VV-1.);
  float B0 = (leftcolor.blue()-mycolor.blue())/(VV-1.);
  //std::cerr << "RO="<<R0<<"  B0="<<B0<<"  G0="<<G0<<"\n";
  for (int i=0; i<VV; i++) {
    int r = (int)(leftcolor.red()-i*R0);
    int g = (int)(leftcolor.green()-i*G0);
    int b = (int)(leftcolor.blue()-i*B0);
      //std::cerr << "r="<<r<<"  g="<<g<<"  b="<<b<<"\n";   
    R.push_back(r);
    G.push_back(g);
    B.push_back(b);
  }
#endif
  std::cerr << "buildVboColor Densmin="<<DENS_MIN<<"\n";
  std::cerr << "buildVboColor Densmax="<<DENS_MAX<<"\n";
  
  float LOGDMIN = log(DENS_MIN);
  float LOGDMAX = log(DENS_MAX);
#if 0  
  DENS_MIN = part_data->getMinRho();
  DENS_MAX = part_data->getMaxRho();CCCC
  LOGDMIN = log(part_data->getMinRho());
  LOGDMAX = log(part_data->getMaxRho());
#endif
  
  float DIFMM   = LOGDMAX - LOGDMIN;
  float INVDIFMM= 1./DIFMM;
  float LOGMPMINRHO = log(part_data->getMinRho());
  float LOGMPMAXRHO = log(part_data->getMaxRho());
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    float log_rho;
    log_rho=0.0;
    if (part_data->rho) {
      float partrhoindex=part_data->rho[index];
      float logpri=log(partrhoindex);
#if 0
      log_rho=((log(part_data->rho[index])-log(part_data->getMinRho()))/
          (log(part_data->getMaxRho())-log(part_data->getMinRho())));
#endif
  
      if (partrhoindex >= DENS_MIN && partrhoindex <= DENS_MAX) {

          // log of the density
        log_rho=((logpri-LOGDMIN)*INVDIFMM);
          
          if (go->reverse_cmap) { // reversed colormap
            log_rho=((LOGDMAX-logpri)*INVDIFMM);
          }
          else {                  // normal colormap

            if (! go->constant_cmap ) {
              // good density colors according to the selection
              log_rho=((logpri- LOGMPMINRHO)/
                  (LOGMPMAXRHO-LOGMPMINRHO));
            } else {
              // keep color map constant between 2 frames
              log_rho=((logpri-LOGDMIN)/
                  (DIFMM));
            }

          }
	}
        else {
           log_rho = 0.;
        }
    }
    GLObject::setColor(po->getColor());
    float powcolor=go->powercolor;
    float powalpha=go->poweralpha;
    float col=mycolor.redF();
    //int cindex=(int) (log_rho*(R.size()-1));
    int cindex=(int) (log_rho*(go->R->size()-1));
    bool ok=true;
    if (part_data->rho && part_data->rho[index] == -1) ok=false;
    if (part_data->rho && ok) {
      colors[nvert_pos*4+0] =  pow((*go->R)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+0] = mycolor.redF();
    }
    col=mycolor.greenF();
    if (part_data->rho && ok) {
      colors[nvert_pos*4+1] =  pow((*go->G)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+1] = mycolor.greenF();
    }
    col=mycolor.blueF();
    if (part_data->rho && ok) {
      colors[nvert_pos*4+2] =  pow((*go->B)[cindex],powcolor);
    } else {
      colors[nvert_pos*4+2] = mycolor.blueF();
    }

    if (part_data->rho) {
      if (part_data->rho[index] != -1) { // Normalize rho for alpha color
#if 1
	//colors[nvert_pos*4+3] = pow(alpha * log_rho,2.);
        colors[nvert_pos*4+3] = pow(log_rho,powalpha);//* log_rho;
#else
        // Normalize rho for alpha color
        float alpha=po->getGazAlpha()/255.;
	colors[nvert_pos*4+3] = alpha+(1-alpha)*(part_data->rho[index]-part_data->getMinRho())/
	    (part_data->getMaxRho()-part_data->getMinRho());
#endif
	if (colors[nvert_pos*4+3] > 1.0 || colors[nvert_pos*4+3] <0.0) {
	  std::cerr << "Erreur alpha = "<<colors[nvert_pos*4+3]<<"\n";
	}
      } 
      else {
        colors[nvert_pos*4+3] = 1.0;
      }
    } 
    else {
      colors[nvert_pos*4+3] = 1.0;
    }
    nvert_pos++;
  }
  vbobench.restart();
  // bind VBO in order to use
  //glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_color);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  //std::cerr << "vbo_pos = " << vbo_pos << "\n";
  assert( nvert_pos <= (po->npart/po->step)+1);
  // upload data to VBO
#if 0
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 4 * sizeof(float), colors, GL_STATIC_DRAW_ARB);
#else
 glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 3 * sizeof(float), nvert_pos * 4 * sizeof(float), colors);
#endif
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

  if (BENCH) qDebug("WALL Time elapsed to build VBO color: %f s", tbench.elapsed()/1000.);
  if (BENCH) qDebug("TRANSFERT Time elapsed to build VBO color: %f sec %f MB/s", vbobench.elapsed()/1000.,
      (nvert_pos * 4 * sizeof(float)/1024./1024.)/(vbobench.elapsed()/1000.));
  
  delete [] colors;
}
// ============================================================================
// buildVboSize                                                                
// Build Vector Buffer Object for size point array                             
void GLObjectParticles::buildVboSize()
{
  QTime tbench;
  tbench.restart();
  nvert_pos=0;
  
  GLfloat * size_neib = new GLfloat[((po->npart/po->step)+1)*3]; // create size_neib array
  QColor mycolor = po->getColor();
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    
    size_neib[nvert_pos*3+1] = 0.;
    size_neib[nvert_pos*3+2] = 0.;
    float log_rho;
    if (part_data->rho) {
      log_rho=((log(part_data->rho[index])-log(part_data->getMinRho()))/
	  (log(part_data->getMaxRho())-log(part_data->getMinRho())));
    }
    if (part_data->rneib) {
      if (part_data->rneib[index] != -1) { // Normalize rho for alpha color

	if (log_rho > 0.0000 ) {
	  size_neib[nvert_pos*3+0] = part_data->rneib[index];
	} else {
	  size_neib[nvert_pos*3+0] = 0.0;
	}
	
      } else {
	size_neib[nvert_pos*3+0] = 1.0;
      }
    } else {
      size_neib[nvert_pos*3+0] = 1.0;
    }
    size_neib[nvert_pos*3+0] *= 2.0;
    nvert_pos++;
  }
  // bind VBO in order to use
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
  //std::cerr << "vbo_pos = " << vbo_pos << "\n";
  assert( nvert_pos <= (po->npart/po->step)+1);
// upload data to VBO
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 3 * sizeof(float), size_neib, GL_STATIC_DRAW_ARB);
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  
  delete [] size_neib;
  if (BENCH) qDebug("Time elapsed to build VBO size: %f s", tbench.elapsed()/1000.);
}
// ============================================================================
// buildVboSize                                                                
// Build Vector Buffer Object for size point array                             
void GLObjectParticles::buildVboSize2()
{
  QTime tbench;
  tbench.restart();
  nvert_pos=0;
  //float LOGDMIN = log(DENS_MIN);
  //float LOGDMAX = log(DENS_MAX);
  //float DIFMM   = log(DENS_MAX)-log(DENS_MIN);
  //float INVDIFMM= 1./DIFMM;

  GLfloat * size_neib = new GLfloat[((po->npart/po->step)+1)]; // create size_neib array
  QColor mycolor = po->getColor();
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].index;
    else                index = po->index_tab[i];
    
    float log_rho;
    if (part_data->rho) {

	if (part_data->rho[index] >= DENS_MIN && part_data->rho[index] <= DENS_MAX) {
          //log_rho=((log(part_data->rho[index])-LOGDMIN))*INVDIFMM;
          log_rho=1.;
	}
        else {
           log_rho = -1.0;
        }
    }
    if (part_data->rneib) {
      if (part_data->rneib[index] != -1) { // Normalize rho for alpha color

        if (log_rho > 0.0000 ) {
          size_neib[nvert_pos] = part_data->rneib[index];
        } else {
          size_neib[nvert_pos] = 0.0;
        }
        
      } else {
        size_neib[nvert_pos] = 1.0;
      }
    } else {
      size_neib[nvert_pos] = 1.0;
    }
    size_neib[nvert_pos] *= 2.0;
    nvert_pos++;
  }
  // bind VBO in order to use
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
  //std::cerr << "vbo_pos = " << vbo_pos << "\n";
  assert( nvert_pos <= (po->npart/po->step)+1);
// upload data to VBO
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * sizeof(float), size_neib, GL_STATIC_DRAW_ARB);
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  
  delete [] size_neib;
  if (BENCH) qDebug("Time elapsed to build VBO size2: %f s", tbench.elapsed()/1000.);
}

// ============================================================================
//
void GLObjectParticles::checkVboAllocation(const int sizebuf)
{
  GLint size;
  glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size); 
  if (size != sizebuf) {
    std::cerr << "FATAL: Pixel Buffer Object allocation failed!";
    std::cerr << "Expected size ["<<sizebuf<<"] != from "<<size<<"\n";
    exit(1);
  }
  //glBindBuffer(GL_ARRAY_BUFFER, 0);
}
// ============================================================================
//                     - Texture and Sprites management -                      
// ============================================================================
// ============================================================================

// ============================================================================
// setTexture()                                                                 
void GLObjectParticles::setTexture(QString name)
{
  assert(0); // we should not go here
  if (! texture) texture = new GLTexture();
  if (texture->load(name,NULL)) {
  } else {
    //QString message=tr("Unable to load texture");
    //QMessageBox::information( this,tr("Warning"),message,"Ok");
    std::cerr << "\n\nUnable to load TEXTURE.....\n\n";
  }
}
// ============================================================================
// setTexture()                                                                 
void GLObjectParticles::setTexture(const int tex)
{
  //assert(tex<3);
  assert(tex<(int)gtv->size());
//  std::cerr << "GLObjectParticles::setTexture["<<tex<<"] = <"<<GLTexture::TEXTURE[tex].toStdString()<<">\n";
//  std::cerr << "GLObjectParticles::setTexture["<<tex<<"] = <"<<((*gtv)[tex].getName()).toStdString()<<">\n";
  //setTexture(GLTexture::TEXTURE[tex]);
  texture = &(*gtv)[tex];
}
// ============================================================================
// setTexture()                                                                 
void GLObjectParticles::setTexture()
{
/*  if (texture) {
    QString name=texture->getName();
    delete texture;
    texture=NULL;
    setTexture(name);
  }*/
  if (!texture) {
    setTexture(0);
  }
}
// ============================================================================
// displayVboPoints                                                            
void GLObjectParticles::displayVboPoints()
{
  GLObject::setColor(po->getColor()); // set the color
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),po->getPartAlpha());
  if (go->render_mode == 0 || go->render_mode == 1) { // Alpha blending accumulation
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
  }
  else
    if (go->render_mode == 2) {  // No Alpha bending accumulation
    glDepthMask(GL_FALSE);
    glDisable(GL_DEPTH_TEST);
      //glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0.0f);
  }

 glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
    // Send vertex object positions
  glEnableClientState(GL_VERTEX_ARRAY); 

  glVertexPointer(3, GL_FLOAT, 0, 0);

  if (go->render_mode == 1 || go->render_mode == 2) { // individual size and color

    // Send vertex object colors
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_FLOAT, 0, (void *) (3*nvert_pos*sizeof(float)));
  }
  // Draw points
  //glPointSize(po->getGazSize()*2.0);
  //GLObject::setColor(po->getColor()); // set the color
#if GLDRAWARRAYS
  glDrawArrays(GL_POINTS, 0, nvert_pos);
#else
#if 0
  glDrawElements(GL_POINTS,nind_sorted,GL_UNSIGNED_INT, indexes_sorted);
  
#else
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,vbo_index);
  glDrawRangeElements(GL_POINTS,0,nind_sorted,nind_sorted,GL_UNSIGNED_INT,0);
#endif
#endif
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);


  if (go->render_mode == 1 || go->render_mode == 2) {
    glDisableClientState(GL_COLOR_ARRAY);
  }
  glDisableClientState(GL_VERTEX_ARRAY);
  
}
// ============================================================================
// displayVboSprites()                                                            
void GLObjectParticles::displayVboSprites(int win_height,const bool front)
{
  static bool zsort=false;
  if (front) {;};
  //  detect if rho exist for the component
  int index;
  bool is_rho=false;
  if (part_data->rho) {
    index = rho[0].index;
    if (part_data->rho[index] != -1) is_rho = true;
  }
  if (go->zsort) { // Z sort particles
      zsort = true;
      sortByDepth();
  } else {
      if (zsort) {
        zsort = false;
        sortByDensity(); // we have to resort by density
      }
  }
  // setup point sprites
  glEnable(GL_POINT_SPRITE_ARB);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
  //float quadratic[] =  { 0.0f, 0.0f, 0.01f };
  //glPointParameterfvARB( GL_POINT_DISTANCE_ATTENUATION_ARB, quadratic );
  glPointSize(po->getGazSize()*win_height);
  GLObject::setColor(po->getColor()); // set the color 
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),po->getGazAlpha());
  //glNormal3f(2.0,1.,1.);
//   GLfloat lightpos[] = {0., 0., 0., 0.};
//   glEnable(GL_LIGHTING);
//   glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
//   glEnable(GL_LIGHT0);
  //glEnable(GL_COLOR_MATERIAL);
  int err;
  if (go->render_mode == 0 || go->render_mode == 1) { // Alpha blending accumulation	
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
    if ((err = glGetError())) { fprintf(stderr,"O et 1 c error %x\n", (unsigned int)err); }
  }
  else 
    if (go->render_mode == 2) {  // No Alpha bending accumulation
#if 1
      glDepthMask(GL_FALSE);
      glDisable(GL_DEPTH_TEST);
      //glEnable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      //glBlendFunc (GL_ONE, GL_ONE);
      //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
      glEnable(GL_ALPHA_TEST);
      glAlphaFunc(GL_GREATER, 0.00f);

      if ((err = glGetError())) { fprintf(stderr,"2 c error %x\n", (unsigned int)err); }
#else
 glEnable( GL_DEPTH_TEST );
 glEnable( GL_MULTISAMPLE_ARB );
 glEnable( GL_SAMPLE_ALPHA_TO_COVERAGE_ARB );
#endif
  }

  glTexEnvi(GL_POINT_SPRITE,GL_COORD_REPLACE,GL_TRUE); 
  glUseProgramObjectARB(GLWindow::m_program);
  if ((err = glGetError())) { fprintf(stderr,"c error %x\n", (unsigned int)err); }
#if 0    
  // Send windows'height data to Vertex Shader
  int h_loc  = glGetUniformLocation(GLWindow::m_program, "win_height");
  if (h_loc == -1) {
    std::cerr << "Error occured when sending \"win_height\" to Vertex shader..\n";
    exit(1);
  }
  glUniform1i(h_loc, win_height); // send windows height
#endif
  // Send alpha channel color to Vertex Shader
  int alpha_loc  = glGetUniformLocation(GLWindow::m_program, "alpha");
  if (alpha_loc == -1) {
    std::cerr << "Error occured when sending \"alpha\" to Vertex shader..\n";
    exit(1);
  }
  float alpha=po->getGazAlpha()/255.;
  if (go->render_mode == 1 || go->render_mode == 2) { // individual size and color
    if (is_rho)  
      glUniform1f(alpha_loc, alpha*alpha); // send alpha channel
    else
      glUniform1f(alpha_loc, alpha); // send alpha channel
  }
  if (go->render_mode == 0 ) { // global size and color
    glUniform1f(alpha_loc, 1.); // send alpha channel
  }

  // Send data to Pixel Shader
  int tex_loc = glGetUniformLocation(GLWindow::m_program, "splatTexture");
  if (tex_loc == -1) {
    std::cerr << "Error occured when sending \"splatTexture\" to Pixel shader..\n";
    exit(1);
  }
  glUniform1i(tex_loc, 0);       // send texture index
  
  glActiveTextureARB(GL_TEXTURE0_ARB);
  texture->glBindTexture();  // bind texture
    //glBindTexture(GL_TEXTURE_2D, m_texture);

  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  // Send vertex object positions
  glEnableClientState(GL_VERTEX_ARRAY);
  //glEnableClientState(GL_COLOR_ARRAY);                
  
  glVertexPointer(3, GL_FLOAT, 0, 0);

  // get attribute location for sprite size
  int a_sprite_size = glGetAttribLocation(GLWindow::m_program, "a_sprite_size");
  glVertexAttrib1f(a_sprite_size,1.0);
  if ( a_sprite_size == -1) {
    std::cerr << "Error occured when getting \"a_sprite_size\" attribute\n";
    exit(1);
  }

  if (go->render_mode == 1 || go->render_mode == 2) { // individual size and color
  
    // Send vertex object colors
    
    //glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_color);
    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(4, GL_FLOAT, 0, (void *) (3*nvert_pos*sizeof(float)));
    // Send vertex object neighbours size
#if 1
    glEnableVertexAttribArray(a_sprite_size);
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
    glVertexAttribPointer(a_sprite_size,1,GL_FLOAT, 0, 0, 0);
#endif
#if 0
    glEnableClientState(GL_NORMAL_ARRAY);                
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
    glNormalPointer(GL_FLOAT,  3*sizeof(float), 0 );
#endif
  }
  // Draw points 
#if GLDRAWARRAYS
  glDrawArrays(GL_POINTS, 0, nvert_pos);
#else
#if 0
  glDrawElements(GL_POINTS,nind_sorted,GL_UNSIGNED_INT, indexes_sorted);
#else
  //glEnableClientState(GL_VERTEX_ARRAY);
  int a,b;
  //glGetIntegerv(GL_MAX_ELEMENTS_VERTICES,&a);
  //glGetIntegerv(GL_MAX_ELEMENTS_INDICES,&b);
  //std::cerr << "Max vert=" << a << "\n";
  //std::cerr << "Max inde=" << b << "\n";
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,vbo_index);
  glDrawElements(GL_POINTS,-1+nind_sorted/2,GL_UNSIGNED_INT,0);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,vbo_index2);
  glDrawElements(GL_POINTS,-1+nind_sorted/2,GL_UNSIGNED_INT,0);

  //glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
  //glDrawRangeElements(GL_POINTS,0,nind_sorted,nind_sorted,GL_UNSIGNED_INT,0);
//glDrawRangeElements(GL_POINTS,0,index_max-index_min,index_max-index_min,GL_UNSIGNED_INT,0);
#endif
#endif
  //glDrawArrays(GL_POINTS, 0, nvert_pos);
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

  glUseProgram(0); // deactivate shaders programs


  if (go->render_mode == 1 || go->render_mode == 2) {
    //glDisableClientState(GL_NORMAL_ARRAY);
    glDisableVertexAttribArray(a_sprite_size);
    glDisableClientState(GL_COLOR_ARRAY);
    
  }
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisable(GL_POINT_SPRITE_ARB);
  glDisable(GL_BLEND);
  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);

}
// ============================================================================
// sortDyDensity                                                               
void GLObjectParticles::sortByDensity()
{
      // sort according to the density
  if (part_data->rho)  {
    sort(rho.begin(),rho.end(),GLObjectIndexTab::compareLow);
  }
  nind_sorted = 0;
    // creates vertex indices array for gLDrawElements
  if (! indexes_sorted || ((int) nind_sorted) < po->npart) {
    if (indexes_sorted) delete [] indexes_sorted;
    indexes_sorted = new GLuint[po->npart];
    nind_sorted = 0;
  }
  index_min = 0;
  index_max = po->npart;
  int cpt=0;
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = rho[i].i_point;
    
    if (part_data->rho) {
        index = rho[i].i_point;
        int index2 = rho[i].index;
    	if (part_data->rho[index2] >= DENS_MIN && part_data->rho[index2] <= DENS_MAX) {
        indexes_sorted[cpt++] = index;
	}
    } else {
      indexes_sorted[cpt++] = i;
    }
  }
nind_sorted=cpt;
  std::cerr << "index_min="<< index_min <<"   index_max=" << index_max << "\n";
      // bind VBO in order to use
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vbo_index);
    // upload data to VBO
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, (nind_sorted/2)*sizeof(GLuint), indexes_sorted, GL_STATIC_DRAW_ARB);
    //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);

  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vbo_index2);
    // upload data to VBO
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, (nind_sorted/2)*sizeof(GLuint), indexes_sorted+nind_sorted/2, GL_STATIC_DRAW_ARB);
    //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);

}
// ============================================================================
// sortDyDepth                                                                 
void GLObjectParticles::sortByDepth()
{
#define MM(row,col)  mModel[col*4+row]
  GLdouble mModel[16];
  // get ModelView Matrix
  glGetDoublev(GL_MODELVIEW_MATRIX, (GLdouble *) mModel);
  nind_sorted=0;
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    if (part_data->rho[index] >= DENS_MIN && part_data->rho[index] <= DENS_MAX) {
    float
        x=part_data->pos[index*3  ],
        y=part_data->pos[index*3+1],
        z=part_data->pos[index*3+2];
         // compute point coordinates according to model view via matrix
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    zdepth[nind_sorted].value=mz;
    zdepth[nind_sorted].i_point=i;
    zdepth[nind_sorted].index = index;
    nind_sorted++;
   }
  }
  // sort according to the Z Depth
  sort(zdepth.begin(),zdepth.begin()+nind_sorted,GLObjectIndexTab::compareLow);
  //nind_sorted = 0;
  // creates vertex indices array for gLDrawElements
  if (! indexes_sorted || nind_sorted < zdepth.size()) {
    if (indexes_sorted) delete [] indexes_sorted;
    indexes_sorted = new GLuint[zdepth.size()];
    //nind_sorted = 0;
  }
  for (unsigned int i=0; i < nind_sorted; i++) {
    int index = zdepth[i].i_point;
    indexes_sorted[i] = index;
  }
      // bind VBO in order to use
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vbo_index);
    // upload data to VBO
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, (nind_sorted)*sizeof(GLuint), indexes_sorted, GL_STATIC_DRAW_ARB);
    //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
}
// ============================================================================
// displaySprites()                                                            
// use billboarding technique to display sprites :                             
// - transforms coordinates points according model view matrix                 
// - draw quad (2 triangles) around new coordinates and facing camera          
void GLObjectParticles::displaySprites(const double * mModel)
{
#define MM(row,col)  mModel[col*4+row]
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  glPushMatrix();
  glLoadIdentity();             // reset opengl state machine
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  glEnable( GL_TEXTURE_2D );
  glEnable(GL_BLEND);
#if 1
  glDisable(GL_DEPTH_TEST);
#else
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0);
#endif
  
  GLObject::setColor(po->getColor()); // set the color 
  glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),po->getGazAlpha());
  
  texture->glBindTexture();  // bind texture
  // uv coordinates
  float uv[4][2] = { {0.0         , 1.0-texture->V()}, {0.0         , 1.0             },
                     {texture->U(), 1.0             }, {texture->U(), 1.0-texture->V()}
                   };

  // some usefull variables
  float rot=0.0;
  float new_texture_size = po->getGazSize()/2.;
  // loop on all particles
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    float
      x=part_data->pos[index*3  ],
      y=part_data->pos[index*3+1],
      z=part_data->pos[index*3+2];
    // compute point coordinates according to model via matrix
    float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3);//*w;
    float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3);//*w;
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    int modulo=2;
    rot=index;

    glPushMatrix();
    glTranslatef(mx,my,mz);        // move to the transform particles
    if (po->isGazRotate())
        glRotatef(rot,0.0,0.0,1.0);// rotate triangles around z axis
    glBegin(GL_TRIANGLES);

    // 1st triangle
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    glVertex3f(new_texture_size , new_texture_size  ,0. );
    modulo = (modulo+1)%4;
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    glVertex3f(new_texture_size , -new_texture_size  ,0. );
    modulo = (modulo+1)%4;
    glTexCoord2f(uv[modulo][0],   uv[modulo][1]);
    glVertex3f(-new_texture_size , -new_texture_size  ,0. );
    // second triangle
    modulo = (modulo+2)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(new_texture_size , new_texture_size  ,0. );
    modulo = (modulo+2)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(-new_texture_size , -new_texture_size  ,0. );
    modulo = (modulo+1)%4;
    glTexCoord2f(uv[modulo][0],  uv[modulo][1]);
    glVertex3f(-new_texture_size , new_texture_size  ,0. );

    glEnd();
    glPopMatrix();
  }
  glDisable(GL_BLEND);
  glDisable( GL_TEXTURE_2D );
  glPopMatrix();
}


}
