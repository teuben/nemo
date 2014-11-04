// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
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
#ifdef _OPENMP
#include <parallel/algorithm>
#include <omp.h>
#endif


#define OLDRENDER 0
#define BENCH 1
#define GLDRAWARRAYS 1
namespace glnemo {
float PHYS_MIN=-1.;
float PHYS_MAX=0.000006;
int index_min, index_max;
int nhisto=10000;

// ============================================================================
// constructor                                                                 
GLObjectParticles::GLObjectParticles(GLTextureVector * _gtv ):GLObject()
{
  dplist_index = glGenLists( 1 );    // get a new display list index
  texture = NULL;                    // no texture yet              
  gtv = _gtv;
  // reserve memory
  index_histo.reserve(nhisto);
  if (GLWindow::GLSL_support) {
    glGenBuffersARB(1,&vbo_pos);
    glGenBuffersARB(1,&vbo_data);
//    glGenBuffersARB(1,&vbo_color);
    glGenBuffersARB(1,&vbo_size);
    glGenBuffersARB(1,&vbo_index);
    glGenBuffersARB(1,&vbo_index2);
  }
  indexes_sorted = NULL;
  nind_sorted = 0;
  hasPhysic = false;
}
// ============================================================================
// constructor                                                                 
GLObjectParticles::GLObjectParticles(const ParticlesData   * _part_data,
                                     ParticlesObject * _po,
                                     const GlobalOptions   * _go,
                                     GLTextureVector * _gtv, CShader * _shader):GLObject()
{
  shader = _shader; // link shader program pointer
  dplist_index = glGenLists( 1 );    // get a new display list index
  vel_dp_list  = glGenLists( 1 );    // get a new display vel list
  orb_dp_list  = glGenLists( 1 );    // get a new display orb list
  // reserve memory
  index_histo.reserve(nhisto);
  if (GLWindow::GLSL_support) {
    glGenBuffersARB(1,&vbo_pos);     // get Vertex Buffer Object
//    glGenBuffersARB(1,&vbo_color);   // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_size);    // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_index);   // get Vertex Buffer Object
    glGenBuffersARB(1,&vbo_index2);
    glGenBuffersARB(1,&vbo_data);    
  }
  indexes_sorted = NULL;
  nind_sorted = 0;
  texture = NULL;                    // no texture yet
  gtv = _gtv;
  assert(gtv->size()>0);
  setTexture(0);
  hasPhysic = false;
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
      if (GLWindow::GLSL_support) 
        displayVboShader(win_height,true);
      else {
        glEnable(GL_BLEND);
        glEnable(GL_POINT_SMOOTH);
        glPointSize((float) po->getPartSize());
        GLObject::updateAlphaSlot(po->getPartAlpha());
        GLObject::setColor(po->getColor());
        GLObject::display();
        glDisable(GL_BLEND);
      }
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
        displayVboShader(win_height,false);      
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
// displayVboShader()                                                            
void GLObjectParticles::displayVboShader(const int win_height, const bool use_point)
{
  static bool zsort=false;
  int err;

  //  detect if rho exist for the component
  int index;
  bool is_rho=false;
  if (phys_select && phys_select->isValid()) {
    index = phys_itv[0].index;
    if (phys_select->data[index] != -1) is_rho = true;
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
  checkGlError("GLObjectParticles::displayVboShader -> Beginning");

  // setup point sprites
  glEnable(GL_POINT_SPRITE_ARB);
  glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);

  if (hasPhysic && go->render_mode==1) {                 // if physic
    GLObject::setColor(Qt::black); // send black color to the shader
  } else {                              // else
    GLObject::setColor(po->getColor()); // send user selected color 
  }
  if (use_point)
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),po->getPartAlpha());
  else
    glColor4ub(mycolor.red(), mycolor.green(), mycolor.blue(),po->getGazAlpha());
  
  //if ((go->render_mode == 0 || go->render_mode == 1) && !hasPhysic) { // Alpha blending accumulation
  if ((go->render_mode == 0 ) ) { // Alpha blending accumulation
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);
    checkGlError("GLObjectParticles::displayVboShader -> Alpha blending accumulation");
  }
  else 
    if (go->render_mode == 1) {  // No Alpha bending accumulation
      glDepthMask(GL_FALSE);
      glDisable(GL_DEPTH_TEST);
      //glEnable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      //glBlendFunc (GL_ONE, GL_ONE);
      //glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE);
      glEnable(GL_ALPHA_TEST);
      glAlphaFunc(GL_GREATER, 0.00f);

      if ((err = glGetError())) { 
        fprintf(stderr,">> 2 c error %x\n", (unsigned int)err); 
      }
  }

  glTexEnvi(GL_POINT_SPRITE,GL_COORD_REPLACE,GL_TRUE); 

  // start shader program
  shader->start();
  
  // process shader color variables
  sendShaderColor(win_height,use_point);  
  
  glActiveTextureARB(GL_TEXTURE0_ARB);
  texture->glBindTexture();  // bind texture
  
  // Send vertex object positions
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  glEnableClientState(GL_VERTEX_ARRAY);
  int start=3*min_index*sizeof(float);
  int maxvert=max_index-min_index+1;
  //std::cerr << "min_index="<<min_index<<" max_index="<<max_index<<" maxvert="<<maxvert<<"\n";
  glVertexPointer(3, GL_FLOAT, 0, (void *) (start));

  // get attribute location for sprite size
  int a_sprite_size = glGetAttribLocation(shader->getProgramId(), "a_sprite_size");
  glVertexAttrib1f(a_sprite_size,1.0);
  if ( a_sprite_size == -1) {
    std::cerr << "Error occured when getting \"a_sprite_size\" attribute\n";
    exit(1);
  }

  // get attribute location for phys data
  int a_phys_data = glGetAttribLocation(shader->getProgramId(), "a_phys_data");
  glVertexAttrib1f(a_phys_data,1.0);
  if ( a_sprite_size == -1) {
    std::cerr << "Error occured when getting \"a_phys_data\" attribute\n";
    exit(1);
  }
  if ((go->render_mode == 1 )) { // individual size and color
    // Send vertex object neighbours size
    if (hasPhysic && phys_select && phys_select->isValid()) { 
      // set back texture_size to one for gas
      //po->setGazSize(1.0);
      //po->setGazSizeMax(1.0);
      glEnableVertexAttribArrayARB(a_sprite_size);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
      start = min_index*sizeof(float);
      glVertexAttribPointerARB(a_sprite_size,1,GL_FLOAT, 0, 0, (void *) (start));
    }
    // Send physical data
    if (hasPhysic && phys_select && phys_select->isValid()) {  
      glEnableVertexAttribArrayARB(a_phys_data);
      glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_data);      
      start = min_index*sizeof(float);
      glVertexAttribPointerARB(a_phys_data,1,GL_FLOAT, 0, 0, (void *) (start));
    }
  } else {
    if (hasPhysic) { // gas only
      //glVertexAttrib1f(a_sprite_size,go->texture_size);
    }
  }
  // Draw points 
#if GLDRAWARRAYS
#if 0
  std::cerr << " hasPhysic ? =" << hasPhysic <<"\n";
  std::cerr << " min_index = " <<min_index << "  max_index = " << max_index << "\n";
  std::cerr <<"maxvert="<<maxvert<< " #part="<<max_index-min_index+1<< " nvert_pos ="<<nvert_pos<<"\n";
#endif
  if (maxvert > 0 && maxvert<=nvert_pos) {
    //std::cerr << ">> rendering...\n";
    glDrawArrays(GL_POINTS, 0, maxvert);
    //std::cerr << "<< rendering...\n";
  }

#else
#if 1

  // marche pas....
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,vbo_index);
  glDrawElements(GL_POINTS,nind_sorted,GL_UNSIGNED_INT,0);
  // marche pas....

  //glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,vbo_index);
  //glDrawRangeElements(GL_POINTS,0,nind_sorted,nind_sorted,GL_UNSIGNED_INT,0);


  // >>> works alone
  //glDrawElements(GL_POINTS,nind_sorted,GL_UNSIGNED_INT, indexes_sorted);
  // <<< works alone
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

  // deactivate shaders programs
  shader->stop();


  if (hasPhysic && ( go->render_mode == 1)) {
    //glDisableClientState(GL_NORMAL_ARRAY);
    if (phys_select && phys_select->isValid()) {  
      glDisableVertexAttribArray(a_sprite_size);
      glDisableVertexAttribArray(a_phys_data);
    }
    //glDisableClientState(GL_COLOR_ARRAY);
    
  }
  glDisableClientState(GL_VERTEX_ARRAY);
  glDisable(GL_POINT_SPRITE_ARB);
  glDisable(GL_BLEND);
  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);

}
// ============================================================================
// update                                                                      
void GLObjectParticles::update( const ParticlesData   * _part_data,
                                ParticlesObject * _po,
                                const GlobalOptions   * _go,
                                const bool update_obj)
{
  // update variables
  part_data = _part_data;
  po        = _po;
  go        = _go;
  min_index = 0;
  // get physical value data array
  phys_select = part_data->getPhysData();
  phys_select_id = part_data->getIpvs();
  hasPhysic = po->hasPhysic();//checkHasPhysic(); // check the object has physic            
  // color
  mycolor   = po->getColor();
  
  if (update_obj) { // force to rebuild VBO and display list
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
      checkGlError("GLObjectParticles::update buildVboPos");
      buildVboPhysData();
      checkGlError("GLObjectParticles::update buildPhysData");
      buildVboHsml();
      checkGlError("GLObjectParticles::update buildVboSize2");
      
#if ! GLDRAWARRAYS
      sortByDensity();
#endif
      updateColormap();
      checkGlError("GLObjectParticles::update updateColormap");
      updateBoundaryPhys();
      checkGlError("GLObjectParticles::update updateBoundaryPhys");
    }
  }
}
// ============================================================================
// updateVbo                                                                   
void GLObjectParticles::updateVbo()
{
  // get physical value data array

  //if (phys_select != part_data->getPhysData()) { // new physical quantity
  if (phys_select_id != part_data->getIpvs()) { // new physical quantity
    phys_select=part_data->getPhysData();
    phys_select_id=part_data->getIpvs();
    hasPhysic = po->hasPhysic();//checkHasPhysic(); // check the object has physic      
    buildVboPos();
  }
  buildVboPhysData();
  if (GLWindow::GLSL_support) { // && po->isGazEnable()) {

    if (go->render_mode == 1) {
      // We must re ordering indexes
      if (phys_select && phys_select->isValid()) {
#if ! GLDRAWARRAYS
        phys_itv.clear();
#endif
        for (int i=0; i < po->npart; i+=po->step) {         
          if (phys_select && phys_select->isValid()) {
#if ! GLDRAWARRAYS
            int index=po->index_tab[i];
            GLObjectIndexTab myrho;
            myrho.index   = index;
            myrho.value   = phys_select->data[index];
            myrho.i_point = i;
            phys_itv.push_back(myrho);
#endif
          }
        }
      }
#if ! GLDRAWARRAYS
      sortByDensity();	
#else // GLDRAWARRAYS
      buildVboHsml();
#endif
      updateBoundaryPhys();
    }
  }
}
// ============================================================================
// updateVbo                                                                   
void GLObjectParticles::updateColorVbo()
{
  if (1||go->render_mode == 1 ) {
    updateBoundaryPhys();
  }
}
// ============================================================================
// updateBoundaryPhys();
// according to slide value from Option Dialog Box, this function
// compute in which range (min,max) particles are selected
void GLObjectParticles::updateBoundaryPhys()
{
  min_index = 0;
  max_index = nvert_pos-1;
  if (hasPhysic && phys_select && phys_select->isValid()) {
    //std::cerr << " Pobj min index="<<po->getMinPercenPhys()
    //    << " max index="<<po->getMaxPercenPhys()<<"\n";
    int permin=po->getMinPercenPhys();
    int permax=po->getMaxPercenPhys();
    assert(permin>=0 && permax<100);
    int imin=permin*nhisto/100;
    min_index = index_histo[imin];
    int imax =permax*nhisto/100;
    max_index = index_histo[imax];
    if (permax==99) {
        max_index = nvert_pos-1; // we take all the vertex
    }
#if 0
    std::cerr << "min ="<<min_index << " max ="<<max_index<<"\n";
    std::cerr << "imin ="<<imin << " imax ="<<imax<<"\n";
    std::cerr << "permin ="<<permin << " permax ="<<permax<<"\n";
#endif
  }

}
// ============================================================================
// buildVboPos                                                                 
// Build Vector Buffer Object for positions array                              
void GLObjectParticles::buildVboPos()
{
  QTime tbench,tbloc;
  tbench.restart();
  nvert_pos=0;
  tbloc.restart();
  
  std::vector <GLfloat> vertices;
  vertices.reserve(((po->npart/po->step)+1)*3);
  
  //rho.clear();      // clear rho density vector
  phys_itv.clear(); // clear ohysical value vector
  rho_itv.clear();
  vindex_sel.clear();   // clear zdepth vector     
  rho_itv.reserve(((po->npart/po->step)+1));
  vindex_sel.reserve(((po->npart/po->step)+1));
  phys_itv.reserve(((po->npart/po->step)+1));
  
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    if (phys_select && phys_select->isValid()) {
      GLObjectIndexTab myphys;
      myphys.index   = index;
      myphys.value   = phys_select->data[index];
      myphys.i_point = i;
      phys_itv.push_back(myphys);
      
      if (po->rhoSorted() && 
          phys_select->getType() != PhysicalData::rho && part_data->rho) {
        GLObjectIndexTab myphys;
        myphys.index   = index;
        myphys.value   = part_data->rho->data[index];
        myphys.i_point = i;
        rho_itv.push_back(myphys);
      }
    }
        
#if 1 // used if z depth test activated
    GLObjectIndexTab myz;
    myz.index = index;
    myz.i_point = i;
    vindex_sel.push_back(myz);
#endif
  }
  if (BENCH) qDebug("Time elapsed to setup PHYSICAL arrays: %f s", tbloc.elapsed()/1000.);
  
#ifdef _OPENMP
  int ntask=omp_get_max_threads(); // return number of task requested
  std::cerr << "#OpenMP TASKS = "<<ntask<<"\n";  
#endif
  // sort by density
#if GLDRAWARRAYS
  tbloc.restart();
  if (po->rhoSorted() &&
      phys_select && phys_select->getType() != PhysicalData::rho && part_data->rho) {
    sort(rho_itv.begin(),rho_itv.end(),GLObjectIndexTab::compareLow);
  } else {
    sort(phys_itv.begin(),phys_itv.end(),GLObjectIndexTab::compareLow);
  }
  if (BENCH) qDebug("Time elapsed to SORT PHYSICAL arrays: %f s", tbloc.elapsed()/1000.);
  //sort(rho.begin(),rho.end(),GLObjectIndexTab::compareHigh);
#endif
  // select vertices
  tbloc.restart();
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (po->rhoSorted() &&
        phys_select && phys_select->getType() != PhysicalData::rho && part_data->rho) {
      index = rho_itv[i].index; // it's temperature/pressure, we sort by density
    }
    else { 
      if (phys_select && phys_select->isValid()) 
        index = phys_itv[i].index; // we sort by physical value
      else                index = po->index_tab[i]; // no physic
    }
    // fill vertices array sorted by density
    vertices.push_back(part_data->pos[index*3  ]);
    vertices.push_back(part_data->pos[index*3+1]);
    vertices.push_back(part_data->pos[index*3+2]);
    nvert_pos++; // one more particle
  }
  if (BENCH) qDebug("Time elapsed to setup VBO arrays: %f s", tbloc.elapsed()/1000.);
  // build first particle index in the histo
  min_index=0; max_index=nvert_pos-1;

  tbloc.restart();
  buildIndexHisto();  
  if (BENCH) qDebug("Time elapsed to build indexes histo : %f s", tbloc.elapsed()/1000.);

  tbloc.restart();
  // bind VBO buffer for sending data
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_pos);
  assert( nvert_pos <= (po->npart/po->step)+1);
  std::cerr << "buildVbo Pos nvert_pos="<<nvert_pos<<"\n";
  
  // upload data to VBO
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, nvert_pos * 3 * sizeof(float), &vertices[0], GL_STATIC_DRAW_ARB);

  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  if (BENCH) qDebug("Time elapsed to transfert VBO arrays to GPU: %f s", tbloc.elapsed()/1000.);
  vertices.clear();
  if (BENCH) qDebug("Time elapsed to build VBO pos: %f s", tbench.elapsed()/1000.);
}
// ============================================================================
// buildVboHsml                                                             
// Build Vector Buffer Object for size point array                             
void GLObjectParticles::buildVboHsml()
{
  QTime tbench;
  tbench.restart();
  
  std::vector <GLfloat> hsml_value;
  hsml_value.reserve(((po->npart/po->step)+1));

  // loop on all the object's particles
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
#if 0
    if (phys_select && phys_select->isValid()) index = phys_itv[i].index;
    else                index = po->index_tab[i];
#endif
    if (po->rhoSorted() &&
        phys_select && phys_select->getType() != PhysicalData::rho && part_data->rho) {
      index = rho_itv[i].index; // it's temperature/pressure, we sort by density
    }
    else { 
      if (phys_select && phys_select->isValid()) 
        index = phys_itv[i].index; // we sort by physical value
      else                index = po->index_tab[i]; // no physic
    }
    if (part_data->rneib) {
      if (part_data->rneib->data[index] != -1) {
        hsml_value.push_back(2.0*part_data->rneib->data[index]);
      } else {
        hsml_value.push_back(2.0);      
      }
    }    
  }
  
  // bind VBO buffer for sending data
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_size);
  
  assert( (int) hsml_value.size() <= (po->npart/po->step)+1);
  // upload data to VBO
  glBufferDataARB(GL_ARRAY_BUFFER_ARB, hsml_value.size() * sizeof(float), &hsml_value[0], GL_STATIC_DRAW_ARB);
  //checkVboAllocation((int) (nvert_pos * 3 * sizeof(float)));
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  std::cerr << "buildVboHsml ="<<hsml_value.size()<<"\n";

  hsml_value.clear();
  if (BENCH) qDebug("Time elapsed to build VBO Hsml: %f s", tbench.elapsed()/1000.);
}
// ============================================================================
// buildVboData                                                                
// Build Vector Buffer Object for physical data  array                             
void GLObjectParticles::buildVboPhysData()
{
  QTime tbench;
  tbench.restart();

  std::vector <GLfloat> phys_data;
  phys_data.reserve(((po->npart/po->step)+1));

  // loop on all the object's particles
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
#if 0
    if (phys_select && phys_select->isValid()) index = phys_itv[i].index;
    else                index = po->index_tab[i];
#endif
    if (po->rhoSorted() &&
        phys_select && phys_select->getType() != PhysicalData::rho && part_data->rho) {
      index = rho_itv[i].index; // it's temperature/pressure, we sort by density
    }
    else { 
      if (phys_select && phys_select->isValid()) 
        index = phys_itv[i].index; // we sort by physical value
      else                index = po->index_tab[i]; // no physic
    }
    if (phys_select && phys_select->isValid()) {
      phys_data.push_back(phys_select->data[index]);
    } else {
      phys_data.push_back(0.0);      
    }
  }
  
  // bind VBO buffer for sending data
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, vbo_data);
  
  assert( (int)phys_data.size() <= (po->npart/po->step)+1);
  // upload data to VBO
  glBufferDataARB(GL_ARRAY_BUFFER_ARB,phys_data.size() * sizeof(float), &phys_data[0], GL_STATIC_DRAW_ARB);
  checkGlError("2222");
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
  checkGlError("3333");
  std::cerr << "Phys_data size="<<phys_data.size()<<"\n";
  phys_data.clear();
  //delete [] phys_data;
  if (BENCH) qDebug("Time elapsed to build VBO data: %f s", tbench.elapsed()/1000.);
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
// updateColormap();                                                                   
void GLObjectParticles::updateColormap()
{
  cmap.clear();
  cmap.reserve(go->R->size()*3);
  float powcolor=go->powercolor;
  for (unsigned int i=0; i<go->R->size(); i++) {    
      cmap.push_back(pow((*go->R)[i],powcolor));
      cmap.push_back(pow((*go->G)[i],powcolor));
      cmap.push_back(pow((*go->B)[i],powcolor));
  }
  //std::cerr << " cmap size="<<cmap.size()/3<<"\n";
}
// ============================================================================
// sendShaderColor();
void GLObjectParticles::sendShaderColor(const int win_height, const bool use_point)
{  
  static bool first=true;
  bool physic=(phys_select && phys_select->isValid())?true:false;
  
  if (first && physic ) {
    first=false;
    //po->setMaxPhys(phys_select->getMax());
    //po->setMinPhys(phys_select->getMin());
  }
  
  float alpha;
  if (use_point) alpha=po->getPartAlpha()/255.;
  else alpha=po->getGazAlpha()/255.;
  if (go->render_mode == 1 ) { // individual size and color
    if (physic)  
      shader->sendUniformf("alpha", alpha*alpha); // send alpha channel
    else
      shader->sendUniformf("alpha", alpha); // send alpha channel      
  }
  if (go->render_mode == 0 ) { // global size and color
    shader->sendUniformf("alpha", 1.0); // send alpha channel      
  }

  // Send texture size factor
  if (use_point) {
    shader->sendUniformf("factor_size",(float) po->getPartSize());
  } else {
    if (go->perspective) {
        shader->sendUniformf("factor_size",po->getGazSize()*win_height);
    } else {
      shader->sendUniformf("factor_size",po->getGazSize()*win_height/fabs(go->zoom));//*win_height);
    }
  }
  if (go->perspective) // perspective  projection
    shader->sendUniformf("zoom",(float) go->zoom);
  else                 // orthographic projection
    shader->sendUniformf("zoom",(float) 0.0);
  // send zoom
  if (go->od_enable)
    shader->sendUniformi("show_zneg",(int) 0);
  else
    shader->sendUniformi("show_zneg",(int) 1);

  // opaque disc radius
  shader->sendUniformf("radius",go->od_radius);

  // coronographe
  shader->sendUniformi("coronograph",(int) go->od_display);

  // Send data to Pixel Shader
  shader->sendUniformi("splatTexture",0);

  int viewporti[4];
  float viewport[4];
  glGetIntegerv(GL_VIEWPORT,viewporti);
  for (int i=0; i<4; i++) {
      viewport[i] = viewporti[i];
  }
  shader->sendUniformXfv("viewport",1,4,&viewport[0]);

  // use point or texture ?
  shader->sendUniformi("use_point",use_point);
  
  // perspective mode ?
  shader->sendUniformi("perspective",go->perspective);
  
  // send colormap stuffs
  shader->sendUniformXfv("colormap",3,cmap.size()/3,&cmap[0]);
  shader->sendUniformi("ncmap",cmap.size()/3);
  //std::cerr << "ncmap ="<< cmap.size()/3<<"\n";
  //shader->sendUniformi("dynamic_cmap",go->dynamic_cmap);
  shader->sendUniformi("reverse_cmap",go->reverse_cmap);
  shader->sendUniformf("powalpha",go->poweralpha);
  
  // send physical quantities stuffs
  if (hasPhysic&& go->render_mode==1) { // physic) {
    if (!go->dynamic_cmap) {
      // send absolute min and max phys of the object
      float moremax=(phys_select->getMax()-phys_select->getMin())/100.0;
      shader->sendUniformf("data_phys_min",log(phys_select->getMin()));
      shader->sendUniformf("data_phys_max",log(phys_select->getMax()+moremax));

      //std::cerr << "!go->dynamic_cmap : log(phys_select->getMin())="<<log(phys_select->getMin())<<"\n";
      //std::cerr << "!go->dynamic_cmap : log(phys_select->getMax())="<<log(phys_select->getMax())<<"\n";
    } else {
      // July 2014, really nasty bug,
      // because of low accuracy with float running on GLSL shaders
      // highest density points were not displayed !!! I fixed this issue
      // by adding 1% more values on max data
      float moremax=(po->getMaxPhys()-po->getMinPhys())/100.0;
      //std::cerr <<"moremax : " << moremax << "\n";
      shader->sendUniformf("data_phys_min",log(po->getMinPhys()));
      shader->sendUniformf("data_phys_max",log(po->getMaxPhys()+moremax));
    }
    
    //int imin=phys_itv[min_index].index;
    //int imax=phys_itv[max_index].index;

    //shader->sendUniformf("osel_phys_min",log(phys_select->data[imin]));
    //shader->sendUniformf("osel_phys_max",log(phys_select->data[imax]));
#if 0
    int imin=phys_itv[min_index].index;
    int imax=phys_itv[max_index].index;
    std::cerr << "imin = "<<imin<<" imax="<<imax<<"\n";
    std::cerr << "object :"<<phys_select->data[imin] << " --- " << phys_select->data[imax] << "\n";
    std::cerr << "data   :"<<phys_select->getMin() << " --- " << phys_select->getMax() << "\n";
#endif
    shader->sendUniformi("data_phys_valid",1);
  } else {
    shader->sendUniformf("data_phys_min",-1.);
    shader->sendUniformf("data_phys_max",-1.);
    //shader->sendUniformf("osel_phys_min",-1.);
    //shader->sendUniformf("osel_phys_max",-1.);
    shader->sendUniformi("data_phys_valid",0);
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
    const float vfactor = po->getVelSize();// / part_data->getMaxVelNorm(); // requested by Peter Teuben
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
// selectParticles();
// select particles and sort them according to their
// physical values
void GLObjectParticles::selectParticles()
{ 
  phys_itv.clear();   // clear physical value vector
  vindex_sel.clear(); // clear vindex vector     
  for (int i=0; i < po->npart; i+=po->step) {
    int index=po->index_tab[i];
    if (phys_select && phys_select->isValid()) {
      GLObjectIndexTab myphys;
      myphys.index   = index;
      myphys.value   = phys_select->data[index];
      myphys.i_point = i;
      phys_itv.push_back(myphys);
    }    
    // used if z depth test activated
    GLObjectIndexTab myz;
    myz.index = index;
    myz.i_point = i;
    vindex_sel.push_back(myz);
  }
  // sort by density
#if GLDRAWARRAYS
  sort(phys_itv.begin(),phys_itv.end(),GLObjectIndexTab::compareLow);
  //sort(rho.begin(),rho.end(),GLObjectIndexTab::compareHigh);
#endif
}

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
// buildIndexHisto()
// store first index particles belonging to each percent of physical array
// index_histo[percentile->0:99%]. It allows to find quickly the index of the
// first min and max particles.
// index_histo array stores the first index of the particle which belong to the
// percentile. Particles physical data must be sorted. Percentile
// vary from 0 to 99% of log of the physical data.
void GLObjectParticles::buildIndexHisto()
{
  // reset index_histo array
  index_histo.clear();
  for (int i=0; i <nhisto;i++) index_histo.push_back(-1);
  

  // compute first particle index in the percentage
  if (phys_select && phys_select->isValid()) {
  
    float log_part_r_max=0., log_part_r_min=0., diff_log_part=0.;
    float log_phys_max,log_phys_min, diff_log_phys;
    
    if (part_data->rho) {
      log_part_r_max = log(part_data->rho->getMax());
      log_part_r_min = log(part_data->rho->getMin());
      diff_log_part  = (nhisto-1.)/(log_part_r_max-log_part_r_min);
    }
    log_phys_max = log(phys_select->getMax());
    log_phys_min = log(phys_select->getMin());
    diff_log_phys= (nhisto-1.)/(log_phys_max-log_phys_min);
    
    int cpt=0;
    // find first index of particle in the percentage
    for (int i=0; i < po->npart; i+=po->step) {
      int index;//phys_itv[i].index;
#if 1
      if (po->rhoSorted() &&
          phys_select && phys_select->getType() != PhysicalData::rho && part_data->rho) {
        index = rho_itv[i].index; // it's temperature/pressure, we sort by density
        float rho_data=part_data->rho->data[index];
        if (rho_data!=0 && rho_data!=-1) {
          //std::cerr << "I="<<i<<" => "<<phys_select->data[index]<<"\n";
          int percen=(log(rho_data)-log_part_r_min)*diff_log_part;
          assert(percen<nhisto && percen>=0);
          if (index_histo[percen]==-1) { // no value yet
            index_histo[percen]=cpt; // store 
            
          }
          cpt++;
        }
      }
      else {        
        index = phys_itv[i].index; // we sort by physical value     
        float phys_data=phys_select->data[index];
        if (phys_data!=0 && phys_data!=-1) {
          //std::cerr << "I="<<i<<" => "<<phys_data<<"\n";
          int percen=(log(phys_data)-log_phys_min)*diff_log_phys;
          assert(percen<nhisto && percen>=0);
          if (index_histo[percen]==-1) { // no value yet
            index_histo[percen]=cpt; // store 
            
          }
          cpt++;
        }
      }
#endif
      
    }
    // fill empty index_histo
    int last=0;
    for (int i=0; i <nhisto; i++) {
      if (index_histo[i]==-1) index_histo[i]=last;
      else last=index_histo[i];
      //std::cerr << nhisto<< " Percentage["<<i<<"%]="<<index_histo[i]<<" quant="<<
      //  phys_select->data[phys_itv[index_histo[i]].index]<<"\n";
    }
    // if no physical quantity for the object
    if (index_histo[nhisto-1]==0) index_histo[nhisto-1] = nvert_pos;
  }
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
  texture = &(*gtv)[tex];
}
// ============================================================================
// setTexture()                                                                 
void GLObjectParticles::setTexture()
{
  if (!texture) {
    setTexture(0);
  }
}
// ============================================================================
// sortDyDensity                                                               
void GLObjectParticles::sortByDensity()
{
      // sort according to the density
  if (part_data->rho)  {
    sort(phys_itv.begin(),phys_itv.end(),GLObjectIndexTab::compareLow);
  }
  nind_sorted = 0;
    // creates vertex indices array for gLDrawElements
  if (! indexes_sorted || ((int) nind_sorted) < po->npart) {
    if (indexes_sorted) delete [] indexes_sorted;
    indexes_sorted = new GLuint[po->npart];
    nind_sorted = 0;
  }
  if ( phys_select && phys_select->isValid()) {
    PHYS_MIN = phys_select->getMin();
    PHYS_MAX = phys_select->getMax();
  }
  index_min = 0;
  index_max = po->npart;
  int cpt=0;
  for (int i=0; i < po->npart; i+=po->step) {
    int index;
    if (part_data->rho) index = phys_itv[i].i_point;
    
    if (part_data->rho) {
        index = phys_itv[i].i_point;
        int index2 = phys_itv[i].index;
    	if (part_data->rho->data[index2] >= PHYS_MIN && part_data->rho->data[index2] <= PHYS_MAX) {
        indexes_sorted[cpt++] = index;
	}
    } else {
      indexes_sorted[cpt++] = i;
    }
  }
  nind_sorted=cpt;
  //std::cerr << "index_min="<< index_min <<"   index_max=" << index_max << " index sorted="<< nind_sorted <<"\n";
      // bind VBO in order to use
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, vbo_index);
    // upload data to VBO
  //glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, (nind_sorted/2)*sizeof(GLuint), indexes_sorted, GL_STATIC_DRAW_ARB);
  glBufferDataARB(GL_ELEMENT_ARRAY_BUFFER_ARB, nind_sorted*sizeof(GLuint), indexes_sorted, GL_STATIC_DRAW_ARB);
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
    if (part_data->rho->data[index] >= PHYS_MIN && part_data->rho->data[index] <= PHYS_MAX) {
    float
        x=part_data->pos[index*3  ],
        y=part_data->pos[index*3+1],
        z=part_data->pos[index*3+2];
         // compute point coordinates according to model view via matrix
    float mz=  MM(2,0)*x + MM(2,1)*y + MM(2,2)*z + MM(2,3);//*w;
    vindex_sel[nind_sorted].value=mz;
    vindex_sel[nind_sorted].i_point=i;
    vindex_sel[nind_sorted].index = index;
    nind_sorted++;
   }
  }
  // sort according to the Z Depth
  sort(vindex_sel.begin(),vindex_sel.begin()+nind_sorted,GLObjectIndexTab::compareLow);
  //nind_sorted = 0;
  // creates vertex indices array for gLDrawElements
  if (! indexes_sorted || nind_sorted < vindex_sel.size()) {
    if (indexes_sorted) delete [] indexes_sorted;
    indexes_sorted = new GLuint[vindex_sel.size()];
    //nind_sorted = 0;
  }
  for (unsigned int i=0; i < nind_sorted; i++) {
    int index = vindex_sel[i].i_point;
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
