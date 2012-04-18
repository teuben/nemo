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
#include "glselection.h"
#include "particlesobject.h"
#include "frustumculling.h"
#include "tools3d.h"
#include "vec3d.h"

namespace glnemo {
#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mModel[col*4+row]
#define MM2(row,col)  mModel2[col*4+row]
// ============================================================================
// constructor                                                                 
GLSelection::GLSelection()
{
  reset();
  zoom        = true;
  anim_zoom   = true;
  total_frame = 25;
  anim_timer  = new QTimer(this);
  connect(anim_timer, SIGNAL(timeout()), this, SLOT(playZoomAnim()));
}

// ============================================================================
// destructor                                                                  
GLSelection::~GLSelection()
{
}
// ============================================================================
// void update
void GLSelection::reset()
{
  enable=false;
  x0    =-1;
}
// ============================================================================
// void update
void GLSelection::update(const GLObjectParticlesVector * _gpv,
                         GlobalOptions   * _go, QMutex * _mutex)
{
  // update variables
  gpv           = _gpv;
  store_options = _go;
  mutex_data    = _mutex;
}
// ============================================================================
//
void GLSelection::getMouse(QMouseEvent * e)
{
  enable=true;
  if (x0 == -1) {
    x0 = e->x();
    y0 = e->y();
  }
  x1 = e->x();
  y1 = e->y();
}
// ============================================================================
//
void GLSelection::display(const int width, const int height)
{
  if (enable) {
    glDisable( GL_DEPTH_TEST );
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    //glOrtho(0.,width,0.,height,-1,1);
    gluOrtho2D(0.,width,0.,height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    //glBlendFunc( GL_SRC_ALPHA, GL_ONE ); // original
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// No Alpha bending accumulation
    glEnable(GL_BLEND);
    // draw blended box
    glColor4f( 1.0f, 1.0f, 1.f,0.4f );
    glBegin(GL_QUADS);
      glVertex2f(x0,height-y0);
      glVertex2f(x1,height-y0);
      glVertex2f(x1,height-y1);
      glVertex2f(x0,height-y1);
    glEnd();
    // draw surounded lines
    glColor4f( 1.0f, 0.f, 0.f,1.f );
    glBegin(GL_LINES);
      glVertex2f(x0,0);
      glVertex2f(x0,height);
      glVertex2f(x1,0);
      glVertex2f(x1,height);
      glVertex2f(0,height-y0);
      glVertex2f(width,height-y0);
      glVertex2f(0,height-y1);
      glVertex2f(width,height-y1);
    glEnd();
    glDisable(GL_BLEND);
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    glEnable( GL_DEPTH_TEST );
  }
}
// ============================================================================
// zoomOnArea                                                                  
// according to the user selection, the selected area will be centered and     
// zoomed in                                                                   
void GLSelection::zoomOnArea(const int nobj, double mProj[16],double mModel[16],
                             const int viewport[4])
{
  if (nobj&&enable) {
    // reordering square selection on X
    if (x1<x0) {
      float xx=x0;
      x0=x1;x1=xx;
    }
    // reordering square selection on Y
    if (y1<y0) {
      float yy=y0;
      y0=y1;y1=yy;
    }
    
    FrustumCulling frustum;
    frustum.getFC(mModel,mProj); // compute frustum  
    Tools3D t3d(mModel,mProj);   // 3D stuffs        

    const ParticlesData * part_data;
    double com[3] = { 0.,0.,0.};  // Center Of Mass   
    int in_area=0;
    list.clear();
    // loop on all the objects
    // to find out all particles in the
    // selected area
    for (int i=0; i<nobj; i++) {
      const ParticlesObject * po = (*gpv)[i].getPartObj();        // object
      float DMIN = po->getMinPhys();
      float DMAX = po->getMaxPhys();

      if (po->isVisible()) {                                   // is visible  
        part_data = (*gpv)[i].getPartData();// get its Data
        // get physical value data array
        PhysicalData * phys_select = part_data->getPhysData();
        // loop on all the particles of the object
        for (int j  = 0; j  <  po->npart; j ++) {
          int jndex= po->index_tab[j];
          float
          x=part_data->pos[jndex*3  ]+store_options->xtrans,
          y=part_data->pos[jndex*3+1]+store_options->ytrans,
          z=part_data->pos[jndex*3+2]+store_options->ztrans;

          // compute the projection on the selected point
          Vec3D v3d=t3d.projPoint(x,y,z);
          // compute screen coordinates
          float winx=viewport[0] + (1 + v3d.x) * viewport[2] / 2;
          float winy=viewport[1] + (1 + v3d.y) * viewport[3] / 2;
          winy = viewport[3]-winy;
          bool indensity=true;
          if (phys_select && phys_select->isValid()) { 
            if (phys_select->data[jndex] >= DMIN && 
                phys_select->data[jndex] <= DMAX) {
                indensity=true;
            } else {
                indensity=false;
            }
          }
          //std::cerr << winx << " " << x0 << " " << x1 << " " << winy << " " << y0 << " " << y1 << " " <<frustum.isPointInside(x,y,z) << "\n"; 
          // is particle visible ?
          if (indensity && winx >= x0 && winx <= x1 && winy >= y0 && winy <= y1 && // in selected area
              frustum.isPointInside(x,y,z)) {                         // in the frustum  
            if (part_data->id.size()>0) {
              //jndex = part_data->id.at(jndex);
              //jndex = part_data->id[jndex];
            }
            list.push_back(jndex); // save particle index
            in_area++;             // one more particle
            com[0]+=x;             // x COM
            com[1]+=y;             // y COM
            com[2]+=z;             // z COM
          }
        }
      }
    }
    reset();
    emit updatePareticlesSelected(in_area); // update Form Option
    if (in_area) { // particles exist
        // normalizing COM
	com[0] /= (float) in_area;
	com[1] /= (float) in_area;
	com[2] /= (float) in_area;        
	
        //std::cerr << "in_area = " << in_area << "  list ="<<list.size()<<"\n";
	//std::cerr << "center of mass:" << com[0] << " " << com[1] << " " << com[2] <<"\n";
        //std::cerr << "zoom in=" << store_options->zoom << "\n";
        // save information
        float zoom1 =store_options->zoom; // zoom value
        float zoomo1=store_options->zoomo; // zoom value
        float ortho1=store_options->ortho_range; // zoom value
        trans_in.set(store_options->xtrans,store_options->ytrans,store_options->ztrans);
        if (zoom) {// best ZOOM on particles inside selected area 
           // centering on COM
           store_options->xtrans -= com[0];
           store_options->ytrans -= com[1];
           store_options->ztrans -= com[2];
           trans_out.set(store_options->xtrans,store_options->ytrans,store_options->ztrans);

           // in following function we compute the best
           // zoom for perspective and orthographic projection
           // BUT for orthographic, best zoom is set to ortho_range
           Tools3D::bestZoomFromList(mProj,mModel,viewport,&list, part_data, store_options);
           if (anim_zoom) {
            float zoom2 =store_options->zoom; // new zoom value            
            // perspective zoom offset
            zoom_dynamic =(zoom2-zoom1)/float(total_frame); // animation zoom value offset
            // orthoraphic zoom offset
            zoomo_dynamic=(store_options->ortho_range-ortho1*zoomo1)/float(total_frame); 
            //std::cerr << "ortho0 =" << store_options->ortho_range << " ortho1="<< ortho1 << "\n";
            //std::cerr << "zoomo_dynmaic = "<< zoomo_dynamic << " zoomo1="<<zoomo1<<"\n";
            // set initial Center
            store_options->xtrans=trans_in[0];
            store_options->ytrans=trans_in[1];
            store_options->ztrans=trans_in[2];
  
            comvec = trans_out-trans_in; // vector director to COM
            Vec3D v;
            v = comvec + comvec;
            store_options->zoom  = zoom1; // set initial zoom value
            store_options->zoomo = 1.; // we want to keep zoomo = 1
            store_options->ortho_range = ortho1*zoomo1; // initial range with zoomo = 1
            frame_counter = 0;
            mutex_data->lock();          // keep priority on data
            anim_timer->start(20);       // start zoom animation 
          }
          else emit updateZoom();          // update GL
        }
        //std::cerr << "zoom out=" << store_options->zoom << "\n";
    }
  }
  else reset();
}

void GLSelection::playZoomAnim()
{
  frame_counter++;                                       // one more frame         
  if (frame_counter<=total_frame) {                      // frame exist            
    float off = float(frame_counter)/float(total_frame); // new displacement offset
    store_options->zoom  += zoom_dynamic;                 // new zoom               
    //store_options->zoomo += zoomo_dynamic;                 // new zoom   
    store_options->ortho_range += zoomo_dynamic;
    store_options->xtrans=trans_in[0]+(off*comvec[0]);   // new x center           
    store_options->ytrans=trans_in[1]+(off*comvec[1]);   // new y center           
    store_options->ztrans=trans_in[2]+(off*comvec[2]);   // new z center           
    emit updateZoom();                                   // update GL
  }
  else {
    
    anim_timer->stop();
    mutex_data->unlock(); // release the data
  }
}
} // glnemo

