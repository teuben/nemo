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
#include "tools3d.h"

namespace glnemo {

Tools3D::Tools3D(double * _m, double * _p)
{
  for (int i=0; i<16; i++) {
    mModel[i] = _m[i];
    mProj[i]  = _p[i];
  }
}


Tools3D::~Tools3D()
{
}
// ============================================================================
// bestZoomFromObject                                                          
// fit all the particles on the screen from perspective view                   
// For each point v(x,y,z), we have to compute its screen coordinates win(x,y) 
// according to its Projection and  ModelView matrix                           
// vp = MProj x MModelView x v                                                 
// winx = viewport[0] + (1 + vpx) * viewport[2] / 2                            
// winy = viewport[1] + (1 + vpy) * viewport[3] / 2                            
//                                                                             
// We have then to resolve the previous equation to figure out the best value  
// for zoom                                                                    
void Tools3D::bestZoomFromObject(double * mProj,double * mModel,
                                  const int * viewport, const ParticlesObjectVector * pov,
                                  const ParticlesData * part_data, GlobalOptions * store_options)
{
  //glGetIntegerv(GL_VIEWPORT,viewport);
  if (pov && pov->size() ) {
    // force ZOOM to fit all particles
    // Zoom is located in ModelView matrix at coordinates MM(2,3)
    double best_zoom;
    
    MM(2,3)  =    -20000000.0;
    best_zoom=     std::numeric_limits<double>::max();//10000;
    std::cerr << "MM = " << MM(2,3)  << " best zoom ="<<best_zoom<<"\n";
    float mid_screenx = (viewport[2]-viewport[0])/2.;
    float mid_screeny = (viewport[3]-viewport[1])/2.;
    float coo[3];
    double absxmax=0.;//fabs(std::numeric_limits<double>::min());
    double absymax=0.;//fabs(std::numeric_limits<double>::min());
    std::cerr << "xy max "<<absxmax << " " << absymax << "\n";
    // loop on all the objects
    for (int i=0; i<(int)pov->size(); i++) {
      //const ParticlesObject * po = gpv[i].getPartObj();        // object
      const ParticlesObject * po = &(*pov)[i];
      if (po->isVisible()) {                                   // is visible  
        //const ParticlesData * part_data = gpv[i].getPartData();// get its Data
        // loop on all the particles of the object
        for (int j  = 0; j  <  po->npart; j ++) {
          int jndex= po->index_tab[j];
          float
          x=part_data->pos[jndex*3  ]+store_options->xtrans,
          y=part_data->pos[jndex*3+1]+store_options->ytrans,
          z=part_data->pos[jndex*3+2]+store_options->ztrans,
          w=1.;
          // do the product Mmodel X point = mxyzw
          float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3)*w;
          float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3)*w;
          float Mmz= MM(2,0)*x + MM(2,1)*y + MM(2,2)*z;
          float mz = Mmz + MM(2,3)*w;
          float mw = MM(3,0)*x + MM(3,1)*y + MM(3,2)*z + MM(3,3)*w;
          // do the product Mproj X mxyzw  = pxyzw
          float Ppx= MP(0,0)*mx + MP(0,1)*my + MP(0,3)*mw;
          float px = Ppx + MP(0,2)*mz;
          float Ppy= MP(1,0)*mx + MP(1,1)*my + MP(1,3)*mw;
          float py = Ppy + MP(1,2)*mz;
          float Ppz= MP(2,0)*mx + MP(2,1)*my + MP(2,3)*mw;
          float pz = Ppz + MP(2,2)*mz;
          float Ppw= MP(3,0)*mx + MP(3,1)*my + MP(3,3)*mw;
          float pw = Ppw + MP(3,2)*mz;
          // normalyze
          px /= pw;
          py /= pw;
          pz /= pw;
          // compute orthographic best zoom
          if (fabs(x)>fabs(absxmax)) absxmax =x;
          if (fabs(y)>fabs(absymax)) absymax =y; 
          
          //std::cerr << px << " " << py << " " << pz << "\n";
          // compute screen coordinates
          float winx=viewport[0] + (1 + px) * viewport[2] / 2;
          float winy=viewport[1] + (1 + py) * viewport[3] / 2;
          // paint particles
          //paint.setPen(red);
          //paint.drawPoint((int) (winx), (int) (winy));
          //paint.drawPoint(90,10);
          // check farest particle
          bool guess_out_zoomx=false;
          bool guess_out_zoomy=false;
          float screen_coo;
  
          // proceed from left to mid-side of the screen
          if (winx >= viewport[0] && winx <= mid_screenx) {
            screen_coo=viewport[0];
            guess_out_zoomx=true;
          }
          else {
            // proceed from right to mid-side of the screen
            if (winx>= mid_screenx && winx <= viewport[2]) {
              screen_coo=viewport[2];
              guess_out_zoomx=true;
            }
          }
          if (guess_out_zoomx) {
            float A=(2.*(screen_coo-viewport[0])/viewport[2])-1.;
            float new_zoomx=-Mmz + (-Ppx+A*Ppw)/(MP(0,2)-A*MP(3,2));
            if (new_zoomx < best_zoom) {
              best_zoom = new_zoomx;
              coo[0] = x; coo[1] = y; coo[2] = z;
            }
          }
          // proceed from bottom to mid-side of the screen
          if (winy <= viewport[3] && winy >= mid_screeny) {
            screen_coo=viewport[3];
            guess_out_zoomy=true;
          }
          else {
            // proceed from top to mid-side of the screen
            if (winy >= viewport[1] && winy <= mid_screeny) {
              screen_coo=viewport[1];
              guess_out_zoomy=true;
            }
          }
          if (guess_out_zoomy) {
            float A=(2*(screen_coo-viewport[1])/viewport[3])-1;
            float new_zoomy=-Mmz + (-Ppy+A*Ppw)/(MP(1,2)-A*MP(3,2));
            if (new_zoomy < best_zoom) {
              best_zoom = new_zoomy;
              coo[0] = x; coo[1] = y; coo[2] = z;
            }
          }
        }
      }
      else { // object not visible
      }
    }
    store_options->ortho_range=std::max(fabs(absxmax),fabs(absymax));
    std::cerr << "ortho_range = " << store_options->ortho_range << "\n";
    //setZoom( best_zoom);
    store_options->zoom  = best_zoom;
    std::cerr << "[" << best_zoom << "] Cordinates for best zoom = "
                      << coo[0] <<" " << coo[1] << " " << coo[2] << "\n";
  }
}
// ============================================================================
// bestZoomFromList                                                            
// fit all the particles on the screen from perspective view                   
// For each point v(x,y,z), we have to compute its screen coordinates win(x,y) 
// according to its Projection and  ModelView matrix                           
// vp = MProj x MModelView x v                                                 
// winx = viewport[0] + (1 + vpx) * viewport[2] / 2                            
// winy = viewport[1] + (1 + vpy) * viewport[3] / 2                            
//                                                                             
// We have then to resolve the previous equation to figure out the best value  
// for zoom                                                                    
void Tools3D::bestZoomFromList(double * mProj,double * mModel,
                              const int * viewport, const std::vector <int> * list,
                              const ParticlesData * part_data, GlobalOptions * store_options)
{

  if (list->size() ) {
    // force ZOOM to fit all particles
    // Zoom is located in ModelView matrix at coordinates MM(2,3)
    float best_zoom;
    
    MM(2,3) = -20000000.0;
    best_zoom=std::numeric_limits<double>::max();//10000;
    
    float mid_screenx = (viewport[2]-viewport[0])/2.;
    float mid_screeny = (viewport[3]-viewport[1])/2.;
    float coo[3];
    double absxmax=0.;//fabs(std::numeric_limits<double>::min());
    double absymax=0.;//fabs(std::numeric_limits<double>::min());
    // loop on all the objects
    for (int i=0; i<(int)list->size(); i++) {

          int jndex= (*list)[i];
          float
          x=part_data->pos[jndex*3  ]+store_options->xtrans,
          y=part_data->pos[jndex*3+1]+store_options->ytrans,
          z=part_data->pos[jndex*3+2]+store_options->ztrans,
          w=1.;
          // ---- compute the projection on the selected point --------
          // do the product Mmodel X point = mxyzw
          float mx = MM(0,0)*x + MM(0,1)*y + MM(0,2)*z + MM(0,3)*w;
          float my = MM(1,0)*x + MM(1,1)*y + MM(1,2)*z + MM(1,3)*w;
          float Mmz= MM(2,0)*x + MM(2,1)*y + MM(2,2)*z;
          float mz = Mmz + MM(2,3)*w;
          float mw = MM(3,0)*x + MM(3,1)*y + MM(3,2)*z + MM(3,3)*w;
          // do the product Mproj X mxyzw  = pxyzw
          float Ppx= MP(0,0)*mx + MP(0,1)*my + MP(0,3)*mw;
          float px = Ppx + MP(0,2)*mz;
          float Ppy= MP(1,0)*mx + MP(1,1)*my + MP(1,3)*mw;
          float py = Ppy + MP(1,2)*mz;
          float Ppz= MP(2,0)*mx + MP(2,1)*my + MP(2,3)*mw;
          float pz = Ppz + MP(2,2)*mz;
          float Ppw= MP(3,0)*mx + MP(3,1)*my + MP(3,3)*mw;
          float pw = Ppw + MP(3,2)*mz;
          // normalyze
          px /= pw;
          py /= pw;
          pz /= pw;
          
          // compute orthographic best zoom
          if (fabs(x)>fabs(absxmax)) absxmax =x;
          if (fabs(y)>fabs(absymax)) absymax =y; 
          
          // compute screen coordinates
          float winx=viewport[0] + (1 + px) * viewport[2] / 2.;
          float winy=viewport[1] + (1 + py) * viewport[3] / 2.;

          // check farest particle
          bool guess_out_zoomx=false;
          bool guess_out_zoomy=false;
          float screen_coo;

          // proceed from left to mid-side of the screen
          if (winx >= viewport[0] && winx <= mid_screenx) {
            screen_coo=viewport[0];
            guess_out_zoomx=true;
          }
          else {
            // proceed from right to mid-side of the screen
            if (winx>= mid_screenx && winx <= viewport[2]) {
              screen_coo=viewport[2];
              guess_out_zoomx=true;
            }
          }
          if (guess_out_zoomx) {
            float A=(2.*(screen_coo-viewport[0])/viewport[2])-1.;
            float new_zoomx=-Mmz + (-Ppx+A*Ppw)/(MP(0,2)-A*MP(3,2));
            if (new_zoomx < best_zoom) {
              best_zoom = new_zoomx;
              coo[0] = x; coo[1] = y; coo[2] = z;
            }
          }
          // proceed from bottom to mid-side of the screen
          if (winy <= viewport[3] && winy >= mid_screeny) {
            screen_coo=viewport[3];
            guess_out_zoomy=true;
          }
          else {
            // proceed from top to mid-side of the screen
            if (winy >= viewport[1] && winy <= mid_screeny) {
              screen_coo=viewport[1];
              guess_out_zoomy=true;
            }
          }
          if (guess_out_zoomy) {
            float A=(2*(screen_coo-viewport[1])/viewport[3])-1;
            float new_zoomy=-Mmz + (-Ppy+A*Ppw)/(MP(1,2)-A*MP(3,2));
            if (new_zoomy < best_zoom) {
              best_zoom = new_zoomy;
              coo[0] = x; coo[1] = y; coo[2] = z;
            }
          }
    }
    store_options->ortho_range=std::max(fabs(absxmax),fabs(absymax));
    store_options->zoom=best_zoom;
  }
}
}
