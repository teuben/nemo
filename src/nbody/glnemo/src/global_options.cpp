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
#include "global_options.h"
#include <iostream>
#include <qnamespace.h>

// ============================================================================
// constructor                                                                 
GlobalOptions::GlobalOptions()
{
  // SET default parameters
  MAX_PARTICLES_SIZE=5.0;
  MAX_TEXTURE_SIZE=1.0;
  // from OpenGL TAB
  show_part=true;
  psize=1.5;
  particles_alpha=255;
  blending=true;
  dbuffer=true;
  perspective=true;
  orthographic=false;
  // from Scene Orientation TAB
  zoom=-1.0;
  zoomo=1.0;
  xrot=yrot=xrot=0.0;
  xtrans=ytrans=ztrans=0.0;
  // from Grids TAB
  show_grid=true;
  mesh_length=1.0;
  nb_meshs=28;
  xy_grid=true;
  yz_grid=xz_grid=false;
  col_x_grid = QColor(136,141,102);
  col_y_grid = QColor(136,141,102);
  col_z_grid = QColor(136,141,102);
  // from HUD TAB
  hud=true;
  hud_title=hud_time=hud_zoom=hud_rot=true;
  hud_trans=hud_data_type=hud_nbody=hud_projection=true;
  background_color=QColor(Qt::black);
  hud_color=QColor(Qt::yellow);
  // from experimental TAB
  show_poly=false;  
  texture_size=0.52;
  texture_alpha_color=125;
  octree_enable=false;
  octree_display=false;
  octree_level=10;
}
// ============================================================================
// destructor                                                                  
GlobalOptions::~GlobalOptions()
{
}
// ============================================================================
// copyTransform()                                                             
// copy transformation data from another object                                
void GlobalOptions::copyTransform(const GlobalOptions &m)
{
  // rot
  xrot    = m.xrot;
  yrot    = m.yrot;
  zrot    = m.zrot;
  // trans
  xtrans = m.xtrans;
  ytrans = m.ytrans;
  ztrans = m.ztrans;
  // zoom
  zoom   = m.zoom;
  zoomo  = m.zoomo;
}
// ============================================================================



