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
#include "global_options.h"
#include <iostream>
#include <qnamespace.h>

// ============================================================================
// constructor                                                                 
GlobalOptions::GlobalOptions()
{
  // SET default parameters
  MAX_PARTICLES_SIZE   = 5.0;
  MAX_TEXTURE_SIZE     = 1.0;
  MAX_VEL_VECTOR_SIZE  = 2.0;
  // from OpenGL TAB
  show_part=true;
  show_vel=false;
  vel_vector_size=MAX_VEL_VECTOR_SIZE/2;
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
  yz_grid=xz_grid= show_cube =false;
  col_x_grid = QColor(136,141,102);
  col_y_grid = QColor(136,141,102);
  col_z_grid = QColor(136,141,102);
  col_cube   = QColor(136,141,102);
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
  // vel 
  vel_req = false;
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
// copyTransform()                                                             
// copy transformation data from another object                                
const GlobalOptions& GlobalOptions::operator=(const GlobalOptions &m)
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
  
    // SET default parameters
  MAX_PARTICLES_SIZE   = m.MAX_PARTICLES_SIZE;
  MAX_TEXTURE_SIZE     = m.MAX_TEXTURE_SIZE;
  MAX_VEL_VECTOR_SIZE  = m.MAX_VEL_VECTOR_SIZE;
  // from OpenGL TAB
  show_part=m.show_part;
  show_vel=m.show_vel;
  vel_vector_size=m.vel_vector_size;
  psize=m.psize;
  particles_alpha=m.particles_alpha;
  blending=m.blending;
  dbuffer=m.dbuffer;
  perspective=m.perspective;
  orthographic=m.orthographic;
  // from Scene Orientation TAB
  
  // from Grids TAB
  show_grid=m.show_grid;
  mesh_length=m.mesh_length;
  nb_meshs=m.nb_meshs;
  xy_grid=m.xy_grid;
  yz_grid=m.yz_grid;
  xz_grid= m.xy_grid;
  show_cube =m.show_cube;
  col_x_grid = m.col_x_grid;
  col_y_grid = m.col_y_grid;
  col_z_grid = m.col_z_grid;
  col_cube   = m.col_cube;
  // from HUD TAB
  hud=m.hud;
  hud_title=m.hud_title;
  hud_time=m.hud_time;
  hud_zoom=m.hud_zoom;
  hud_rot=m.hud_rot;
  hud_trans=m.hud_trans;
  hud_data_type=m.hud_data_type;
  hud_nbody=m.hud_nbody;
  hud_projection=m.hud_projection;
  background_color=m.background_color;
  hud_color=m.hud_color;
  // from experimental TAB
  show_poly=m.show_poly;
  texture_size=m.texture_size;
  texture_alpha_color=m.texture_alpha_color;
  octree_enable=m.octree_enable;
  octree_display=m.octree_display;
  octree_level=m.octree_level;
  // vel 
  vel_req = m.vel_req;
  return *this;
}
// ============================================================================



