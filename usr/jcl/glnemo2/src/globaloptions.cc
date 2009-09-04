// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2009                                  
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
#include <QtGui>
#include "globaloptions.h"

namespace glnemo {
  float GlobalOptions::MAX_PARTICLES_SIZE   = 15.0;
  float GlobalOptions::MAX_TEXTURE_SIZE     = 1.0;
  float GlobalOptions::MAX_VEL_VECTOR_SIZE  = 2.0;
  bool  GlobalOptions::rho_exist         = false;
  bool  GlobalOptions::temperature_exist = false;
  bool  GlobalOptions::pressure_exist    = false;
// ============================================================================
// constructor                                                                 
GlobalOptions::GlobalOptions()
{
  // from PLAY options TAB
  auto_play_screenshot = false;
  auto_gl_screenshot   = false;
  auto_com             = false;
  play_fps             = 25;
  frame_index          = 0;
  base_frame_name      = "frame";
  base_frame_ext       = "png";
  auto_screen_size_index = 1;
  frame_width          = 1280;
  frame_height         = 720;
  // from OpenGL TAB
  show_part=true;
  show_points=false;
  show_vel=false;
  vel_vector_size=MAX_VEL_VECTOR_SIZE/2;
  psize=2.0;
  particles_alpha=255;
  blending=true;
  dbuffer=true;
  perspective=true;
  orthographic=false;
  use_shaders = true;
  render_mode = 2;
  init_glsl   = true;
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
  // from camera tab
  cam_mode = false;
  cam_display_ctrl = false;
  cam_display_spline = false;
  cam_play = false;
  // from HUD TAB
  show_osd=true;
  osd_title=osd_time=osd_zoom=osd_rot=true;
  osd_trans=osd_data_type=osd_nbody=osd_projection=true;
  background_color=QColor(Qt::black);
  osd_color=QColor(Qt::yellow);
  // from experimental TAB
  show_poly=false;
  texture_size=1.;
  texture_alpha_color=255;
  octree_enable=false;
  octree_display=false;
  octree_level=10;
  // vel 
  vel_req = false;
  new_frame = false;
  // memory
  duplicate_mem = true;
  octree_enable = false;
  // density
  phys_local = true;
  phys_min_glob = -1.;
  phys_max_glob = -1.;
  // colormap
  colormap=108;
  reverse_cmap = false;
  constant_cmap=true;
  powercolor = 0.8;
  poweralpha = 1.2;
  // Z sorting
  zsort = false;
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
  // from options PLAY TAB
  auto_play_screenshot = m.auto_play_screenshot;
  auto_gl_screenshot   = m.auto_gl_screenshot;
  auto_com             = m.auto_com;
  play_fps             = m.play_fps;
  frame_index          = m.frame_index;
  base_frame_name      = m.base_frame_name;
  base_frame_ext       = m.base_frame_ext;
  auto_screen_size_index = m.auto_screen_size_index;
  frame_width          = m.frame_width;
  frame_height         = m.frame_height;
  // from OpenGL TAB
  show_part=m.show_part;
  show_points=m.show_points;
  show_vel=m.show_vel;
  vel_vector_size=m.vel_vector_size;
  psize=m.psize;
  particles_alpha=m.particles_alpha;
  blending=m.blending;
  dbuffer=m.dbuffer;
  perspective=m.perspective;
  orthographic=m.orthographic;
  use_shaders =m.use_shaders;
  render_mode =m.render_mode;
  init_glsl   =m.init_glsl;
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
  show_osd=m.show_osd;
  osd_title=m.osd_title;
  osd_time=m.osd_time;
  osd_zoom=m.osd_zoom;
  osd_rot=m.osd_rot;
  osd_trans=m.osd_trans;
  osd_data_type=m.osd_data_type;
  osd_nbody=m.osd_nbody;
  osd_projection=m.osd_projection;
  background_color=m.background_color;
  osd_color=m.osd_color;
  // from camera tab
  cam_display_ctrl   = m.cam_display_ctrl;
  cam_mode           = m.cam_mode;
  cam_display_spline = m.cam_display_spline;
  cam_play           = m.cam_play;
  // from experimental TAB
  show_poly=m.show_poly;
  texture_size=m.texture_size;
  texture_alpha_color=m.texture_alpha_color;
  octree_enable=m.octree_enable;
  octree_display=m.octree_display;
  octree_level=m.octree_level;
  // vel 
  vel_req = m.vel_req;
  new_frame = m.new_frame;
  // memory
  duplicate_mem = m.duplicate_mem;
  // physical value
  phys_local = m.phys_local;
  phys_max_glob = m.phys_max_glob;
  phys_min_glob = m.phys_min_glob;
  // colormap
  colormap = m.colormap;
  reverse_cmap = m.reverse_cmap;
  constant_cmap=m.constant_cmap;
  R = m.R;
  G = m.G;
  B = m.B;
  powercolor = m.powercolor;
  poweralpha = m.poweralpha;
  // z sorting
  zsort = m.zsort;
  return *this;
}
// ============================================================================
}
