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
#ifndef GLOBAL_OPTIONS_H
#define GLOBAL_OPTIONS_H

#include <qgl.h>
#include <qcolor.h>
class GlobalOptions{
public:
    GlobalOptions();

    ~GlobalOptions();

    // from OpenGL TAB
    bool show_part;
    float psize;
    int particles_alpha;
    bool blending;
    bool dbuffer;
    bool perspective;
    bool orthographic;
    // from Scene Orientation TAB
    float zoom,zoomo;
    float xrot,yrot,zrot;
    float xtrans,ytrans,ztrans;
    // from Grids TAB
    bool show_grid;
    float mesh_length;
    int nb_meshs;
    bool xy_grid, yz_grid, xz_grid;
    QColor col_x_grid, col_y_grid, col_z_grid;
    // from HUD TAB
    bool hud;
    bool hud_title, hud_time, hud_zoom, hud_rot,
      hud_trans, hud_data_type, hud_nbody, hud_projection;
    QColor background_color, hud_color;
    // from experimental TAB
    bool show_poly;
    float texture_size;
    int texture_alpha_color;
    bool octree_enable;
    bool octree_display;
    int octree_level;
    // network;
    int port;
    // const
    float MAX_PARTICLES_SIZE;
    float MAX_TEXTURE_SIZE;
    
    // method
    void copyTransform(const GlobalOptions &m);
};
#endif
// ============================================================================
