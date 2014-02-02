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
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@lam.fr>
 */
#ifndef GLNEMOGLOBALOPTIONS_H
#define GLNEMOGLOBALOPTIONS_H

#include <QColor>
#include <vector>
namespace glnemo {
 
class GlobalOptions{
public:
    GlobalOptions();
    ~GlobalOptions();
    const GlobalOptions& operator=(const GlobalOptions& m);
	// Network stuff
	std::string  network_host;
	int  network_port;
     bool list_type;
    // from play option TAB
    bool auto_play_screenshot;
    bool auto_gl_screenshot;
    bool auto_com;
    bool play_forward;
    int jump_frame;
    int play_fps;
    int frame_index;
    int frame_height;
    int frame_width;
    QString base_frame_name;
    QString base_frame_ext;
    int auto_screen_size_index;
    // from OpenGL TAB
    bool show_part;
    bool show_points;
    bool show_vel;
    float psize;
    float vel_vector_size;
    int particles_alpha;
    bool blending;
    bool dbuffer;
    bool perspective;
    bool orthographic;
    bool use_shaders;
    int render_mode;
    bool init_glsl;
    // from camera tab
    bool cam_mode;
    bool cam_display_ctrl;
    bool cam_display_spline;
    bool cam_play;
    // from Scene Orientation TAB
    float zoom,zoomo;
    float xrot,yrot,zrot;
    float urot,vrot,wrot;
    float ixrot,iyrot,izrot;
    float iurot,ivrot,iwrot;
    bool xbrot,ybrot,zbrot,ubrot,vbrot,wbrot;
    float xtrans,ytrans,ztrans;
    // from Grids TAB
    bool show_grid;
    float mesh_length;
    int nb_meshs;
    bool xy_grid, yz_grid, xz_grid, show_cube;
    QColor col_x_grid, col_y_grid, col_z_grid, col_cube;    
    // from OSD TAB
    bool show_osd;
    bool osd_title, osd_time, osd_zoom, osd_rot,
    osd_trans, osd_data_type, osd_nbody, osd_projection;
    QColor background_color, osd_color;
    QString osd_font_name, osd_title_name;
    float osd_font_size;
    // from experimental TAB
    bool show_poly;
    float texture_size, texture_alpha;
    bool octree_enable;
    bool octree_display;
    int octree_level;
    // from rotation/axis tab
    // axes
    bool  axes_enable;
    int   axes_loc; // 0 bottom, 1 center
    float axes_psize; // percentage of windows width size
    bool rotate_screen; // 0 world, 1 screen
    
    // Colorbar
    bool   gcb_enable, gcb_logmode;
    int    gcb_min, gcb_max; // threshold
    int    gcb_ndigits;      // #digits
    int    gcb_orientation;  // North,Est,South,West
    float  gcb_pwidth, gcb_pheight; // factor size
    float  gcb_font_size;
    int    gcb_offset;       // colorbar offset in pixels
    QString gcb_font_name;
    QColor  gcb_color;
    // network;
    int port;
    // const
    static float MAX_PARTICLES_SIZE;
    static float MAX_TEXTURE_SIZE;
    static float MAX_VEL_VECTOR_SIZE;
    static QString RESPATH;
    // velocity
    bool vel_req;
    // cod
    bool cod;
    // method
    void copyTransform(const GlobalOptions &m);
    // frame
    bool new_frame;
    // memory
    bool duplicate_mem;
    // density
    bool  phys_local;
    float phys_max_glob;
    float phys_min_glob;
    static bool rho_exist;
    // temperature
    static bool temperature_exist;
    // pressure
    static bool pressure_exist;
    // colormap
    int colormap;
    bool reverse_cmap; // reverse colormap
    std::vector <float> *R,*G,*B;
    float powercolor, poweralpha;
    bool dynamic_cmap;
    // z sorting
    bool zsort;
    bool enable_gui;
    // texture size
    bool auto_texture_size;
    // ramses stuffs
    float xmin,xmax,ymin,ymax,zmin,zmax;
    int lmax,lmin;
    float scale;
    // select
    std::string select_time;
    std::string select_part;
    bool auto_render;
    // orthographic range
    float ortho_range;
    // opaque disc
    bool od_enable;
    bool od_display;
    float od_radius;
};

}

#endif
