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
/**
	@author Jean-Charles Lambert <Jean-Charles.Lambert@oamp.fr>
 */
#ifndef GLNEMOMAINWINDOW_H
#define GLNEMOMAINWINDOW_H

// QT's includes
#include <QMainWindow>
#include <QTimer>
#include <QMutex>

// project's includes

#include "globaloptions.h"
#include "particlesobject.h"
#include "loadingthread.h"
#include "formobjectcontrol.h"
#include "formabout.h"
#include "formhelp.h"
#include "formscreenshot.h"
#include "formselectpart.h"
#include "formoptions.h"
#include "formconnect.h"
#include "componentrange.h"
#include "userselection.h"
#include "colormap.h"
#include "camera.h"

// prototypes
class QAction;
class QMenu;


namespace glnemo {
  // prototypes of class in glnemo namespace
  class GLWindow;
  class PluginsManage;
  class SnapshotInterface;
  class UserSelection;
class MainWindow : public QMainWindow {
  Q_OBJECT

  public:
    MainWindow(std::string);
    void start(std::string);
    ~MainWindow();
  private slots:
    // action
    void actionFullScreen();
    void actionBestZoom();
    void actionCenterToCom(const bool ugl=true);
    void actionRenderMode();
    void actionScreenshot();
    void actionZSorting();
    void actionPlay();
    void actionGrid();
    void actionFormObjectControl();
    void actionFormOptions();
    void actionFormOptionsShowMovie();
    void actionMenuFileOpen();
	void actionMenuFileConnect();
    bool actionMenuFileConnect2(std::string ip, int port, bool velocities, bool dencities, bool interactivSelect);
    void actionReset();
    void actionReload();
    void actionPrint();
    void actionAutoRotate(const int);
    void actionRotateX (const bool b=false);
    void actionRotateRX(const bool b=false);
    void actionRotateY (const bool b=false);
    void actionRotateRY(const bool b=false);
    void actionRotateZ (const bool b=false);
    void actionRotateRZ(const bool b=false);
    void actionRotateU (const bool b=false);
    void actionRotateRU(const bool b=false);
    void actionRotateV (const bool b=false);
    void actionRotateRV(const bool b=false);
    void actionRotateW (const bool b=false);
    void actionRotateRW(const bool b=false);    
    void actionTranslateX();
    void actionTranslateY();
    void actionTranslateZ();
    void actionToggleOsd();
    void actionToggleRotationScreen();
    void actionToggleCamera();
    void actionToggleProjection();
    void actionEmpty() {;}
    void actionQuit();
    void actionAutoScreenshot();
    void actionGLAutoScreenshot();
    void playEvent();
    void pressedKeyMouse(const bool, const bool);
    void uploadNewFrame();
    void takeScreenshot(const int, const int, std::string name="");
    void selectPart(const std::string, const bool);
    void startBench(const bool);
    void updateBenchFrame();
    void startAutoScreenshot();
    void saveIndexList();
    void createObjFromIndexList();
    void updateIpvs(const int ipvs=-1) {
       if (current_data) {
         current_data->part_data->setIpvs(ipvs);
       }
    }
    // OSD
    void updateOsd(bool ugl=false);
    //
    void toggleRotateScreen();
  private:

    // ------------------------------
    // Graphic User Interface        
    //-------------------------------
    void createForms();
    void createActions();
    void createMenus();
    void createToolBars();
    void createStatusBar();
    void createDockWindows();
    // -------------------------------
    // method                         
    // -------------------------------
    void interactiveSelect(std::string s="", const bool first_snaphot=false);
    void loadNewData(const std::string, const std::string ,
                            const bool , const bool, const bool first=false);
    void killPlayingEvent();
    
    // Menus
    QMenu *file_menu;
    QMenu *help_menu;
    // ToolsBars
    QToolBar * icons_tool_bar;
    // File menu actions
    QAction * open_file_action;
    QAction * connect_file_action;
    QAction * print_file_action;
    QAction * quit_file_action;
    // Help menu actions
    QAction * doc_action;
    QAction * about_action;
    QAction * about_qt;
    // Toolbar actions
    QAction * fullscreen_action;
    QAction * com_action;
    QAction * render_mode_action;
    QAction * reset_action;
    QAction * fit_particles_action;
    QAction * toggle_trans_action;
    QAction * toggle_grid_action;
    QAction * particles_form_action;
    QAction * options_form_action;
    QAction * toggle_play_action;
    QAction * screenshot_action;
    QAction * movie_form_action;
    QAction * reload_action;
    QAction * toggle_gas_action;
    QAction * auto_screenshot_action;
    QAction * auto_gl_screenshot_action;
    QAction * toggle_osd_action;
    QAction * toggle_camera_action;
    QAction * toggle_rotation_screen_action;
    // status bar
    QStatusBar * status_bar;
    // Projection ortho/perps
    QAction * proj_action;
    // GL actions
    QAction * next_cmap_action;
    QAction * prev_cmap_action;
    QAction * reverse_cmap_action;
    QAction * dynamic_cmap_action;
    // auto rotate around SCREEN axis
    QAction * rotatex_action, * rotatey_action,* rotatez_action,
    * rotatexr_action, * rotateyr_action,* rotatezr_action;
    // auto rotate around SCENE/Object axis
    QAction * rotateu_action, * rotatev_action,* rotatew_action,
    * rotateur_action, * rotatevr_action,* rotatewr_action;
    // auto translate
    QAction * transx_action, * transy_action, * transz_action;
    // Z sorting
    QAction * zsorting_action;
    // OpenGL setting dialog box
    //DlgGLSettings * box_gl_settings;
    // Forms
    FormObjectControl * form_o_c;
    FormAbout         * form_about;
    FormHelp          * form_help;
    FormScreenshot    * form_sshot;
    FormSelectPart    * form_spart;
    FormOptions       * form_options;
    FormConnect       * form_connect;
    // DockWidgets
    QDockWidget *dock_o_c, * dock_options;
    //------------------------------
    // Objects
    //------------------------------
    GLWindow * gl_window;

    void parseNemoParameters();
    GlobalOptions * store_options;
    PluginsManage * plugins;
    ParticlesObjectVector pov,pov2;
    ComponentRangeVector * crv;
    void listObjects(ParticlesObjectVector& ppov);
    void setDefaultParamObject(ParticlesObjectVector&);
    // NEMO parameters
    std::string snapshot,select,server, s_times,version;
    bool keep_all;
    // SnapshotInterface
    SnapshotInterface * current_data;
    // Playing snapshot/network
    bool play_animation, play;
    QTimer * play_timer;
    LoadingThread * loading_thread;
    // auto rotate timers arounf SCREEN axis
    QTimer * auto_rotx_timer, * auto_roty_timer, * auto_rotz_timer;
    QTimer * auto_transx_timer, * auto_transy_timer, * auto_transz_timer;
    // auto rotate timers arounf SCENE/Object axis
    QTimer * auto_rotu_timer, * auto_rotv_timer, * auto_rotw_timer;
    // bench
    QTimer * bench_gup_timer, * bench_nframe_timer;
    int total_frame;
    // divers stuff
    void initVariables();
    
    // timers
    void startTimers();
    // mouse and keys and zoom
    bool is_key_pressed;
    bool is_mouse_pressed;
    bool bestzoom;
    // mutex
    QMutex mutex_loading;
    QMutex  * mutex_data;
    QMutex sel_part;
    // user selection class
    UserSelection * user_select;
    bool reload;
    // Colormap
    Colormap * colormap;
    // Camera
    Camera * camera;
    // Selected Physical quantitie from the command line
    int selphys;
    signals:
    void endOfSnapshot();
};
}

#endif
