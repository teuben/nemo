// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
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

//#include <QGLWidget>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QLayout>
#include <assert.h>
#include <iomanip>
#include <sstream>
#include "glwindow.h"
#include "mainwindow.h"
#include "pluginsmanage.h"
#include "particlesselectrange.h"
#include "snapshotinterface.h"
#include "globjectparticles.h"
#include "snapshotnetwork.h"

#include "ftmio.h"
#include "nemo.h"
namespace glnemo {
#define ICONSIZE 25
  
// -----------------------------------------------------------------------------
// MainWindow constructor                                                       
// -----------------------------------------------------------------------------
MainWindow::MainWindow(std::string _ver)
{
  user_select = NULL;
  version = _ver;
  mutex_data = new QMutex(QMutex::Recursive); // Recursive: a thread can lock a mutex more than
                                              // once time, but mustunlock it as much as it    
                                              // has been locked                               
   
  // set Windows color
//  QPalette mp;
//  mp.setColor(QPalette::Window,QColor(224,212,247));
//  setPalette(mp);
  status_bar = statusBar();
  
  // Plugins
  plugins = new PluginsManage();

  // parse parameters
  parseNemoParameters();

  // Camera
  camera = new  Camera();  
  // ------- openGL object ---------
  gl_window = new glnemo::GLWindow(this,store_options,mutex_data, camera);
  camera->init(GlobalOptions::RESPATH.toStdString()+"/camera/path_01");
  // colormap object
  colormap  = new Colormap(store_options);
  
  // ----- build GUI ------------
  createForms();
  createDockWindows();
  createActions();
  createMenus();
  createToolBars();
  createStatusBar();

  setWindowTitle(tr("glnemo2"));
  setWindowIcon(QIcon(GlobalOptions::RESPATH+"/images/glnemo2.png"));

  // create SIGNAL/SLOTS connexions
  connect(gl_window, SIGNAL(sigKeyMouse(const bool, const bool)),
          this,    SLOT(pressedKeyMouse(const bool, const bool)));
  connect(form_o_c,SIGNAL(objectSettingsChanged()),gl_window,SLOT(updateGL()));
  connect(form_o_c,SIGNAL(objectUpdateVel(const int)),gl_window,SLOT(updateVel(const int)));
  connect(form_o_c,SIGNAL(objectUpdate()),gl_window,SLOT(update()));
  connect(form_o_c,SIGNAL(textureObjectChanged(const int, const int)),
          gl_window,SLOT(setTextureObject(const int, const int)));
  connect(form_o_c,SIGNAL(changeBoundaryPhys(const int)),
          gl_window,SLOT(updateBoundaryPhys(const int)));
  connect(form_o_c,SIGNAL(gazAlphaObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  connect(form_o_c,SIGNAL(gazSizeObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  connect(form_o_c,SIGNAL(gazColorObjectChanged(const int)),
          gl_window,SLOT(updateColorVbo(const int)));
  connect(form_o_c,SIGNAL(densityProfileObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  connect(form_o_c,SIGNAL(updateIpvs(int)),
          gl_window,SLOT(updateIpvs(int)));
  connect(form_o_c,SIGNAL(updateIpvs(int)),
          this,SLOT(updateIpvs(int)));
  // Interactive select/ glselect connections
  connect(gl_window->gl_select, SIGNAL(updatePareticlesSelected(const int)),
          form_options,SLOT(updateParticlesSelect(const int)));
  connect(form_options,SIGNAL(select_and_zoom(const bool)),
          gl_window->gl_select,SLOT(setZoom(bool)));
  connect(form_options,SIGNAL(save_selected()),this,SLOT(saveIndexList()));
  connect(form_options,SIGNAL(create_obj_selected()),this,SLOT(createObjFromIndexList()));
  // Camera
  connect(camera,SIGNAL(updateGL()),gl_window,SLOT(updateGL()));
  connect(form_options,SIGNAL(setCamDisplay(bool,bool)),camera,SLOT(setCamDisplay(bool,bool)));
  connect(form_options,SIGNAL(setSplineParam(int,double)),camera,SLOT(setSplineParam(int,double)));
  connect(form_options,SIGNAL(startStopPlay()),camera,SLOT(startStopPlay()));
  // colormap
  connect(colormap,SIGNAL(newColorMap()),gl_window,SLOT(changeColorMap()));
  //connect(colormap,SIGNAL(newColorMap()),gl_window,SLOT(reverseColorMap()));
  connect(colormap,SIGNAL(newColorMap()),form_o_c,SLOT(changeColorMap()));
  connect(form_o_c,SIGNAL(nextColorMap()),colormap,SLOT(next()));
  connect(form_o_c,SIGNAL(prevColorMap()),colormap,SLOT(prev()));
  connect(form_o_c,SIGNAL(constantColorMap(bool)),colormap,SLOT(constant(bool)));
  connect(form_o_c,SIGNAL(reverseColorMap(bool)),colormap,SLOT(reverse(bool)));
  // options play tab
  connect(form_options,SIGNAL(playPressed()),this,SLOT(actionPlay()));
  connect(this,SIGNAL(endOfSnapshot()),form_options,SLOT(on_play_pressed()));
  // leaveEvent to pass focus to gl_window
  connect(form_o_c,SIGNAL(leaveEvent()),gl_window,SLOT(setFocus()));
  connect(form_options,SIGNAL(leaveEvent()),gl_window,SLOT(setFocus()));
  // options grid tab
  connect(form_options,SIGNAL(update_grid()),gl_window,SLOT(updateGrid()));
  connect(form_options,SIGNAL(rebuild_grid()),gl_window,SLOT(rebuildGrid()));
  // options osd tab
  connect(form_options,SIGNAL(update_osd(bool)),this,SLOT(updateOsd(bool)));
  connect(form_options,SIGNAL(update_osd_font()),gl_window,SLOT(changeOsdFont()));
  connect(form_options,SIGNAL(update_gl()),gl_window,SLOT(updateGL()));
  // options GL colorbar tab
  connect(form_options,SIGNAL(update_gcb_font()),gl_window->gl_colorbar,SLOT(updateFont()));
  
  // --------- init some stuffs
  initVariables();
  startTimers();
  loading_thread = NULL;
  current_data = NULL;
  reload = false;
  //mainLayout->addWidget(qgl);
  setCentralWidget(gl_window);
}
// -----------------------------------------------------------------------------
// Start                                                                        
// -----------------------------------------------------------------------------
void MainWindow::start(std::string shot)
{
  const int wsize=getiparam((char *)"wsize");
  const int hsize=getiparam((char* )"hsize");
  if ( hasvalue((char *) "shot_ext") ) {
    store_options->base_frame_ext=getparam((char* )"shot_ext");
  }
  // try to load a snapshot
  user_select = new UserSelection();
  if (hasvalue((char*)"in")) {
    bool exist=hasvalue((char*)"select");
    current_data = plugins->getObject(snapshot);
    if (current_data) {
      connect(current_data,SIGNAL(stringStatus(QString)),status_bar, SLOT(showMessage(QString)));
      current_data->part_data->setIpvs(selphys);
      if (! exist ) {
        current_data->initLoading(store_options); 
        if (shot == "") interactiveSelect("",true);
      }
      else 
  	loadNewData(select,store_options->select_time,keep_all, false,true);
      if (!store_options->rho_exist) {
        store_options->render_mode = 0;
      }
      
      gl_window->updateGL();
    }
  }
  else if (hasvalue((char*)"server")) {
    std::string ip=getparam((char*)"server");
    int port=getiparam((char*)"port");
    bool vel=getbparam((char *)"vel");
    std::string select=getparam((char *)"select");
    if (select.length() == 0) {
        actionMenuFileConnect2(ip, port, vel, false, true);
    }
    else {
        if(actionMenuFileConnect2(ip, port, vel, false, false)) {
            crv = current_data->getSnapshotRange();
            selectPart(select, true);
        }
    }
  }
  if (shot != "" && play) {
      store_options->enable_gui=false;
      store_options->auto_play_screenshot=true;
      store_options->frame_width=wsize;
      store_options->frame_height=hsize;
      store_options->base_frame_name=QString(shot.c_str());
      startAutoScreenshot();
    } else {
      if (shot != "") {
        // get suffix
        QString suffix = ((QString(shot.c_str())).section('.', -1)).toUpper(); // 
        if (suffix != "PNG" && suffix != "JPG" && suffix != "JPEG") {
          // add selected extension as suffix
          shot = shot + "." + store_options->base_frame_ext.toStdString();
        }
        takeScreenshot(wsize,hsize,shot);
      }
    }
  
  if (play) {
      actionPlay(); // start playing time step
  }
  
  gl_window->setFocus();
  //actionMenuFileConnect();
}
// -----------------------------------------------------------------------------
// MainMainWindow destructor                                                    
// -----------------------------------------------------------------------------
MainWindow::~MainWindow()
{
  std::cerr << ">>  MainWindow::~MainWindow()\n";
  std::cerr << "delete form object control\n";
  delete form_o_c;
  std::cerr << "del store_options\n";
  delete store_options;
  std::cerr << "del plugins \n";
  delete plugins;
  std::cerr << "delete object\n";
  pov.clear();
  pov2.clear();
  std::cerr << "delete current data\n";
  delete current_data;
   std::cerr << "<<  MainWindow::~MainWindow()\n";
  if (user_select) delete user_select;

}
// -----------------------------------------------------------------------------
// create menus                                                                 
// -----------------------------------------------------------------------------
void MainWindow::createMenus()
{
  // file menu
  file_menu = menuBar()->addMenu(tr("&File"));
  file_menu->addAction(open_file_action);
  file_menu->addAction(connect_file_action);
  file_menu->addAction(print_file_action);
  file_menu->addSeparator();
  file_menu->addAction(quit_file_action);
  // help menu
  help_menu = menuBar()->addMenu(tr("&Help"));
  help_menu->addAction(doc_action);
  help_menu->addSeparator();
  help_menu->addAction(about_action);
  help_menu->addAction(about_qt);
}
// -----------------------------------------------------------------------------
// create tool bars                                                             
// -----------------------------------------------------------------------------
void MainWindow::createToolBars()
{
  icons_tool_bar = addToolBar(tr("Icons"));
  icons_tool_bar->addAction(fullscreen_action);
  icons_tool_bar->addAction(reset_action);
  icons_tool_bar->addAction(fit_particles_action);
  icons_tool_bar->addAction(toggle_grid_action);
  icons_tool_bar->addAction(particles_form_action);
  icons_tool_bar->addAction(options_form_action);
  //icons_tool_bar->addAction(toggle_trans_action);
  icons_tool_bar->addAction(toggle_play_action);
  icons_tool_bar->addAction(reload_action);
  icons_tool_bar->addAction(screenshot_action);
  //icons_tool_bar->addAction(toggle_gas_action);
  icons_tool_bar->addAction(print_file_action);
  icons_tool_bar->addAction(movie_form_action);
  //icons_tool_bar->addAction(com_action);
  icons_tool_bar->addAction(toggle_rotation_screen_action);
  
  QSize icons;
  icons.scale(ICONSIZE,ICONSIZE,Qt::KeepAspectRatio);
  icons_tool_bar->setIconSize(icons);
  icons_tool_bar->setAllowedAreas(Qt::AllToolBarAreas);
  icons_tool_bar->setOrientation(Qt::Horizontal);
  addToolBar(Qt::LeftToolBarArea,icons_tool_bar); // move to the left
  icons_tool_bar->setFocus();
}
// -----------------------------------------------------------------------------
// create status bars                                                           
// -----------------------------------------------------------------------------
void MainWindow::createStatusBar()
{
  statusBar()->showMessage(tr("Ready"));
}
// -----------------------------------------------------------------------------
// create Forms windows                                                         
// -----------------------------------------------------------------------------
void MainWindow::createForms()
{
  form_about  = new FormAbout(this);
  form_help   = new FormHelp(this);
  form_sshot  = new FormScreenshot(this);
  form_spart  = new FormSelectPart(this);
  form_o_c    = new FormObjectControl(this);
  form_options= new FormOptions(store_options,this);
  form_connect = new FormConnect(this);
  // sig/slot
  connect(form_sshot,SIGNAL(screenshot(const int, const int)),
	  this,SLOT(takeScreenshot(const int, const int)));
  connect(form_spart,SIGNAL(selectPart(const std::string,const bool)),
	  this,SLOT(selectPart(const std::string, const bool)));
  connect(form_options,SIGNAL(start_bench(const bool)),this,SLOT(startBench(const bool)));
  connect(form_connect, SIGNAL(newConnect(std::string, int, bool, bool, bool)), this, SLOT(actionMenuFileConnect2(std::string, int, bool, bool, bool)));
  // some init
  form_about->setVersion(QString(version.c_str()));
  form_o_c->init(mutex_data);
}
// -----------------------------------------------------------------------------
// create docking windows                                                       
// -----------------------------------------------------------------------------
void MainWindow::createDockWindows()
{

  // Docking area
  dock_o_c = new QDockWidget(tr("Object settings"),
                             this,Qt::Widget|
                                           Qt::WindowStaysOnTopHint|
                                           Qt::X11BypassWindowManagerHint);
  dock_o_c->setAllowedAreas(Qt::LeftDockWidgetArea|Qt::RightDockWidgetArea);
  // Create a scroll area to put the widget
  QScrollArea * scrollArea = new QScrollArea(this);
  scrollArea->setWidget(form_o_c);
  scrollArea->setWidgetResizable(true);
  // dock object
  dock_o_c->setWidget(scrollArea);
  dock_o_c->close();
  addDockWidget(Qt::RightDockWidgetArea,dock_o_c);

  //
  dock_options = new QDockWidget(tr("Options"),this,Qt::Widget|
		     Qt::WindowStaysOnTopHint|
		     Qt::X11BypassWindowManagerHint);
  dock_options->setAllowedAreas(Qt::AllDockWidgetAreas);
  // Create a scroll area to put the widget
  scrollArea = new QScrollArea(this);
  scrollArea->setWidget(form_options);
  scrollArea->setWidgetResizable(true);
  // dock object
  dock_options->setWidget(scrollArea);
  dock_options->close();
  addDockWidget(Qt::BottomDockWidgetArea,dock_options);

}
// -----------------------------------------------------------------------------
// create actions                                                               
// -----------------------------------------------------------------------------
void MainWindow::createActions()
{
  // ------ File menu actions --------
  // open
  open_file_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/folder_open.png"),tr("Open &File"),this);
  open_file_action->setShortcut(tr("Ctrl+O"));
  open_file_action->setStatusTip(tr("Open File from disk"));
  connect(open_file_action, SIGNAL(triggered()), this, SLOT(actionMenuFileOpen()));

  // connexion
  connect_file_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/connect.png"),tr("&Connect"),this);
  //connect_file_action->setShortcut(tr("Ctrl+O"));
  connect_file_action->setStatusTip(tr("Connect to a simulation server"));
  connect(connect_file_action, SIGNAL(triggered()), this, SLOT(actionMenuFileConnect()));
  // print
  print_file_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/fileprint.png"),tr("&Print"),this);
  print_file_action->setShortcut(tr("Ctrl+P"));
  print_file_action->setStatusTip(tr("Print OpenGL buffer"));
  connect(print_file_action, SIGNAL(triggered()), this, SLOT(actionEmpty()));

  // quit
  quit_file_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/exit.png"),tr("&Quit"),this);
  quit_file_action->setShortcut(tr("Ctrl+Q"));
  quit_file_action->setStatusTip(tr("Quit"));
  connect(quit_file_action, SIGNAL(triggered()), this, SLOT(actionQuit()));
  
  // ------- Help menu actions ---------
  // Documentation
  doc_action = new QAction(QIcon(""),tr("Documentation"),this);
  doc_action->setStatusTip(tr("Documentation"));
  connect(doc_action, SIGNAL(triggered()), form_help, SLOT(show()));
  // about glnemo
  about_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/glnemo2.png"),tr("About glnemo2"),this);
  about_action->setStatusTip(tr("About glnemo2"));
  connect(about_action, SIGNAL(triggered()), form_about, SLOT(show()));
  // about Qt
  about_qt = new QAction(tr("About &Qt"),this);
  about_qt->setStatusTip(tr("Show the Qt library's About box"));
  connect(about_qt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));

  // ------- options toolbar actions -------
  //
  // Full Screen
  fullscreen_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/window_fullscreen.png"),
                                  tr("Toggle Full Screen"),this);
  fullscreen_action->setShortcut(tr("F"));
  fullscreen_action->setStatusTip(tr("Toogle Full Screen"));
  connect( fullscreen_action, SIGNAL(activated()), this, SLOT(actionFullScreen()) );

  // Reset 
  reset_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/home-mdk.png"),tr("Reset to initial positions"),this);
  reset_action->setShortcut(tr("Ctrl+R"));
  reset_action->setStatusTip(tr("Reset to initial positions"));
  connect( reset_action, SIGNAL( activated() ), this, SLOT( actionReset() ) );

  // Center to COM 
  com_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/home-mdk.png"),tr("Center to COM"),this);
  com_action->setShortcut(tr("C"));
  com_action->setStatusTip(tr("Center to COM"));
  connect( com_action, SIGNAL( activated() ), this, SLOT( actionCenterToCom() ) );
  addAction(com_action);
  
  // render mode 
  render_mode_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/home-mdk.png"),tr("Rendering mode"),this);
  render_mode_action->setShortcut(tr("M"));
  render_mode_action->setStatusTip(tr("Rendering mode"));
  connect( render_mode_action, SIGNAL( activated() ), this, SLOT( actionRenderMode() ) );
  addAction(render_mode_action);

  // Fit all particles
  fit_particles_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/zoom-best-fit.png"),tr("Fit all particles on screen"),
                                     this);
  fit_particles_action->setShortcut(tr("Ctrl+A"));
  fit_particles_action->setStatusTip(tr("Fit all particles on screen"));
  connect( fit_particles_action, SIGNAL( activated() ), this, SLOT(actionBestZoom()) );

  // Grid 
  toggle_grid_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/grid.png"),tr("Toggle Grid"),this);
  toggle_grid_action->setShortcut(tr("G"));
  toggle_grid_action->setStatusTip(tr("Toggle Grid"));
  connect(toggle_grid_action, SIGNAL( activated() ),  this, SLOT(actionGrid()) );

  // Particles range & color
  particles_form_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/colors.png"),
                                      tr("set Particles range and color"),this);
  particles_form_action->setShortcut(tr("R"));
  particles_form_action->setStatusTip(tr("set Particles range and color"));
  connect( particles_form_action, SIGNAL( activated() ),
           this, SLOT(actionFormObjectControl() ) );

  // options
  options_form_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/options.png"),tr("Options dialog box"),this);
  options_form_action->setShortcut(tr("O"));
  options_form_action->setStatusTip(tr("Options dialog box"));
  connect(options_form_action, SIGNAL( activated() ), this, SLOT( actionFormOptions() ) );

  // Translation
  toggle_trans_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/move.png"),tr("Toggle translation"),this);
  toggle_trans_action->setShortcut(tr("T"));
  toggle_trans_action->setStatusTip(tr("Toggle translation"));
  connect(toggle_trans_action, SIGNAL( activated() ), this, SLOT(actionEmpty()) );

  // play simulation
  toggle_play_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/player_play.png"),tr("Play next snapshot"),this);
  toggle_play_action->setShortcut(tr("p"));
  toggle_play_action->setStatusTip(tr("Play next snapshot"));
  //connect(toggle_play_action, SIGNAL( activated() ), this, SLOT(actionPlay()) );
  connect(toggle_play_action, SIGNAL( activated() ), form_options, SLOT(on_play_pressed()));
  
  // reload
  reload_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/reload.png"),tr("Reload snaphot"),this);
  reload_action->setShortcut(tr("l"));
  reload_action->setStatusTip(tr("Reload snapshot"));
  connect(reload_action, SIGNAL( activated() ), this, SLOT(actionReload()) );
  
  // screenshot
  screenshot_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/camera.png"),tr("Take a screenshot"),this);
  screenshot_action->setShortcut(tr("S"));
  screenshot_action->setStatusTip(tr("Take a screenshot"));
  connect(screenshot_action, SIGNAL( activated() ), this, SLOT(actionScreenshot()) );
  
  // automatic screenshot during Play
  auto_screenshot_action = new QAction(QIcon(""),tr("Take a screenshot during play event"),this);
  auto_screenshot_action->setShortcut(tr("ctrl+S"));
  auto_screenshot_action->setStatusTip(tr("Take a screenshot during play event"));
  connect(auto_screenshot_action, SIGNAL( activated() ), this, SLOT(actionAutoScreenshot()) );
  addAction(auto_screenshot_action);

  // automatic screenshot during GL event
  auto_gl_screenshot_action = new QAction(QIcon(""),tr("Take a screenshot during OpenGL event"),this);
  auto_gl_screenshot_action->setShortcut(tr("ctrl+G"));
  auto_gl_screenshot_action->setStatusTip(tr("Take a screenshot during OpenGL event"));
  connect(auto_gl_screenshot_action, SIGNAL(activated()), this, SLOT(actionGLAutoScreenshot()) );
  addAction(auto_gl_screenshot_action);

  // toggle OSD 
  toggle_osd_action = new QAction(this);
  toggle_osd_action->setShortcut(tr("alt+t"));
  connect( toggle_osd_action, SIGNAL( activated() ), this, SLOT( actionToggleOsd() ) );
  addAction(toggle_osd_action);
  
  // print
  print_file_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/fileprint.png"),tr("Print OpenGL window"),this);
  print_file_action->setShortcut(tr(""));
  print_file_action->setStatusTip(tr("Print OpenGL window"));
  connect(print_file_action, SIGNAL( activated() ), this, SLOT(actionPrint()) );
  
  // movie
  movie_form_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/video_section.png"),tr("Make a movie"),this);
  movie_form_action->setShortcut(tr(""));
  movie_form_action->setStatusTip(tr("Make a movie"));
  connect(movie_form_action, SIGNAL( activated() ), this, SLOT(actionEmpty()) );
  // Next colormap
  next_cmap_action = new QAction(this);
  next_cmap_action->setShortcut(tr("Alt+Shift+n"));
  connect( next_cmap_action, SIGNAL( activated() ), colormap, SLOT( next() ) );
  addAction(next_cmap_action);
  // prev colormap
  prev_cmap_action = new QAction(this);
  prev_cmap_action->setShortcut(tr("Alt+Shift+p"));
  connect( prev_cmap_action, SIGNAL( activated() ), colormap, SLOT( prev() ) );
  addAction(prev_cmap_action);
  // reverse colormap
  reverse_cmap_action = new QAction(this);
  reverse_cmap_action->setShortcut(tr("Alt+Shift+i"));
  connect( reverse_cmap_action, SIGNAL( activated() ), colormap, SLOT( reverse() ) );
  addAction(reverse_cmap_action);
  // constant colormap
  dynamic_cmap_action = new QAction(this);
  dynamic_cmap_action->setShortcut(tr("Alt+Shift+c"));
  connect( dynamic_cmap_action, SIGNAL( activated() ), colormap, SLOT( constant() ) );
  addAction(dynamic_cmap_action);

  // Z sorting
  zsorting_action = new QAction(this);
  zsorting_action->setShortcut(tr("Z"));
  connect( zsorting_action, SIGNAL( activated() ), this, SLOT( actionZSorting() ) );
  addAction(zsorting_action);

  // Toggle rotation screen
  toggle_rotation_screen_action = new QAction(QIcon(GlobalOptions::RESPATH+"/images/fitscreen.png"),tr("Toggle rotation mode around axes screen/world"),this);
  toggle_rotation_screen_action->setShortcut(tr("Ctrl+L"));
  connect( toggle_rotation_screen_action, SIGNAL( activated() ), gl_window, SLOT( toggleRotateScreen()) );
  addAction(toggle_rotation_screen_action);
  
  // Auto rotate around X 
  rotatex_action = new QAction(this);
  rotatex_action->setShortcut(tr("Ctrl+X"));
  connect( rotatex_action, SIGNAL( activated() ), this, SLOT( actionRotateX() ) );
  addAction(rotatex_action);
  
  // Auto rotate reverse around X 
  rotatexr_action = new QAction(this);
  rotatexr_action->setShortcut(tr("Ctrl+Shift+X"));
  connect( rotatexr_action, SIGNAL( activated() ), this, SLOT( actionRotateRX() ) );
  addAction(rotatexr_action);
  
  // Auto rotate around Y 
  rotatey_action = new QAction(this);
  rotatey_action->setShortcut(tr("Ctrl+Y"));
  connect( rotatey_action, SIGNAL( activated() ), this, SLOT( actionRotateY() ) );
  addAction(rotatey_action);
  
  // Auto rotate reverse around Y 
  rotateyr_action = new QAction(this);
  rotateyr_action->setShortcut(tr("Ctrl++Shift+Y"));
  connect( rotateyr_action, SIGNAL( activated() ), this, SLOT( actionRotateRY() ) );
  addAction(rotateyr_action);
  
  // Auto rotate around Z 
  rotatez_action = new QAction(this);
  rotatez_action->setShortcut(tr("Ctrl+Z"));
  connect( rotatez_action, SIGNAL( activated() ), this, SLOT( actionRotateZ() ) );
  addAction(rotatez_action);
  
  // Auto rotate reverse around Z 
  rotatezr_action = new QAction(this);
  rotatezr_action->setShortcut(tr("Ctrl+Shift+Z"));
  connect( rotatezr_action, SIGNAL( activated() ), this, SLOT( actionRotateRZ() ) );
  addAction(rotatezr_action);
  
  // Auto rotate around U 
  rotateu_action = new QAction(this);
  rotateu_action->setShortcut(tr("Ctrl+U"));
  connect( rotateu_action, SIGNAL( activated() ), this, SLOT( actionRotateU() ) );
  addAction(rotateu_action);
  
  // Auto rotate reverse around U
  rotateur_action = new QAction(this);
  rotateur_action->setShortcut(tr("Ctrl+Shift+U"));
  connect( rotateur_action, SIGNAL( activated() ), this, SLOT( actionRotateRU() ) );
  addAction(rotateur_action);
  
  // Auto rotate around V 
  rotatev_action = new QAction(this);
  rotatev_action->setShortcut(tr("Ctrl+V"));
  connect( rotatev_action, SIGNAL( activated() ), this, SLOT( actionRotateV() ) );
  addAction(rotatev_action);
  
  // Auto rotate reverse around V 
  rotatevr_action = new QAction(this);
  rotatevr_action->setShortcut(tr("Ctrl+Shift+V"));
  connect( rotatevr_action, SIGNAL( activated() ), this, SLOT( actionRotateRV() ) );
  addAction(rotatevr_action);
  
  // Auto rotate around W 
  rotatew_action = new QAction(this);
  rotatew_action->setShortcut(tr("Ctrl+W"));
  connect( rotatew_action, SIGNAL( activated() ), this, SLOT( actionRotateW() ) );
  addAction(rotatew_action);
  
  // Auto rotate reverse around W
  rotatewr_action = new QAction(this);
  rotatewr_action->setShortcut(tr("Ctrl+Shift+W"));
  connect( rotatewr_action, SIGNAL( activated() ), this, SLOT( actionRotateRW() ) );
  addAction(rotatewr_action);
  
  // Auto translate along X
  transx_action = new QAction(this);
  transx_action->setShortcut(tr("Alt+X"));
  connect( transx_action, SIGNAL( activated() ), this, SLOT( actionTranslateX() ) );
  addAction(transx_action);
  // Auto translate along Y
  transy_action = new QAction(this);
  transy_action->setShortcut(tr("Alt+Y"));
  connect( transy_action, SIGNAL( activated() ), this, SLOT( actionTranslateY() ) );
  addAction(transy_action);
  // Auto translate along Z
  transz_action = new QAction(this);
  transz_action->setShortcut(tr("Alt+Z"));
  connect( transz_action, SIGNAL( activated() ), this, SLOT( actionTranslateZ() ) );
  addAction(transz_action);
}
// -----------------------------------------------------------------------------
// interactiveSelect                                                            
// Lauch select particles dialog box                                            
// -----------------------------------------------------------------------------
void MainWindow::interactiveSelect(std::string _select, const bool first_snapshot)
{
  if (current_data) {
    if (!reload) {
      crv = current_data->getSnapshotRange();
    }
    form_spart->update(current_data,&current_data->crv_first,_select, first_snapshot);
    form_spart->show();
    //ComponentRange::list(crv);
  }
}
// -----------------------------------------------------------------------------
// selectPart                                                                   
// Slots connected to "select particles dialog box" to allow to load particles  
// according to the user's selection.                                           
// -----------------------------------------------------------------------------
void MainWindow::selectPart(const std::string _select, const bool first_snapshot)
{
  select = _select;
  store_options->select_part = select;
  if (reload && current_data) {// reload action requested   
    store_options->phys_max_glob = store_options->phys_min_glob = -1; // reset for colobar display
    current_data->close();     // close the current snapshot
    delete current_data;       // delete previous object    
    current_data = plugins->getObject(snapshot); // connect
    connect(current_data,SIGNAL(stringStatus(QString)),status_bar, SLOT(showMessage(QString)));
    current_data->initLoading(store_options);
    crv = current_data->getSnapshotRange();    
    //ComponentRange::list(crv);
    //ComponentRange::list(&current_data->crv_first);
  } else {
    actionReset();             // reset view if menu file open
  }
  current_data->setSelectPart(select);
  std::cerr << "MainWindow::selectPart store_options->select_time = " << store_options->select_time << "\n";
  loadNewData(select,store_options->select_time,  // load data
	      keep_all,true,first_snapshot);
}
// -----------------------------------------------------------------------------
// loadNewData                                                                  
// -----------------------------------------------------------------------------
void MainWindow::loadNewData(const std::string select,
                             const std::string s_time,
                             const bool keep_all, const bool interact,
                             const bool first)
{
  ParticlesObjectVector povold;
  store_options->select_part = select;
  if (keep_all) {;};
  if (current_data) {
    std::cerr << "MainWindow::loadNewData s_time = " << s_time << "\n";
    if (!interact) current_data->initLoading(store_options);
    // get snapshot component ranges
    //ComponentRangeVector * crv = current_data->getSnapshotRange();
    if (!interact) crv = current_data->getSnapshotRange();
    //crv = current_data->getSnapshotRange();
    assert(crv);
    assert(crv->size());
    ComponentRange::list(crv);
    user_select->setSelection(select,crv,&pov); // fill pov according user sel and crv
    
    // load from disk
    mutex_data->lock();
    QTime tbench;
    tbench.restart();
    
    if (current_data->nextFrame(user_select->getIndexesTab(),user_select->getNSel())) {
      qDebug("Time elapsed to load snapshot: %d s", tbench.elapsed()/1000);
      store_options->new_frame=true;
      mutex_data->unlock();
      listObjects(pov);
      //listObjects(pov2);
      if (reload) { // backup old pov2 properties to povold
        povold.clear();
        ParticlesObject::backupVVProperties(pov2,povold,pov.size());
        //std::cerr << "POVOLD list\n";
        //listObjects(povold);        
      }
      ParticlesObject::initOrbitsVectorPOV(pov);
      pov2 = pov;   // copy new pov object to pov2
      ParticlesObject::clearOrbitsVectorPOV(pov2); // clear orbits vectors if present
      if (reload) { // copy back povold properties to pov2
        //std::cerr << "POVOLD list\n";
        //listObjects(povold);
        //std::cerr << "POV2 list before\n";
        //listObjects(pov2);
        ParticlesObject::backupVVProperties(povold,pov2,pov.size());
        //std::cerr << "POV2 list after\n";
        //listObjects(pov2);
      }
      if (first) {
        if (store_options->auto_texture_size && !store_options->rho_exist) {
          store_options->texture_size = current_data->part_data->getMaxSize()/100.;
          if (store_options->texture_size>1.0) store_options->texture_size=1.0;
          //store_options->texture_size = pow(current_data->part_data->getMaxSize(),3)/(*(current_data->part_data->nbody));
          std::cerr << "Resampled Texture Size = "<< store_options->texture_size <<"\n";
        }
        setDefaultParamObject(pov2); // set some default parameter if first loading
      }
      
      form_o_c->update( current_data->part_data, &pov2,store_options);
      updateOsd();
      tbench.restart();
      
      if (interact && !reload && store_options->rho_exist) {
        store_options->render_mode = 2; // density mode
      }
      if (interact && !reload && !store_options->rho_exist) {
        store_options->render_mode = 0; // alpha blending accumulation mode
      }
      // 9 dec 2011 force to density mode if rho exist
      if (store_options->rho_exist) {
        store_options->render_mode = 2; // density mode
      } else {
        store_options->render_mode = 0; // alpha blending accumulation mode
      }
      if (! store_options->auto_render) {
        store_options->render_mode=0;
      }
      if (store_options->auto_com) {
        actionCenterToCom(false);
      }
      gl_window->update( current_data->part_data, &pov2,store_options);
      qDebug("Time elapsed to update GL with new data: %d s", tbench.elapsed()/1000);
      if (!reload && bestzoom) gl_window->bestZoomFit();
      statusBar()->showMessage(tr("Snapshot loaded."));
    }
    else {
      mutex_data->unlock();
    }
    // it may be necessary to delete "current_data"
  }
}
// -----------------------------------------------------------------------------
// startTimers                                                                  
// -----------------------------------------------------------------------------
void MainWindow::startTimers()
{
  play_timer = new QTimer(this);
  connect(play_timer, SIGNAL(timeout()), this, SLOT(playEvent()));
  auto_rotx_timer = new QTimer(this);
  connect(auto_rotx_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundX()));
  auto_roty_timer = new QTimer(this);
  connect(auto_roty_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundY()));
  auto_rotz_timer = new QTimer(this);
  connect(auto_rotz_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundZ()));
  auto_transx_timer = new QTimer(this);
  connect(auto_transx_timer, SIGNAL(timeout()), gl_window, SLOT(translateX()));
  auto_transy_timer = new QTimer(this);
  connect(auto_transy_timer, SIGNAL(timeout()), gl_window, SLOT(translateY()));
  auto_transz_timer = new QTimer(this);
  connect(auto_transz_timer, SIGNAL(timeout()), gl_window, SLOT(translateZ()));

  auto_rotu_timer = new QTimer(this);
  connect(auto_rotu_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundU()));
  auto_rotv_timer = new QTimer(this);
  connect(auto_rotv_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundV()));
  auto_rotw_timer = new QTimer(this);
  connect(auto_rotw_timer, SIGNAL(timeout()), gl_window, SLOT(rotateAroundW()));
  
  bench_gup_timer = new QTimer(this);
  connect(bench_gup_timer, SIGNAL(timeout()), gl_window, SLOT(updateGL()));
  bench_nframe_timer = new QTimer(this);
  connect(bench_nframe_timer, SIGNAL(timeout()),this, SLOT(updateBenchFrame()));

}
// -----------------------------------------------------------------------------
// List Objects                                                                 
// -----------------------------------------------------------------------------
void MainWindow::listObjects(ParticlesObjectVector& ppov)
{
  for (int i=0; i<(int)ppov.size();i++) {
    std::cerr << "---------------------------------------\n";
    std::cerr << "Object #["<<i<<"]: ";
    std::cerr << "#npart = "<<ppov[i].npart << " -- "
              << "first  = "<<ppov[i].first << " -- "
              << "last   = "<<ppov[i].last  << " -- "
              << "step   = "<<ppov[i].step  << "\n";
  }
}
// -----------------------------------------------------------------------------
// Set default parameters to all the object                                     
// -----------------------------------------------------------------------------
void MainWindow::setDefaultParamObject(ParticlesObjectVector & pov){
  for (int i=0; i<(int)pov.size();i++) {
    pov[i].setPartSize(store_options->psize);
    pov[i].setPart(store_options->show_points);
    pov[i].setGaz(store_options->show_poly);
    pov[i].setGazSize(store_options->texture_size);
    pov[i].setGazAlpha(store_options->texture_alpha*255);
    pov[i].setGazSizeMax(store_options->texture_size);
    pov[i].setVel(store_options->show_vel);
    if (store_options->phys_min_glob!=-1) {
      pov[i].setMinPhys(store_options->phys_min_glob);
    }
    if (store_options->phys_max_glob!=-1) {
      pov[i].setMaxPhys(store_options->phys_max_glob);
    }
    
  }
}
// -----------------------------------------------------------------------------
// parse Nemo parameters                                                        
// -----------------------------------------------------------------------------
void MainWindow::parseNemoParameters()
{
  // instantiate store_options object
  store_options = new GlobalOptions();

  // Initialyze NEMO parameters
  snapshot                = getparam((char *) "in");
  server                  = getparam((char *) "server");
  select                  = getparam((char *) "select");
  store_options->select_part = select;
  store_options->select_time = getparam((char *) "times");
  keep_all                = getbparam((char *) "keep_all");
  store_options->vel_req  = getbparam((char *) "vel");
  store_options->show_vel = getbparam((char *) "disp_vel");
  store_options->blending = getbparam((char *) "blending");
  store_options->dbuffer  = getbparam((char *) "dbuffer");
  store_options->show_grid= getbparam((char *) "grid");
  store_options->mesh_length=getdparam((char *) "mesh_size");
  store_options->nb_meshs = getiparam((char *) "nb_meshs");
  store_options->xy_grid  = getbparam((char *) "xyg");
  store_options->xz_grid  = getbparam((char *) "xzg");
  store_options->yz_grid  = getbparam((char *) "yzg");
  store_options->show_cube= getbparam((char *) "cube");
  // On screen Display
  store_options->show_osd = getbparam((char *) "osd");
  store_options->osd_time = getbparam((char *) "osdtime");
  store_options->osd_nbody= getbparam((char *) "osdnbody");
  store_options->osd_zoom = getbparam((char *) "osdzoom");
  store_options->osd_rot  = getbparam((char *) "osdrot");
  store_options->osd_trans= getbparam((char *) "osdtrans");
  store_options->osd_title= getbparam((char *) "osdtitle");
  store_options->osd_data_type = getbparam((char *) "osddata");
  if (hasvalue((char *) "osd_set_title") ) {
    store_options->osd_title_name = getparam((char *) "osd_set_title");
  }
  store_options->osd_font_size = getdparam((char *) "osdfs");
  // Color Bar
  store_options->gcb_enable      = getbparam((char *) "cb");
  store_options->gcb_logmode     = getbparam((char *) "cblog");
  store_options->gcb_orientation = getiparam((char *) "cbloc");
  store_options->gcb_ndigits     = getiparam((char *) "cbdigits");
  store_options->gcb_offset      = getiparam((char *) "cboffset");
  store_options->gcb_pwidth      = getdparam((char *) "cbpw");
  store_options->gcb_pheight     = getdparam((char *) "cbph");
  store_options->gcb_font_size   = getdparam((char *) "cbfs");
  
  store_options->perspective=getbparam((char *) "perspective");
  store_options->orthographic = !store_options->perspective;
  play                   = getbparam((char *) "play");
  store_options->init_glsl=getbparam((char *) "glsl");
  bestzoom           = getbparam((char *) "bestzoom");
  store_options->xrot     = getdparam((char *) "xrot");
  store_options->yrot     = getdparam((char *) "yrot");
  store_options->zrot     = getdparam((char *) "zrot");
  store_options->xtrans   = getdparam((char *) "xtrans");
  store_options->ytrans   = getdparam((char *) "ytrans");
  store_options->ztrans   = getdparam((char *) "ztrans");
  store_options->zoom     = getdparam((char *) "zoom");
  store_options->psize    = getdparam((char *) "psize");
  store_options->port     = getiparam((char *) "port");
  store_options->show_points= getbparam((char *) "point");
  store_options->show_poly= getbparam((char *) "texture");
  store_options->xmin     = getdparam((char *) "xmin");
  store_options->xmax     = getdparam((char *) "xmax");
  store_options->ymin     = getdparam((char *) "ymin");
  store_options->ymax     = getdparam((char *) "ymax");
  store_options->zmin     = getdparam((char *) "zmin");
  store_options->zmax     = getdparam((char *) "zmax");
  store_options->lmin     = getiparam((char *) "lmin");
  store_options->lmax     = getiparam((char *) "lmax");
  store_options->scale    = getdparam((char *) "scale");
  // auto rendering mode
  store_options->auto_render = getbparam((char *) "auto_render");
  
  // select physical quantity to display
  selphys = getiparam((char* )"selphys");  
  // min/max physical value
  if ( hasvalue((char *) "minphys") ) {
      store_options->phys_min_glob=getdparam((char *) "minphys");
      store_options->phys_local=false;
  } 
  if ( hasvalue((char *) "maxphys") ) {
      store_options->phys_max_glob=getdparam((char *) "maxphys");
      store_options->phys_local=false;
  }
  // color map
  store_options->colormap += getiparam((char *) "cmapindex");
  
  store_options->auto_com           =getbparam((char *) "com");
  // textures
  store_options->auto_texture_size  =getbparam((char *) "auto_ts");
  store_options->texture_size       =getdparam((char *) "texture_s");
  store_options->texture_alpha      =getdparam((char *) "texture_a");
  
  store_options->duplicate_mem = getbparam((char *) "smooth_gui");
  // ortho
  float range_ortho;
  if (store_options->orthographic) {
    range_ortho=getdparam((char *) "ortho_range");
  }
  if (store_options->port) {;} // do nothing (remove compiler warning)
  
  //                         finish NEMO
}
// -----------------------------------------------------------------------------
// initVariables()                                                              
void MainWindow::initVariables()
{
  static bool first=true;
  if (first) {
    play_animation=false;
    first=false;
  }
  is_key_pressed   = false;
  is_mouse_pressed = false;
  gl_window->setMouseRot(store_options->xrot,store_options->yrot,store_options->zrot);
}
// -----------------------------------------------------------------------------
// killPlayingEvent                                                             
// -----------------------------------------------------------------------------
void MainWindow::killPlayingEvent()
{
  if (play_animation) {
    play_animation = false;
    play_timer->stop();
    if (loading_thread) { // a thread is running........
      while ( ! loading_thread->wait()) {
        std::cerr << "AD_THREAD still running, waiting....\n";
      }
    }
  }
}
// -----------------------------------------------------------------------------
// actionMenuFileOpen()
//------------------------------------------------------------------------------
void MainWindow::actionMenuFileOpen()
{
  static QString menudir("");
  killPlayingEvent();       // wait the end of loading thread
  QString fileName = QFileDialog::getOpenFileName(this,tr("Select snapshot"),menudir);
  if (!fileName.isEmpty()) {
    menudir = fileName;
    snapshot = fileName.toStdString();
    SnapshotInterface * new_data = plugins->getObject(snapshot);
    if (new_data)  { // valid object
      mutex_loading.lock();     // protect area
      if (current_data)
        delete current_data;      // free memory                   
      current_data = new_data;  // link new_data   
      connect(current_data,SIGNAL(stringStatus(QString)),status_bar, SLOT(showMessage(QString)));
      current_data->part_data->setIpvs(selphys);
//       loadNewData("all","all",  // load data
//           keep_all,store_options->vel_req,true); //
      reload=false;
      bestzoom = true;
      store_options->select_part="";
      current_data->initLoading(store_options);
      interactiveSelect("",true);
      mutex_loading.unlock();   // release area                  
    }
  }
}
// -----------------------------------------------------------------------------
// actionMenuFileConnect()
//------------------------------------------------------------------------------
void MainWindow::actionMenuFileConnect()
{
  killPlayingEvent();       // wait the end of loading thread
  form_connect->show(); // show connect form
}
// -----------------------------------------------------------------------------
// actionMenuFileConnect2()
//------------------------------------------------------------------------------
bool MainWindow::actionMenuFileConnect2(std::string adresseIP, int Port, bool velocities, bool densities, bool interactivSelect)
{
  if (densities) {;} // remove compiler warning
  SnapshotNetwork * net = new SnapshotNetwork();
  std::string adresseip = adresseIP;
  int port = Port;
  SnapshotInterface * new_data;
  new_data = net->newObject(adresseip,port);
  if (new_data->isValidData()) {
    mutex_loading.lock();     // protect area
    delete current_data;      // free memory
    current_data = new_data;  // link new_data
    //     loadNewData("all","all",  // load data
    //           keep_all,store_options->vel_req,true); //
    store_options->vel_req = velocities;
    current_data->initLoading(store_options);
    reload=false;
    bestzoom = true;
    if (interactivSelect) interactiveSelect(select,true);
    mutex_loading.unlock();   // release area
    return true;
  }
  else {
    return false;
  }
}
// -----------------------------------------------------------------------------
// actionQuit()                                                                 
void MainWindow::actionQuit()
{
  if (current_data) {
    mutex_loading.lock();      // protect area
    killPlayingEvent();        // wait the end of loading thread
    mutex_loading.unlock();    // release area
  }
  if (store_options->enable_gui)
    close();  
  else
    //QCoreApplication::quit();
    exit(1);
}
// -----------------------------------------------------------------------------
// actionReload()                                                               
void MainWindow::actionReload()
{
  if (! current_data) return;
  mutex_loading.lock();      // protect area
  killPlayingEvent();        // wait the end of loading thread
  
//   current_data->close();     // close the current snapshot
//   delete current_data;
//   current_data = plugins->getObject(snapshot);
  
//   loadNewData(select,store_options->select_time,// load data
//         keep_all,store_options->vel_req,true); //
  reload=true;
  interactiveSelect(select);
  mutex_loading.unlock();    // release area
}
// -----------------------------------------------------------------------------
// actionFormObjectControl()                                                    
void MainWindow::actionFormObjectControl()
{
  bool show= ! dock_o_c->isVisible();
  if (show) dock_o_c->show();
  else      dock_o_c->close();
}
// -----------------------------------------------------------------------------
// actionFormOptions()                                                    
void MainWindow::actionFormOptions()
{
  bool show= ! dock_options->isVisible();
  if (show) dock_options->show();
  else      dock_options->close();
}
// -----------------------------------------------------------------------------
// actionFullScreen()                                                           
void MainWindow::actionFullScreen()
{
  static bool full_screen=true;

  if (full_screen)  {
    fullscreen_action->setIcon(QIcon(GlobalOptions::RESPATH+"/images/window_nofullscreen.png"));
    showFullScreen();
  }
  else {
    setWindowIcon(QIcon(GlobalOptions::RESPATH+"/images/glnemo2.png"));
    fullscreen_action->setIcon(QIcon(GlobalOptions::RESPATH+"/images/window_fullscreen.png"));
    showNormal();
  }
  full_screen=!full_screen;
  statusBar()->showMessage(tr("Ready"));
}
// -----------------------------------------------------------------------------
// actionBestactionBestZoom()                                                   
void MainWindow::actionBestZoom()
{
  gl_window->bestZoomFit();
}
// -----------------------------------------------------------------------------
// actionRenderMode()                                                           
void MainWindow::actionRenderMode()
{
  store_options->render_mode = (store_options->render_mode+1)%3;
  if (store_options->render_mode==1) store_options->render_mode=2; // giveup mode 1
  store_options->auto_render=true;
  gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// actionCenterToCom()                                                   
void MainWindow::actionCenterToCom(const bool ugl)
{
  double com[3] = {0., 0., 0.};
  int np=0;
  mutex_data->lock();
  //mutex_loading.lock();
  if (current_data ) {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // September 2009, 4th                               
    // Change "pov2" by "pov", seems to fix              
    // a bug regarding bad COM when interactive centering
    // while loading snapshot...                         
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i=0; i<(int)pov2.size();i++) {
	 // loop on all the objects
	const ParticlesObject * po = &(pov2[i]);        // object
	if (po->isVisible()) {                                   // is visible  
	  ParticlesData * part_data = current_data->part_data;// get its Data
	  // loop on all the particles of the object
	  for (int j  = 0; j  <  po->npart; j ++) {
	    np++;
	    int jndex= po->index_tab[j];
	    com[0] +=part_data->pos[jndex*3  ];
	    com[1] +=part_data->pos[jndex*3+1];
	    com[2] +=part_data->pos[jndex*3+2];
	  }
	}
    }
    store_options->xtrans = -(com[0]/np);
    store_options->ytrans = -(com[1]/np);
    store_options->ztrans = -(com[2]/np);
    if (ugl) gl_window->updateGL();
  }
  mutex_data->unlock();
  //mutex_loading.unlock();
}
// -----------------------------------------------------------------------------
// actionReset()                                                                
void MainWindow::actionReset()
{
  gl_window->resetView();
}
// -----------------------------------------------------------------------------
// actionGrid()                                                                 
void MainWindow::actionGrid()
{
  store_options->show_grid = !store_options->show_grid;
  form_options->update();
  gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// actionPrint()                                                                
void MainWindow::actionPrint()
{
  QPrinter printer;
  QPrintDialog *dlg = new QPrintDialog(&printer, this);
  if (dlg->exec() != QDialog::Accepted) return;
  gl_window->updateGL();
  QImage img=gl_window->grabFrameBuffer();
  QPainter painter( &printer );
  painter.drawImage(0,0,img);
}
// -----------------------------------------------------------------------------
// actionScreenshot()                                                           
void MainWindow::actionScreenshot()
{
  form_sshot->show();
}
// -----------------------------------------------------------------------------
// actionAutoScreenshot()                                                       
void MainWindow::actionAutoScreenshot()
{
  store_options->auto_play_screenshot = !store_options->auto_play_screenshot; 
}
// -----------------------------------------------------------------------------
// actionGLAutoScreenshot()                                                     
void MainWindow::actionGLAutoScreenshot()
{
  store_options->auto_gl_screenshot = !store_options->auto_gl_screenshot;
}
// -----------------------------------------------------------------------------
// actionToggleRotationScreen
void MainWindow::actionToggleRotationScreen()
{
}
// -----------------------------------------------------------------------------
// actionToggleOsd()
void MainWindow::actionToggleOsd()
{
  store_options->show_osd = !store_options->show_osd;
  form_options->update();
  gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// actionToggleCamera()
void MainWindow::actionToggleCamera()
{
  store_options->cam_mode = !store_options->cam_mode;
  //gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// actionZSorting()
void MainWindow::actionZSorting()
{
  store_options->zsort = !store_options->zsort;
  gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// startAutoScreenshot()                                                       
void MainWindow::startAutoScreenshot()
{
  std::ostringstream stm1;
  // string index
  stm1 << std::setfill('0')<< std::setw(5)<<store_options->frame_index;
  // frame name (jpg)
  std::string framename = store_options->base_frame_name.toStdString()+"."+stm1.str()+"."+
                          store_options->base_frame_ext.toStdString();
  std::cerr << "-----> framename = "<<framename<<"\n";
  //gl_window->rotateAroundAxis(1);
  std::string mess="Offscreen rendering : "+framename;
  statusBar()->showMessage(tr(mess.c_str()));
  takeScreenshot(store_options->frame_width,store_options->frame_height,framename);
  store_options->frame_index++;
}
// -----------------------------------------------------------------------------
// takeScreenshot()                                                             
void MainWindow::takeScreenshot(const int width, const int height,  std::string name)
{
    if ( width!=0 && height!=0) { // valid dimensions
        QSize size,sizegl;
        sizegl=gl_window->size();   // get the current Ogl windows's size
        if (width != -1) { // standard resolution or custom
            size.setWidth(width);
            size.setHeight(height);
        }
        else {           // screen resolution
            size  = sizegl;
        }
        gl_window->resizeOsd(size.width(),size.height());

        QRect geom = gl_window->geometry(); // save the current geometry of the GL's window
        std::cerr << "MainWindow::takeScreenshot call resizen";
        gl_window->resize(size.width(),size.height());
        // !!! activate the following line if you want to see OSD
#if 0
        gl_window->setGeometry(geom.x(),geom.y(),             // give to the widget the new size
                               size.width(),size.height());   // bc width() and height() are used by
        // renderText
#endif

        //std::cerr << "GLWINDOW width = " << gl_window->width() << "\n";
        gl_window->setFBO(true);                              // activate Frame Buffer Object
        gl_window->setFBOSize(size.width(),size.height());    // set the offscreen rendering size

        gl_window->updateGL();                                // draw in FBO

        QImage img=(((gl_window->grabFrameBufferObject()).mirrored()).rgbSwapped()); // convert FBO to img

        gl_window->resize(sizegl.width(),sizegl.height());    // revert to the previous Ogl windows's size
        gl_window->updateGL();

        // !!! activate the following line if you want to see OSDr
#if 0

        gl_window->setGeometry(geom.x(),geom.y(),
                               sizegl.width(),sizegl.height());
#else

#endif
        //gl_window->setFixedSize(sizegl);                   // revert to the previous Ogl Widget's size
        //
        if (name == "") { // interactive screenshot
            QString fileName = QFileDialog::getSaveFileName(this, tr("Save File"),"",
                                                            tr("Images (*.png *.jpg)"));
            if (!fileName.isEmpty()) {
                img.save(fileName);
                //resize(width(),height());
            }
        } else {          // screenshot from the command line
            int quality=-1;
            if (store_options->base_frame_ext=="jpg") {
                quality=95;
            }
            std::cerr << "takescreenshot name ="<<name<<"\n";
            std::cerr << "base_frame_ext="<<store_options->base_frame_ext.toStdString()<<"\n";
            //img.save(QString(name.c_str()),(store_options->base_frame_ext.toStdString()).c_str(),quality);
            img.save(QString(name.c_str()),0,quality);
        }
    }
}
// -----------------------------------------------------------------------------
// actionRotate around SCREEN axis
// -----------------------------------------------------------------------------
// actionRotateX() starts timer to rotate around screen x axis                                                       
void MainWindow::actionRotateX()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->ixrot = 1;
  if (rot)  auto_rotx_timer->start(20);
  else      auto_rotx_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRX() starts timer to rotate around screen x axis (reverse)                                   
void MainWindow::actionRotateRX()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->ixrot = -1;
  if (rot)  auto_rotx_timer->start(20);
  else      auto_rotx_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateY() starts timer to rotate around screen y axis              
void MainWindow::actionRotateY()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iyrot = 1;
  if (rot)  auto_roty_timer->start(20);
  else      auto_roty_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRY() starts timer to rotate around screen y axis (reverse)                                   
void MainWindow::actionRotateRY()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iyrot = -1;
  if (rot)  auto_roty_timer->start(20);
  else      auto_roty_timer->stop();
}

// -----------------------------------------------------------------------------
// actionRotateZ() starts timer to rotate around screen z axis                              
void MainWindow::actionRotateZ()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->izrot = 1;
  if (rot)  auto_rotz_timer->start(20);
  else      auto_rotz_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRZ() starts timer to rotate around screen z axis (reverse)                                   
void MainWindow::actionRotateRZ()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->izrot = -1;
  if (rot)  auto_rotz_timer->start(20);
  else      auto_rotz_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotate around SCENE/Object axis
// -----------------------------------------------------------------------------
// actionRotateU() starts timer to rotate around scene u axis                                                       
void MainWindow::actionRotateU()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iurot = 1;
  if (rot)  auto_rotu_timer->start(20);
  else      auto_rotu_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRU() starts timer to rotate around scene u axis (reverse)                                   
void MainWindow::actionRotateRU()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iurot = -1;
  if (rot)  auto_rotu_timer->start(20);
  else      auto_rotu_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateV() starts timer to rotate around scene v axis                                                       
void MainWindow::actionRotateV()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->ivrot = 1;
  if (rot)  auto_rotv_timer->start(20);
  else      auto_rotv_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRV() starts timer to rotate around scene v axis (reverse)                                   
void MainWindow::actionRotateRV()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->ivrot = -1;
  if (rot)  auto_rotv_timer->start(20);
  else      auto_rotv_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateW() starts timer to rotate around scene w axis                                                       
void MainWindow::actionRotateW()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iwrot = 1;
  if (rot)  auto_rotw_timer->start(20);
  else      auto_rotw_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateRW() starts timer to rotate around scene w axis (reverse)                                   
void MainWindow::actionRotateRW()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  store_options->iwrot = -1;
  if (rot)  auto_rotw_timer->start(20);
  else      auto_rotw_timer->stop();
}
// -----------------------------------------------------------------------------
// actionTranlateX()                                                            
void MainWindow::actionTranslateX()
{
  static bool trans=false;
  trans = !trans; // toggle translation
  if (trans)  auto_transx_timer->start(50);
  else        auto_transx_timer->stop();
}
// -----------------------------------------------------------------------------
// actionTranlateZ()                                                            
void MainWindow::actionTranslateY()
{
  static bool trans=false;
  trans = !trans; // toggle translation
  if (trans)  auto_transy_timer->start(50);
  else        auto_transy_timer->stop();
}
// -----------------------------------------------------------------------------
// actionTranlateZ()                                                            
void MainWindow::actionTranslateZ()
{
  static bool trans=false;
  trans = !trans; // toggle translation
  if (trans)  auto_transz_timer->start(50);
  else        auto_transz_timer->stop();
}
// -----------------------------------------------------------------------------
// actionPlay()                                                                 
void MainWindow::actionPlay()
{
  if ( ! current_data ) {
    QString message=tr("No Data loaded");
    QMessageBox::information( this,tr("Warning"),message,"Ok");
  }
  else {
    play_animation = !play_animation;
    if (play_animation) {
      if ( current_data->isEndOfData() ) {
        if (store_options->enable_gui) {
            std::cerr << "store_options->enable_gui.......\n";
            QMessageBox::information( this,tr("Warning"),
                                      current_data->endOfDataMessage(),"Ok");
            emit endOfSnapshot();
        }
        else {
            //killPlayingEvent();
            actionQuit();
            close();
        }
        //play_animation = false;
        //emit endOfSnapshot();
      }
      else { // start timer
        play_timer->start( 20 );
      }
    }
    else {
      play_timer->stop();
      gl_window->updateGL(); // flush openGL buffer
      //if (! anim_engine->record->isActivated()) 
      //  glbox->setHudActivate(GLHudObject::Loading, FALSE);
    }
  }
}
// -----------------------------------------------------------------------------
// playEvent()                                                                  
void MainWindow::playEvent()
{
  mutex_loading.lock();
  if ( ! loading_thread) { // no active thread
    //pov = pov2; // modif orbits
    loading_thread = new LoadingThread(current_data,user_select,&pov,select,mutex_data,store_options);
    //connect(loading_thread,SIGNAL(finished()),this,SLOT(uploadNewFrame()));
    loading_thread->start();
  }
  else {
    //std::cerr << "loading_thread->isFinished() ?\n" << loading_thread->isFinished() << "\n";
    if (loading_thread->isFinished() && // loading complete            
        !is_key_pressed              && // no interactive user request 
        !is_mouse_pressed               // (mouse, keyboard)           
       ) {
      uploadNewFrame();
      delete loading_thread;
      loading_thread=NULL;
    }
    else { // could implement a blinking flag
    }
  }
  mutex_loading.unlock();
}
// -----------------------------------------------------------------------------
// uploadNewFrame                                                               
void MainWindow::uploadNewFrame()
{
  static bool first=true;
  if (loading_thread->isValidLoading()) {
    if (  current_data->getInterfaceType() == "List of Ftm"      ||
	  current_data->getInterfaceType() == "List of Gadget 2" ||
          current_data->getInterfaceType() == "List of PhiGRAPE" ||
          current_data->getInterfaceType() == "List of Nemo") { 
      //mutex_data->lock();
      //pov2=pov;
      ParticlesObject::initOrbitsVectorPOV(pov);
      ParticlesObject::copyVVkeepProperties(pov,pov2,user_select->getNSel()); 
      form_o_c->update( current_data->part_data, &pov2,store_options,false); // update Form
      //mutex_data->unlock();
    } else {
      //pov2=pov; // modif orbits
    }
    updateOsd();
#if 0
    if (store_options->rho_exist) {
      store_options->render_mode = 2; // density mode
    }
    if (!store_options->rho_exist) {
      store_options->render_mode = 0; // alpha blending accumulation mode
    }
#endif
    if (store_options->auto_com) {
       actionCenterToCom(false);
    } 
    gl_window->update( current_data->part_data, &pov2,store_options);
    if (first && bestzoom) {
      first=false;
      //gl_window->bestZoomFit();
    }
    //gl_window->bestZoomFit();
    if (store_options->auto_play_screenshot && !store_options->auto_gl_screenshot) {
      startAutoScreenshot();
    }
  }
  else {
    if ( current_data->isEndOfData()) {
      //std::cerr << "current_data is end of data\n";
      play_animation = false;
      play_timer->stop();
      if (store_options->enable_gui)
          QMessageBox::information( this,tr("Warning"),current_data->endOfDataMessage(),"Ok");
      else {
          mutex_loading.unlock();
          //killPlayingEvent();
          actionQuit();
      }
    }
  }
}
// -----------------------------------------------------------------------------
// pressedKeyMouse()                                                            
void MainWindow::pressedKeyMouse(const bool k, const bool m)
{
  is_key_pressed   = k;
  is_mouse_pressed = m;
}
// -----------------------------------------------------------------------------
// updateOsd()                                                                  
void MainWindow::updateOsd(bool ugl)
{
  if (current_data) {
    GlobalOptions * g = store_options;
    gl_window->setOsd(GLObjectOsd::Nbody,
                      *current_data->part_data->nbody,g->osd_nbody,false);
    gl_window->setOsd(GLObjectOsd::Time,
                      *current_data->part_data->timu,g->osd_time,false);
    
    std::string title = g->osd_title_name.toStdString();
    if (title=="") {
      title = current_data->getFileName();
    }
    gl_window->setOsd(GLObjectOsd::Title,
                      QString(title.c_str()),
                      g->osd_title,false);
    gl_window->setOsd(GLObjectOsd::Getdata,
                      QString((current_data->getInterfaceType()).c_str()),
                      g->osd_data_type,false);
    gl_window->setOsd(GLObjectOsd::Zoom,(const float) store_options->zoom,
                      g->osd_zoom,false);
    gl_window->setOsd(GLObjectOsd::Rot,(const float) store_options->xrot,
                      (const float) store_options->yrot,
                      (const float) store_options->zrot, g->osd_rot,false);
    gl_window->setOsd(GLObjectOsd::Trans,(const float) store_options->xtrans,
                      (const float) store_options->ytrans,
                      (const float) store_options->ztrans,g->osd_trans,false);
  }
  if (ugl) {
    gl_window->updateGL();
  }
  
}
// -----------------------------------------------------------------------------
// startBench()                                                                 
void MainWindow::startBench(const bool start)
{
  if (start) {
    total_frame=0;
    gl_window->resetFrame();
    bench_gup_timer->start(1);       // update GL every 5ms
    bench_nframe_timer->start(500);  // update tt frame every 500 ms
  }
  else {
    bench_gup_timer->stop();
    bench_nframe_timer->stop();
  }
}
// -----------------------------------------------------------------------------
// updateBenchFrame()                                                           
void MainWindow::updateBenchFrame()
{
  int frame=gl_window->getFrame();
  total_frame += frame;
  //form_options->updateFrame(1000*frame/500,total_frame);
  form_options->updateFrame(frame,total_frame);
  gl_window->resetFrame();
  
}
// -----------------------------------------------------------------------------
// createObjFromIndexList()                                                                
void MainWindow::createObjFromIndexList()
{
  std::vector <int> * list = gl_window->gl_select->getList();
  if (list->size()) { //&& current_data->part_data->id.size()>0) {
    std::vector <int> indexes;
    indexes.reserve(list->size());
    for (std::vector<int>::iterator i=list->begin(); i<list->end(); i++) {     
      if (current_data->part_data->id.size()>0) { // save by ids
        //indexes.push_back(current_data->part_data->id[*i]);
        indexes.push_back(*i);
      } else {
        indexes.push_back(*i);
      }
    }
    ParticlesObject * po = new ParticlesObject(ParticlesObject::Range); // new object
    po->buildIndexList(indexes);
    pov.push_back(*po);
    delete po;
    listObjects(pov);
    ParticlesObject::copyVVkeepProperties(pov,pov2,user_select->getNSel()); 
    form_o_c->update( current_data->part_data, &pov2,store_options,false); // update Form
    updateOsd();
    gl_window->update( current_data->part_data, &pov2,store_options,false);
  }
}

// -----------------------------------------------------------------------------
// saveIndexList                                                                
void MainWindow::saveIndexList()
{
  std::vector <int> * list = gl_window->gl_select->getList();
  if (list->size()) {// && current_data->part_data->id.size()>0) {
    QString fileName = QFileDialog::getSaveFileName(this,tr("Save list of indexes"));
    if (!fileName.isEmpty()) {
      QFile file(fileName);
      if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << "#glnemo_index_list\n";
        for (std::vector<int>::iterator i=list->begin(); i<list->end(); i++) {   
          //std::cerr << "iterator  *i="<< *i <<"\n";           
          if (current_data->part_data->id.size()>0) { // save by ids
            out << current_data->part_data->id[*i] << "\n";
          } else {                                    // save by indexes
            out << *i << "\n";
          }
        }
        file.close();
      }
    }
  }
}
} // glnemo namespace }
