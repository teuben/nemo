// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2008                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           P�le de l'Etoile, site de Ch�teau-Gombert                         
//           38, rue Fr�d�ric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include <QtGui>

#include <QGLWidget>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QLayout>
#include <assert.h>
#include <iomanip>
#include <sstream>
#include "mainwindow.h"
#include "glwindow.h"
#include "pluginsmanage.h"
#include "particlesselectrange.h"
#include "snapshotinterface.h"

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
  QPalette mp;
  mp.setColor(QPalette::Window,QColor(224,212,247));
  setPalette(mp);
  
  // Plugins
  plugins = new PluginsManage();

  // parse parameters
  parseNemoParameters();

  // Camera
  camera = new  Camera();  
  // ------- openGL object ---------
  gl_window = new glnemo::GLWindow(this,store_options,mutex_data, camera);
  camera->init(":/camera/path_01");
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
  setWindowIcon(QIcon(":/images/glnemo2.png"));
  //readColormaps();
  // create SIGNAL/SLOTS connexions
  connect(gl_window, SIGNAL(sigKeyMouse(const bool, const bool)),
          this,    SLOT(pressedKeyMouse(const bool, const bool)));
  connect(form_o_c,SIGNAL(objectSettingsChanged()),gl_window,SLOT(updateGL()));
  connect(form_o_c,SIGNAL(objectUpdateVel(const int)),gl_window,SLOT(updateVel(const int)));
  connect(form_o_c,SIGNAL(objectUpdate()),gl_window,SLOT(update()));
  connect(form_o_c,SIGNAL(textureObjectChanged(const int, const int)),
          gl_window,SLOT(setTextureObject(const int, const int)));
  connect(form_o_c,SIGNAL(gazAlphaObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  connect(form_o_c,SIGNAL(gazSizeObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  connect(form_o_c,SIGNAL(gazColorObjectChanged(const int)),
          gl_window,SLOT(updateColorVbo(const int)));
  connect(form_o_c,SIGNAL(densityProfileObjectChanged(const int)),
          gl_window,SLOT(updateVbo(const int)));
  // Interactive select/ glselect connections
  connect(gl_window->gl_select, SIGNAL(updatePareticlesSelected(const int)),
          form_options,SLOT(updateParticlesSelect(const int)));
  connect(form_options,SIGNAL(select_and_zoom(const bool)),
          gl_window->gl_select,SLOT(setZoom(bool)));
  connect(form_options,SIGNAL(save_selected()),this,SLOT(saveIndexList()));
  // Camera
  connect(camera,SIGNAL(updateGL()),gl_window,SLOT(updateGL()));
  connect(form_options,SIGNAL(setCamDisplay(bool,bool)),camera,SLOT(setCamDisplay(bool,bool)));
  connect(form_options,SIGNAL(setSplineParam(int,double)),camera,SLOT(setSplineParam(int,double)));
  connect(form_options,SIGNAL(startStopPlay()),camera,SLOT(startStopPlay()));
  // colormap
  connect(colormap,SIGNAL(newColorMap()),gl_window,SLOT(changeColorMap()));
  connect(colormap,SIGNAL(newColorMap()),form_o_c,SLOT(changeColorMap()));
  // options play tab
  connect(form_options,SIGNAL(playPressed()),this,SLOT(actionPlay()));
  connect(this,SIGNAL(endOfSnapshot()),form_options,SLOT(on_play_pressed()));
  // leaveEvent to pass focus to gl_window
  connect(form_o_c,SIGNAL(leaveEvent()),gl_window,SLOT(setFocus()));
  connect(form_options,SIGNAL(leaveEvent()),gl_window,SLOT(setFocus()));
  
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
  // try to load a snapshot
  user_select = new UserSelection();
  if (hasvalue((char*)"in")) {
    bool exist=hasvalue((char*)"select");
    current_data = plugins->getObject(snapshot);
    if (current_data) {
      if (! exist ) {
        current_data->initLoading(store_options->vel_req, s_times); 
        if (shot == "") interactiveSelect();
      }
      else 
  	loadNewData(select,s_times,keep_all,store_options->vel_req, false);
      if (!store_options->rho_exist) {
        store_options->render_mode = 0;
      }
      gl_window->updateGL();
    }
  }
  if (play) actionPlay(); // start playing time step
  if (shot != "") takeScreenshot(wsize,hsize,shot);
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
  icons_tool_bar->addAction(toggle_trans_action);
  icons_tool_bar->addAction(toggle_play_action);
  icons_tool_bar->addAction(reload_action);
  icons_tool_bar->addAction(screenshot_action);
  //icons_tool_bar->addAction(toggle_gas_action);
  icons_tool_bar->addAction(print_file_action);
  icons_tool_bar->addAction(movie_form_action);
  //icons_tool_bar->addAction(com_action);
  
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
  form_sshot  = new FormScreenshot(this);
  form_spart  = new FormSelectPart(this);
  form_o_c    = new FormObjectControl(this);
  form_options= new FormOptions(store_options,this);
  // sig/slot
  connect(form_sshot,SIGNAL(screenshot(const int, const int)),
	  this,SLOT(takeScreenshot(const int, const int)));
  connect(form_spart,SIGNAL(selectPart(const std::string)),
	  this,SLOT(selectPart(const std::string)));
  connect(form_options,SIGNAL(start_bench(const bool)),this,SLOT(startBench(const bool)));
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
  open_file_action = new QAction(QIcon(":/images/folder_open.png"),tr("Open &File"),this);
  open_file_action->setShortcut(tr("Ctrl+O"));
  open_file_action->setStatusTip(tr("Open File from disk"));
  connect(open_file_action, SIGNAL(triggered()), this, SLOT(actionMenuFileOpen()));

  // connexion
  connect_file_action = new QAction(QIcon(":/images/connect.png"),tr("&Connect"),this);
  //connect_file_action->setShortcut(tr("Ctrl+O"));
  connect_file_action->setStatusTip(tr("Connect to a simulation server"));
  connect(connect_file_action, SIGNAL(triggered()), this, SLOT(actionEmpty()));
  // print
  print_file_action = new QAction(QIcon(":/images/fileprint.png"),tr("&Print"),this);
  print_file_action->setShortcut(tr("Ctrl+P"));
  print_file_action->setStatusTip(tr("Print OpenGL buffer"));
  connect(print_file_action, SIGNAL(triggered()), this, SLOT(actionEmpty()));

  // quit
  quit_file_action = new QAction(QIcon(":/images/exit.png"),tr("&Quit"),this);
  quit_file_action->setShortcut(tr("Ctrl+Q"));
  quit_file_action->setStatusTip(tr("Quit"));
  connect(quit_file_action, SIGNAL(triggered()), this, SLOT(actionQuit()));
  
  // ------- Help menu actions ---------
  // about glnemo
  about_action = new QAction(QIcon(":/images/glnemo2.png"),tr("About glnemo"),this);
  about_action->setStatusTip(tr("About glnemo2"));
  connect(about_action, SIGNAL(triggered()), form_about, SLOT(show()));
  // about Qt
  about_qt = new QAction(tr("About &Qt"),this);
  about_qt->setStatusTip(tr("Show the Qt library's About box"));
  connect(about_qt, SIGNAL(triggered()), qApp, SLOT(aboutQt()));

  // ------- options toolbar actions -------
  //
  // Full Screen
  fullscreen_action = new QAction(QIcon(":/images/window_fullscreen.png"),
                                  tr("Toggle Full Screen"),this);
  fullscreen_action->setShortcut(tr("F"));
  fullscreen_action->setStatusTip(tr("Toogle Full Screen"));
  connect( fullscreen_action, SIGNAL(activated()), this, SLOT(actionFullScreen()) );

  // Reset 
  reset_action = new QAction(QIcon(":/images/home-mdk.png"),tr("Reset to initial positions"),this);
  reset_action->setShortcut(tr("Ctrl+R"));
  reset_action->setStatusTip(tr("Reset to initial positions"));
  connect( reset_action, SIGNAL( activated() ), this, SLOT( actionReset() ) );

  // Center to COM 
  com_action = new QAction(QIcon(":/images/home-mdk.png"),tr("Center to COM"),this);
  com_action->setShortcut(tr("C"));
  com_action->setStatusTip(tr("Center to COM"));
  connect( com_action, SIGNAL( activated() ), this, SLOT( actionCenterToCom() ) );
  addAction(com_action);
  
  // render mode 
  render_mode_action = new QAction(QIcon(":/images/home-mdk.png"),tr("Rendering mode"),this);
  render_mode_action->setShortcut(tr("M"));
  render_mode_action->setStatusTip(tr("Rendering mode"));
  connect( render_mode_action, SIGNAL( activated() ), this, SLOT( actionRenderMode() ) );
  addAction(render_mode_action);

  // Fit all particles
  fit_particles_action = new QAction(QIcon(":/images/zoom-best-fit.png"),tr("Fit all particles on screen"),
                                     this);
  fit_particles_action->setShortcut(tr("Ctrl+A"));
  fit_particles_action->setStatusTip(tr("Fit all particles on screen"));
  connect( fit_particles_action, SIGNAL( activated() ), this, SLOT(actionBestZoom()) );

  // Grid 
  toggle_grid_action = new QAction(QIcon(":/images/grid.png"),tr("Toggle Grid"),this);
  toggle_grid_action->setShortcut(tr("G"));
  toggle_grid_action->setStatusTip(tr("Toggle Grid"));
  connect(toggle_grid_action, SIGNAL( activated() ),  this, SLOT(actionGrid()) );

  // Particles range & color
  particles_form_action = new QAction(QIcon(":/images/colors.png"),
                                      tr("set Particles range and color"),this);
  particles_form_action->setShortcut(tr("R"));
  particles_form_action->setStatusTip(tr("set Particles range and color"));
  connect( particles_form_action, SIGNAL( activated() ),
           this, SLOT(actionFormObjectControl() ) );

  // options
  options_form_action = new QAction(QIcon(":/images/options.png"),tr("Options dialog box"),this);
  options_form_action->setShortcut(tr("O"));
  options_form_action->setStatusTip(tr("Options dialog box"));
  connect(options_form_action, SIGNAL( activated() ), this, SLOT( actionFormOptions() ) );

  // Translation
  toggle_trans_action = new QAction(QIcon(":/images/move.png"),tr("Toggle translation"),this);
  toggle_trans_action->setShortcut(tr("T"));
  toggle_trans_action->setStatusTip(tr("Toggle translation"));
  connect(toggle_trans_action, SIGNAL( activated() ), this, SLOT(actionEmpty()) );

  // play simulation
  toggle_play_action = new QAction(QIcon(":/images/player_play.png"),tr("Play next snapshot"),this);
  toggle_play_action->setShortcut(tr("p"));
  toggle_play_action->setStatusTip(tr("Play next snapshot"));
  //connect(toggle_play_action, SIGNAL( activated() ), this, SLOT(actionPlay()) );
  connect(toggle_play_action, SIGNAL( activated() ), form_options, SLOT(on_play_pressed()));
  
  // reload
  reload_action = new QAction(QIcon(":/images/reload.png"),tr("Reload snaphot"),this);
  reload_action->setShortcut(tr("l"));
  reload_action->setStatusTip(tr("Reload snapshot"));
  connect(reload_action, SIGNAL( activated() ), this, SLOT(actionReload()) );
  
  // screenshot
  screenshot_action = new QAction(QIcon(":/images/camera.png"),tr("Take a screenshot"),this);
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
  print_file_action = new QAction(QIcon(":/images/fileprint.png"),tr("Print OpenGL window"),this);
  print_file_action->setShortcut(tr(""));
  print_file_action->setStatusTip(tr("Print OpenGL window"));
  connect(print_file_action, SIGNAL( activated() ), this, SLOT(actionPrint()) );
  
  // movie
  movie_form_action = new QAction(QIcon(":/images/video_section.png"),tr("Make a movie"),this);
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
  // inverse colormap
  reverse_cmap_action = new QAction(this);
  reverse_cmap_action->setShortcut(tr("Alt+Shift+i"));
  connect( reverse_cmap_action, SIGNAL( activated() ), gl_window, SLOT( reverseColorMap() ) );
  addAction(reverse_cmap_action);
  // constant colormap
  constant_cmap_action = new QAction(this);
  constant_cmap_action->setShortcut(tr("Alt+Shift+c"));
  connect( constant_cmap_action, SIGNAL( activated() ), colormap, SLOT( constant() ) );
  addAction(constant_cmap_action);

  // Z sorting
  zsorting_action = new QAction(this);
  zsorting_action->setShortcut(tr("Z"));
  connect( zsorting_action, SIGNAL( activated() ), this, SLOT( actionZSorting() ) );
  addAction(zsorting_action);

  // Auto rotate around X 
  rotatex_action = new QAction(this);
  rotatex_action->setShortcut(tr("Ctrl+X"));
  connect( rotatex_action, SIGNAL( activated() ), this, SLOT( actionRotateX() ) );
  addAction(rotatex_action);
  // Auto rotate around Y 
  rotatey_action = new QAction(this);
  rotatey_action->setShortcut(tr("Ctrl+Y"));
  connect( rotatey_action, SIGNAL( activated() ), this, SLOT( actionRotateY() ) );
  addAction(rotatey_action);
  // Auto rotate around Z 
  rotatez_action = new QAction(this);
  rotatez_action->setShortcut(tr("Ctrl+Z"));
  connect( rotatez_action, SIGNAL( activated() ), this, SLOT( actionRotateZ() ) );
  addAction(rotatez_action);
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
void MainWindow::interactiveSelect(std::string _select)
{
  if (current_data) {
    //crv = current_data->getSnapshotRange();
#if 1
    if (!reload) {
        
      crv = current_data->getSnapshotRange();
    }

    form_spart->update(current_data,&current_data->crv_first,_select);
    form_spart->show();
#else
  selectPart("all");
#endif
  }
}
// -----------------------------------------------------------------------------
// selectPart                                                                   
// Slots connected to "select particles dialog box" to allow to load particles  
// according to the user's selection.                                           
// -----------------------------------------------------------------------------
void MainWindow::selectPart(const std::string _select)
{
  select = _select;
  if (reload && current_data) {// reload action requested   
    current_data->close();     // close the current snapshot
    delete current_data;       // delete previous object    
    current_data = plugins->getObject(snapshot); // connect
    current_data->initLoading(store_options->vel_req, s_times);
    crv = current_data->getSnapshotRange();

  } else {
    actionReset();             // reset view if menu file open
  }
  std::cerr << "MainWindow::selectPart s_times = " << s_times << "\n";
  loadNewData(select,s_times,  // load data
	      keep_all,store_options->vel_req,true);
}
// -----------------------------------------------------------------------------
// loadNewData                                                                  
// -----------------------------------------------------------------------------
void MainWindow::loadNewData(const std::string select,
                             const std::string s_time,
                             const bool keep_all, const bool load_vel, const bool interact)
{
  ParticlesObjectVector povold;
  
  if (keep_all) {;};
  if (current_data) {
    std::cerr << "MainWindow::loadNewData s_time = " << s_time << "\n";
    current_data->initLoading(load_vel, s_time);
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
    current_data->nextFrame(user_select->getIndexesTab(),user_select->getNSel());
    qDebug("Time elapsed to load snapshot: %d s", tbench.elapsed()/1000);
    store_options->new_frame=true;
    mutex_data->unlock();
    listObjects();
    if (reload) { // backup old pov2 properties to povold
      povold.clear();
      ParticlesObject::backupVVProperties(pov2,povold,pov.size());
    }
    pov2 = pov;   // copy new pov object to pov2
    ParticlesObject::clearOrbitsVectorPOV(pov2); // clear orbits vectors if present
    if (reload) { // copy back povold properties to pov2
      ParticlesObject::backupVVProperties(povold,pov2,pov.size());
    }
    form_o_c->update( current_data->part_data, &pov2,store_options);
    updateOsd();
    tbench.restart();
    gl_window->update( current_data->part_data, &pov2,store_options);
    qDebug("Time elapsed to update GL with new date: %d s", tbench.elapsed()/1000);
    if (!reload && bestzoom) gl_window->bestZoomFit();
    statusBar()->showMessage(tr("Snapshot loaded."));

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

  bench_gup_timer = new QTimer(this);
  connect(bench_gup_timer, SIGNAL(timeout()), gl_window, SLOT(updateGL()));
  bench_nframe_timer = new QTimer(this);
  connect(bench_nframe_timer, SIGNAL(timeout()),this, SLOT(updateBenchFrame()));

}
// -----------------------------------------------------------------------------
// List Objects                                                                 
// -----------------------------------------------------------------------------
void MainWindow::listObjects()
{
  for (int i=0; i<(int)pov.size();i++) {
    std::cerr << "---------------------------------------\n";
    std::cerr << "Object #["<<i<<"]: ";
    std::cerr << "#npart = "<<pov[i].npart << " -- "
              << "first  = "<<pov[i].first << " -- "
              << "last   = "<<pov[i].last  << " -- "
              << "step   = "<<pov[i].step  << "\n";
  }
}
// -----------------------------------------------------------------------------
// parse Nemo parameters                                                        
// -----------------------------------------------------------------------------
void MainWindow::parseNemoParameters()
{
  // instantiate store_options object
  store_options = new GlobalOptions();

  //                         Initialyze NEMO engine
  snapshot                = getparam((char *) "in");
  server                  = getparam((char *) "server");
  select                  = getparam((char *) "select");
  ::string select_list;
  if ( hasvalue((char *) "select_list") )
    select_list           = getparam((char *) "select_list");
  else
    select_list = NULL;
  s_times                 = getparam((char *) "times");
  keep_all                = getbparam((char *) "keep_all");
  store_options->vel_req  = getbparam((char *) "vel");
  store_options->show_vel = getbparam((char *) "disp_vel");
  //bool range_visib        = getbparam((char *) "range_visib");
  store_options->blending = getbparam((char *) "blending");
  store_options->dbuffer  = getbparam((char *) "dbuffer");
  store_options->show_grid= getbparam((char *) "grid");
  store_options->show_osd = getbparam((char *) "osd");
  store_options->perspective=getbparam((char *) "perspective");
  store_options->orthographic = !store_options->perspective;
  play                   = getbparam((char *) "play");
  bestzoom           = getbparam((char *) "bestzoom");
  store_options->xrot     = getdparam((char *) "xrot");
  store_options->yrot     = getdparam((char *) "yrot");
  store_options->zrot     = getdparam((char *) "zrot");
  store_options->xtrans   = getdparam((char *) "xtrans");
  store_options->ytrans   = getdparam((char *) "ytrans");
  store_options->ztrans   = getdparam((char *) "ztrans");
  store_options->zoom     = getdparam((char *) "zoom");
  store_options->psize    = getiparam((char *) "psize");
  store_options->port     = getiparam((char *) "port");
  store_options->show_poly= getbparam((char *) "gas");
  store_options->texture_size       =getdparam((char *) "texture_s");
  store_options->texture_alpha_color=getiparam((char *) "texture_ac");
  store_options->duplicate_mem = getbparam((char *) "smooth_gui");
  // Animation variables
  bool anim_bench=true;
  bool anim_play= true;
  QString anim_file="";
  if ( hasvalue((char *) "anim_file")) {
    
    anim_bench = getbparam((char *) "anim_bench");
    anim_play  = getbparam((char *) "anim_play");
    anim_file  = getparam((char *) "anim_file");
  }
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
void MainWindow::actionMenuFileOpen()
{
  killPlayingEvent();       // wait the end of loading thread
  QString fileName = QFileDialog::getOpenFileName(this);
  if (!fileName.isEmpty()) {
    snapshot = fileName.toStdString();
    SnapshotInterface * new_data = plugins->getObject(snapshot);
    if (new_data)  { // valid object
      mutex_loading.lock();     // protect area                  
      delete current_data;      // free memory                   
      current_data = new_data;  // link new_data                 
//       loadNewData("all","all",  // load data
//           keep_all,store_options->vel_req,true); //
      reload=false;
      bestzoom = true;
      interactiveSelect();
      mutex_loading.unlock();   // release area                  
    }
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
  close();
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
  
//   loadNewData(select,s_times,// load data
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
    fullscreen_action->setIcon(QIcon(":/images/window_nofullscreen.png"));
    showFullScreen();
  }
  else {
    setWindowIcon(QIcon(":/images/glnemo2.png"));
    fullscreen_action->setIcon(QIcon(":/images/window_fullscreen.png"));
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
  gl_window->updateGL();
}
// -----------------------------------------------------------------------------
// actionCenterToCom()                                                   
void MainWindow::actionCenterToCom(const bool ugl)
{
  double com[3] = {0., 0., 0.};
  int np=0;
  mutex_data->lock();
  if (current_data ) {
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
  std::cerr << "grid " << store_options->show_grid << "\n";
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
// actionToggleOsd()
void MainWindow::actionToggleOsd()
{
  store_options->show_osd = !store_options->show_osd;
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
  std::string framename = store_options->base_frame_name.toStdString()+"."+stm1.str()+".jpg";
  
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
    //gl_window->resizeOsd(size.width(),size.height());
    
    QRect geom = gl_window->geometry(); // save the current geometry of the GL's window

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
    
    QImage img=((gl_window->grabFrameBufferObject()).mirrored()).rgbSwapped(); // convert FBO to img

    gl_window->resize(sizegl.width(),sizegl.height());    // revert to the previous Ogl windows's size

     // !!! activate the following line if you want to see OSDr
#if 1
    
    gl_window->setGeometry(geom.x(),geom.y(),             
                           sizegl.width(),sizegl.height());
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
          img.save(QString(name.c_str()),"jpg",95);
     }
    }
}
// -----------------------------------------------------------------------------
// actionRotateX()                                                              
void MainWindow::actionRotateX()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  if (rot)  auto_rotx_timer->start(20);
  else      auto_rotx_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateY()                                                              
void MainWindow::actionRotateY()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  if (rot)  auto_roty_timer->start(20);
  else      auto_roty_timer->stop();
}
// -----------------------------------------------------------------------------
// actionRotateZ()                                                              
void MainWindow::actionRotateZ()
{
  static bool rot=false;
  rot = !rot; // toggle rotation
  if (rot)  auto_rotz_timer->start(20);
  else      auto_rotz_timer->stop();
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
        QMessageBox::information( this,tr("Warning"),current_data->endOfDataMessage(),"Ok");
        //play_animation = false;
        emit endOfSnapshot();
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
  if (loading_thread->isValidLoading()) {
    if (  current_data->getInterfaceType() == "List of Ftm"      ||
	  current_data->getInterfaceType() == "List of Gadget 2" ||
          current_data->getInterfaceType() == "List of Nemo") { 
      //mutex_data.lock();
      //pov2=pov;
      ParticlesObject::copyVVkeepProperties(pov,pov2,user_select->getNSel()); 
      form_o_c->update( current_data->part_data, &pov2,store_options,false); // update Form
      //mutex_data.unlock();
    } else {
      //pov2=pov; // modif orbits
    }
    updateOsd();
    if (store_options->auto_com) {
       actionCenterToCom(false);
    } 
    gl_window->update( current_data->part_data, &pov2,store_options);
    if (store_options->auto_play_screenshot && !store_options->auto_gl_screenshot) {
      startAutoScreenshot();
    }
  }
  else {
    if ( current_data->isEndOfData()) {
      //std::cerr << "current_data is end of data\n";
      play_animation = false;
      play_timer->stop();
      QMessageBox::information( this,tr("Warning"),current_data->endOfDataMessage(),"Ok");
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
void MainWindow::updateOsd()
{
  gl_window->setOsd(GLObjectOsd::Nbody,*current_data->part_data->nbody,false);
  gl_window->setOsd(GLObjectOsd::Time,*current_data->part_data->timu,false);
  gl_window->setOsd(GLObjectOsd::Title,QString((current_data->getFileName()).c_str()),false);
  gl_window->setOsd(GLObjectOsd::Getdata,
		    QString((current_data->getInterfaceType()).c_str()),false);
  gl_window->setOsd(GLObjectOsd::Zoom,(const float) store_options->zoom);
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
// saveIndexList                                                                
void MainWindow::saveIndexList()
{
  std::vector <int> * list = gl_window->gl_select->getList();
  if (list->size()) {
    QString fileName = QFileDialog::getSaveFileName(this,tr("Save list of indexes"));
    if (!fileName.isEmpty()) {
      QFile file(fileName);
      if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&file);
        out << "#glnemo_index_list\n";
        for (std::vector<int>::iterator i=list->begin(); i<list->end(); i++) {
           out << (*i) << "\n";
        }
        file.close();
      }
    }
  }
}
} // glnemo namespace }
