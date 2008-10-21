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
//                                                                             
// implementation of the main Widget Class                                     
//                                                                             
// Draw the API, manage all the events                                         
// ============================================================================

#include <qpushbutton.h>
#include <qaction.h>
#include <qlayout.h>
#include <qframe.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qtoolbar.h>
#include <qapplication.h>
#include <qfiledialog.h>
#include <qkeycode.h>
#include <qvbox.h>
#include <qimage.h>
#include <qpixmap.h>
#include <qworkspace.h>
#include <iostream>
#include <qstatusbar.h>
#include <string.h>
#include <qmessagebox.h>
#include <qcombobox.h>
#include <qcursor.h>
#include <qlineedit.h>
#include <qlistbox.h>
#include <qcheckbox.h>
#include <qgroupbox.h>
#include <qtimer.h>
#include <qprinter.h>
#include <qpainter.h>

#include <iomanip>
#include <assert.h>
#include <unistd.h>

// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "globjwin.h"
#include "glbox.h"

#include "network_data.h"
#include "movie_thread.h"
#include "filelistdata.h"

// FORM designed with 'designer'
#include "about_form.h"
#include "hostname_select_form.h"
#include "crun_server_form.h"
#include "select_nbody_form.h"
#include "particles_select_form.h"
#include "animation_form.h"

#define NETWORK_BROWSE 1
#define RECORD_FRAME 1
#define DRAWBOX 0  // DrawBox deactivated
// Include XPM images
#include "images/glnemo35.xpm"
#include "images/arrowsb.xpm"
#include "images/target1.xpm"
#include "images/screen2.xpm"
#include "images/grid.xpm"
#include "images/options_particles_range.xpm"
#include "images/tools.xpm"
#include "images/arw01rt.xpm"
#include "images/camera5.xpm"
#include "images/comm.xpm"
#include "images/allfit.xpm"
#include "images/file_print.xpm"
#if NETWORK_BROWSE
  #include "images/broadcast.xpm"
#endif  
#include "images/reload.xpm"
#include "images/recorder.xpm"

using namespace std;
#define LINE cerr << "Line: " << __LINE__ << "\n";

#define LOCAL_DEBUG 0
#include "print_debug.h"


// ============================================================================
// Constructor                                                                 
GLObjectWindow::GLObjectWindow( QWidget* parent, const char* name)
  : QMainWindow( parent, name)
{
  // instantiate store_options object
  store_options = new GlobalOptions();
  
  //                         Initialyze NEMO engine
  //initparam(argv,defv);
  ::string in             = getparam("in");
  QString  server         = getparam("server");
  ::string select         = getparam("select");

  if ( hasvalue("select_list") )
    select_list           = getparam("select_list");
  else
    select_list = NULL;
  ::string s_time         = getparam("times");
  store_options->vel_req  = getbparam("vel");
  store_options->show_vel = getbparam("disp_vel");
  range_visib             = getbparam("range_visib");
  store_options->blending = getbparam("blending");
  store_options->dbuffer  = getbparam("dbuffer");
  store_options->show_grid= getbparam("grid");
  store_options->perspective=getbparam("perspective");
  store_options->orthographic = !store_options->perspective;
  bool play                   = getbparam("play");
  bool bestzoom           = getbparam("bestzoom");
  store_options->xrot     = getdparam("xrot");
  store_options->yrot     = getdparam("yrot");
  store_options->zrot     = getdparam("zrot");
  store_options->xtrans   = getdparam("xtrans");
  store_options->ytrans   = getdparam("ytrans");
  store_options->ztrans   = getdparam("ztrans");
  store_options->zoom     = getdparam("zoom");  
  store_options->psize    = getiparam("psize");
  store_options->port     = getiparam("port");
  store_options->show_poly= getbparam("gas");
  store_options->texture_size       =getdparam("texture_s");
  store_options->texture_alpha_color=getiparam("texture_ac");
  // Animation variables
  bool anim_bench=true;
  bool anim_play= true;
  QString anim_file="";
  if ( hasvalue("anim_file")) {
    
    anim_bench = getbparam("anim_bench");
    anim_play  = getbparam("anim_play");
    anim_file  = getparam("anim_file");
  }
  float range_ortho;
  if (store_options->orthographic) {
    range_ortho=getdparam("ortho_range");
  } 
  if (store_options->port) ; // do nothing (remove compiler warning)
  
  //                         finish NEMO

  // no bodies yet ;)
  //nbody = 0;
  is_mouse_pressed   = FALSE;
  is_key_pressed     = FALSE;
  new_virtual_object = FALSE;
  // allocate part_data structure
  //part_data = new ParticlesData();
  part_data = NULL;
  
  // Initialyse mutex
  int status = pthread_mutex_init(&mutex_timer,NULL);
  if (status != 0 ) {
    std::cerr <<"error during [pthread_mutex_init], aborted...\n";
    exit(1);
  }
  status = pthread_mutex_init(&mutex_data,NULL);
  if (status != 0 ) {
    std::cerr <<"error during [pthread_mutex_init], aborted...\n";
    exit(1);
  }
  // Initialyze Variable
  initStuff();

  // icon main window
  setIcon( QPixmap( glnemo35_xpm ) );
  
  // create menu-entry fileOpen
  QAction* fileOpenAction;
  fileOpenAction = new QAction( this, "fileOpenAction" );
  fileOpenAction->setText( tr( "Open File" ) );
  connect( fileOpenAction, SIGNAL( activated() ) , this, SLOT( selectFileOpen() ) );
  fileOpenAction->setIconSet( QIconSet( glnemo35_xpm ) );
  
  // create menu-entry connectionOpen
  QAction* connexionOpenAction;
  connexionOpenAction = new QAction( this, "connexionOpenAction" );
  connexionOpenAction->setText( tr( "Open Network Connexion" ) );
  connect( connexionOpenAction, SIGNAL( activated() ) , this, SLOT(selectConnexionOpen()));
  connexionOpenAction->setIconSet( QIconSet( comm ) );
  
  // create menu-entry printBuffer
  QAction* printBufferAction;
  printBufferAction = new QAction( this, "printBufferAction");
  printBufferAction->setText( tr( "Print OpenGL buffer"));
  connect( printBufferAction, SIGNAL( activated() ) , this, SLOT(optionsPrintBuffer()));
  printBufferAction->setIconSet( QIconSet( file_print));
  
  // Create file menu
  QPopupMenu *file = new QPopupMenu( this );
  fileOpenAction->addTo( file );
  connexionOpenAction->addTo( file );
  printBufferAction->addTo( file );
  file->insertItem( "Exit",  qApp, SLOT(quit()), CTRL+Key_Q );

  // Create a Help menu
  AboutForm * about_form = new AboutForm(this);
  QPopupMenu *help_menu = new QPopupMenu( this );
  help_menu->insertItem( "About Glnemo",  about_form, SLOT(show()), QKeySequence());//CTRL+Key_B );
  
  // Create a menu bar
  QMenuBar *m = new QMenuBar( this );
  m->setSeparator( QMenuBar::InWindowsStyle );
  m->insertItem("&File", file );           // insert File menu
  m->insertItem("&Help", help_menu );      // insert Help menu

  //                        Create butons widgets for the toolbar options
  // Full Screen
  QAction * optionsToggleFullScreenAction = 
    new QAction("Toggle Full Screen",QPixmap( screen2 ),"&FullScreen", Key_F,
		this,"FullScreen");
  connect( optionsToggleFullScreenAction, SIGNAL( activated() ),
	   this, SLOT( optionsToggleFullScreen() ) );
		
  // Reset 
  QAction * optionsResetAction = 
    new QAction("Reset to initial positions",
		QPixmap( target1 ),"&Reset", CTRL+Key_R,this,"Reset");
  connect( optionsResetAction, SIGNAL( activated() ),
	   this, SLOT( optionsReset() ) );
  // Fit all particles
  QAction * optionsFitAllPartOnScreenAction = 
    new QAction("Fit all particles on screen",
		QPixmap( allfit ),"&AllFit", CTRL+Key_A,this,
                "Fit all particles on screen");
  connect( optionsFitAllPartOnScreenAction, SIGNAL( activated() ),
          this, SLOT( optionsFitAllPartOnScreen() ) );           
  // Translation
  QAction * optionsToggleTranslationAction = 
    new QAction("Toggle translation",
		QPixmap( arrowsb ),"&Translation", Key_T,this,"Translation");
  connect( optionsToggleTranslationAction, SIGNAL( activated() ),
	   this, SLOT( optionsToggleTranslation() ) );

  // Grid 
  QAction * optionsToggleGridAction = 
    new QAction("Toggle Grid",
		QPixmap( grid ),"&Grid", Key_G,this,"Grid");
  connect( optionsToggleGridAction, SIGNAL( activated() ),
	   this, SLOT( optionsToggleGrid() ) );

  // Particles range & color
  QAction * optionsParticlesRangeAction = 
    new QAction("set Particles range and color",
		QPixmap( options_particles_range ),
		"P&rc", Key_R,this,"Particles range and color");
  connect( optionsParticlesRangeAction, SIGNAL( activated() ),
	   this, SLOT( optionsParticlesRange() ) );

  // Blending options
  QAction * optionsBlendingAction = 
    new QAction("set options",
		QPixmap( tools ),
		"P&rc", Key_O,this,"set options");
  connect( optionsBlendingAction, SIGNAL( activated() ),
	   this, SLOT( optionsBlending() ) );
  
  // Playing action
  QAction * optionsTogglePlayingAction = 
    new QAction("Toogle play",
		QPixmap( arw01rt),
		"P&rc", Key_P,this,"Play");
  connect( optionsTogglePlayingAction, SIGNAL( activated() ),
	   this, SLOT( optionsTogglePlay() ) );
  //optionsTogglePlayingAction->setToggleAction( TRUE );        
  // Camera action
  QAction * optionsCameraAction = 
    new QAction("take a screenshot",
		QPixmap( camera5 ),
		"P&rc", Key_S,this,"take a screenshot");
  QString qs("");              
  connect( optionsCameraAction, SIGNAL( activated() ),
	   this, SLOT( optionsCamera()) );
  
  // Record action
  QAction * optionsRecordAction = 
    new QAction("Record a movie",
		QPixmap( recorder ),
		"P&rc", Key_M,this,"Record a movie");
  connect( optionsRecordAction, SIGNAL( activated() ),
	   this, SLOT( optionsRecord()) );
  // Animation action
  QAction * optionsAnimationAction =     
    new QAction("Animation",
		QPixmap( recorder ),
		"P&rc", Key_M,this,"Animation");
  connect( optionsAnimationAction, SIGNAL( activated() ),
	   this, SLOT( optionsAnimation()) );        
           
  //optionsRecordAction->setToggleAction( TRUE );         
#if NETWORK_BROWSE           
  // Look for network server
  QAction * optionsLookForNetworkServerAction = 
    new QAction("Look for network simulation server",
		QPixmap( broadcast ),
		"P&rc", ALT+Key_B,this,
                "Look for network simulation server");
  connect( optionsLookForNetworkServerAction, SIGNAL( activated() ),
	   this, SLOT( optionsLookForNetworkServer() ) );
#endif  
  // reload snapshot
  QAction * optionsReloadSnapshotAction = 
    new QAction("Reload snapshot",
		QPixmap( reload ),
		"P&rc", Key_L,this,"reload snapshot");
  connect( optionsReloadSnapshotAction, SIGNAL( activated() ),
	   this, SLOT( optionsReloadSnapshot() ) );
            
  // gaz action
  QAction * optionsToggleGazParticlesAction = 
    new QAction("Toggle gaz like particles",
		QPixmap( glnemo35_xpm ),
		"P&rc", ALT+Key_G,this,"Toggle gaz like particles");
  connect( optionsToggleGazParticlesAction, SIGNAL( activated() ),
	   this, SLOT( optionsToggleGazParticles() ) );
        
  // rotate around X axe
  QAction * optionsRotateAroundXAction = 
    new QAction("RX",
		QPixmap( glnemo35_xpm ),
		"P&rc", CTRL+Key_X,this,"RX");
  connect( optionsRotateAroundXAction, SIGNAL( activated() ),
	   this, SLOT( optionsRotateAroundX() ) );
#if 0           
  // 1/2 rotate around X axe
  QAction * optionsHalfRotateAroundXAction = 
    new QAction("RX1/2",
		QPixmap( glnemo35_xpm ),
		"P&rc", CTRL+SHIFT+Key_X,this,"RX1/2");
  connect( optionsHalfRotateAroundXAction, SIGNAL( activated() ),
	   this, SLOT( optionsRotateAroundX(180.) ) );
#endif           
  // rotate around Y axe
  QAction * optionsRotateAroundYAction = 
    new QAction("RY",
		QPixmap( glnemo35_xpm ),
		"P&rc", CTRL+Key_Y,this,"RY");
  connect( optionsRotateAroundYAction, SIGNAL( activated() ),
	   this, SLOT( optionsRotateAroundY() ) );

  // rotate around Z axe
  QAction * optionsRotateAroundZAction = 
    new QAction("RZ",
		QPixmap( glnemo35_xpm ),
		"P&rc", CTRL+Key_Z,this,"RZ");
  connect( optionsRotateAroundZAction, SIGNAL( activated() ),
	   this, SLOT( optionsRotateAroundZ() ) );

  //optionsToggleGazParticlesAction->setToggleAction( TRUE );    
  // print action
  QAction * optionsPrintBufferAction =
    new QAction("Print OpenGL buffer",
                QPixmap( file_print),
                "P&rc", CTRL+Key_P,this,"Print OpenGL buffer");
  connect( optionsPrintBufferAction, SIGNAL( activated() ),
	   this, SLOT(  optionsPrintBuffer()));             
  //                     Create Tool bar options
  QToolBar *optionsTools = new QToolBar( this, "options operations" );
  optionsTools->setLabel( "Options Operations" );
  optionsTools->setOrientation(Horizontal);
  moveDockWindow(optionsTools,DockLeft);

  // add buttons to the tool bar options
  optionsToggleFullScreenAction->addTo(optionsTools);
  optionsResetAction->addTo(optionsTools);
  optionsFitAllPartOnScreenAction->addTo(optionsTools);
  optionsToggleGridAction->addTo(optionsTools);
  optionsParticlesRangeAction->addTo(optionsTools);
  optionsBlendingAction->addTo(optionsTools);
  optionsToggleTranslationAction->addTo(optionsTools);
  optionsTogglePlayingAction->addTo(optionsTools);
  optionsReloadSnapshotAction->addTo(optionsTools);
  optionsCameraAction->addTo(optionsTools);
  optionsToggleGazParticlesAction->addTo(optionsTools);
  optionsPrintBufferAction->addTo(optionsTools);
  optionsAnimationAction->addTo(optionsTools);
#if  RECORD_FRAME
//  optionsRecordAction->addTo(optionsTools);
#endif  
#if  NETWORK_BROWSE 
//  optionsLookForNetworkServerAction->addTo(optionsTools);
#endif

#if DRAWBOX
  // create a Draw box
  draw_box = new DrawBox();
  draw_box->show();
#endif
  // Create a Vertical Box to store the Display
  glframe =  new QHBox( this );
  glframe->setFrameStyle( QFrame::StyledPanel | QFrame::Sunken );
  setCentralWidget( glframe );

  
  // Create AnimationEngine object
  anim_engine = new AnimationEngine();
  
  // Create an OpenGL widget inside the vertical box
  glbox = new GLBox( glframe, "glbox",store_options,anim_engine,0);
  
  // Create an OptionForm box
  options_form = new OptionsForm(this);
  // option_form download global options
  options_form->downloadOptions(store_options);

  // divers options
  connect(options_form,SIGNAL(updateGL()),glbox,SLOT(updateGL()));      
  // size Grids
  connect(options_form,SIGNAL(resizeGrid(const float,const int)),
          glbox,SLOT(resizeGrid(const float,const int)));   
  // color Grids & cube
  connect(options_form,SIGNAL(changeColorGridXY(const QColor)),
          glbox,SLOT(changeColorGridX(const QColor )));
  connect(options_form,SIGNAL(changeColorGridYZ(const QColor)),
          glbox,SLOT(changeColorGridY(const QColor )));   
  connect(options_form,SIGNAL(changeColorGridXZ(const QColor)),
          glbox,SLOT(changeColorGridZ(const QColor )));   
  connect(options_form,SIGNAL(changeColorCube(const QColor)),
          glbox,SLOT(changeColorCube(const QColor )));   
  // velocity vector
  connect(options_form,SIGNAL(sigVelVectorFactor()),glbox,SLOT(updateVelVectorFactor()));
         
  // toggle Grids & cube
  connect(options_form,SIGNAL(toggleGrid()),this,SLOT(optionsToggleGrid()));
  connect(options_form,SIGNAL(toggleGridX()),glbox,SLOT(toggleGridX()));
  connect(options_form,SIGNAL(toggleGridY()),glbox,SLOT(toggleGridY()));
  connect(options_form,SIGNAL(toggleGridZ()),glbox,SLOT(toggleGridZ()));
  connect(options_form,SIGNAL(toggleCube()),glbox,SLOT(toggleCube()));
  // gaz like particles
  connect(options_form,SIGNAL(setTextureSize(const float)),
          glbox,SLOT(setTextureSize(const float)));
  connect(options_form,SIGNAL(changeTextureAlphaColor(const int)),
          glbox,SLOT(changeTextureAlphaColor(const int)));
  connect(options_form,SIGNAL(setProjection()),this,SLOT(setProjection()));
  // octree
  connect(options_form,SIGNAL(sigUpdateTree()),
          glbox,SLOT(treeUpdate()));
  // hud
  connect(options_form,SIGNAL(setHudActivate()),glbox,SLOT(setHudActivate()));
  connect(options_form,SIGNAL(changeColorBackground()),glbox,SLOT(updateGL()));
  connect(options_form,SIGNAL(changeColorHUD(const QColor)),
          glbox,SLOT(changeColorHUD(const QColor )));
  // alpha particles
  //connect(options_form,SIGNAL(updateAlphaSignal(int)),
  // animation engine
  connect(this,SIGNAL(allowRecord()),anim_engine->record,SLOT(allowRecord()));
  connect(this,SIGNAL(newTime(const float )),anim_engine->record,SLOT(updateFrameTime(const float )));
  
  //resizeEvent(NULL);
  
  // acquire_data_thread
  ad_thread = NULL;
  // initialyze status bar
  statusBar()->message("Ready");
  //
  // TRY TO LOAD A NEMO SNAPSHOT
  //
  psv.clear();   // clear particles range vectors
  VirtualParticlesSelect::nb_select = 0;
  if (hasvalue("in")) {
    // try Snapshotdata
    virtual_data = new SnapshotData(in,select,s_time,store_options->vel_req,&mutex_data);
    if (! virtual_data->isValidData()) {
      delete virtual_data;

      // try filelistdata
      virtual_data = new FileListData(in,select,s_time,store_options->vel_req,&mutex_data);

      if (! virtual_data->isValidData()) {
        cerr << "File [" << in << "] is neither a NEMO snapshot nor a FileList file, aborting...\n";
        delete virtual_data;
        exit(1);
      }
    }
    connect(virtual_data,SIGNAL(loadedData(const ParticlesData *,
                                        ParticlesSelectVector * )),
            glbox,SLOT(getData(const ParticlesData*, ParticlesSelectVector *)));
    connect(virtual_data,SIGNAL(infoMessage(std::string)),
            this,SLOT(infoMessage(std::string )));
    // animation_engine
    //connect(virtual_data,SIGNAL(newTime(const float )),
    //anim_engine->record,SLOT(updateFrameTime(const float )));
    if ( ! virtual_data->loadPos(&psv, store_options->vel_req)) {
      std::cerr << "error nemo loading....\n";
    }
    else {
      //setObjectVisible(range_visib, 1);
      selectListIndex();
      part_data = virtual_data->getParticlesData();
      options_form->setVelBox((bool) (part_data->vel));
      glbox->getData(part_data,&psv);
      glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
      glbox->setHud(GLHudObject::Time,*part_data->timu);
      glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
      glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());

      if (bestzoom) optionsFitAllPartOnScreen();
      statusBar()->message("Snapshot loaded.");
    }
  } 
  else {
    //
    // TRY TO CONNECT ON A SIMULATION SERVER
    //
    if (hasvalue("server")) { // simulation server selected
      //
      // Instantiate a new NetworkData object
      //
      virtual_data=  new NetworkData(server,store_options->port);
      // check connexion
      if (!virtual_data->isConnected()) { // not connected ?
        QString message="["+server+
                        "] is not a running simulation server\n";
        QMessageBox::information( this,"Warning",message,"Ok");      
        delete virtual_data;
        virtual_data=NULL;
      } 
      else { //>> successfull connexion

        // Get NBODY
        virtual_data->getNbody();
        //PRINT_D std::cerr << "I got nbody : " << nbody << "\n";
                    
        psv.clear();   // clear particles range vectors
        VirtualParticlesSelect::nb_select = 0;

        // Set particle range
        virtual_data->setSelectedRange(select);
        
        // establish Signal connexion
        connect(virtual_data,SIGNAL(loadedData(const ParticlesData *,ParticlesSelectVector * )),
                glbox,SLOT(getData(const ParticlesData *,ParticlesSelectVector *)));
        // animation_engine 
        //connect(virtual_data,SIGNAL(newTime(const float )),anim_engine->record,SLOT(updateFrameTime(const float )));
 
              
        // load positions
        if ( ! virtual_data->loadPos(&psv, store_options->vel_req)) {
            QString message="error during Gyrfalcon loading";
            QMessageBox::information( this,"Warning",message,"Ok");
            std::cerr << "error during Gyrfalcon loading....\n";
            delete virtual_data;
            virtual_data=NULL;
        } 
        else { //>> true loadPos
          selectListIndex();
          part_data = virtual_data->getParticlesData();
          options_form->setVelBox((bool) (part_data->vel));
          glbox->getData(part_data,&psv);
          glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
          glbox->setHud(GLHudObject::Time,*part_data->timu);
          glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
          glbox->setHud(GLHudObject::Getdata,
                        virtual_data->getDataType()+server);
          if (bestzoom) optionsFitAllPartOnScreen();
          statusBar()->message("Network data loaded.");
          
        } //<< true load pos
      } //<< successfull connexion
    } //<< if simulation_server...
    else { // no snapshot and no server in command line
      virtual_data = NULL;
    }
  }
  
  setProjection(range_ortho);
  // automatically play next snapshot if selected
  if (play) {
    optionsTogglePlay();
  }
  // Animation variables
  if ( anim_file != "") {
    optionsAnimation( false ); // build animation form
    anim_engine->loadFrameData(anim_file);
    if (anim_bench) {
      anim_engine->benchSlot();
    }
  }
  PRINT_D cerr << "End of globwin\n";
}
// ============================================================================
// Destructor                                                                  
GLObjectWindow::~GLObjectWindow()
{
  delete store_options;
  //delete part_data;
  delete anim_engine;
  delete virtual_data;
#if DRAWBOX  
  delete draw_box;
#endif
}
// ============================================================================
// GLObjectWindow::setObjectVisible(I                                          
//                                                                             
void GLObjectWindow::setObjectVisible(bool b, int v_type) {
  for (unsigned int i=0; i< psv.size(); i++) {
    if (psv[i].vps->v_type == v_type) {
      psv[i].vps->is_visible = b;
    }
  }
}
// ============================================================================
// GLObjectWindow::selectListIndex()                                           
// create particles select vector according to list of indexes                 
void GLObjectWindow::selectListIndex()
{
  if (select_list) {
    VirtualParticlesSelect * vps = new VirtualParticlesSelect();
    try { // try to load data from particle list
      //vps->storeParticlesList(&psv,*part_data->nbody,select_list);
      vps->storeParticlesList(&psv,1000,select_list);
    }// try
    catch (int n) {
      switch (n) {
          case -1: 
              break;
          case -2:
              break;	
          default:
              assert(1);            
      } //switch
      infoMessage(vps->error_message);
    }// catch
    delete vps;
    setObjectVisible(range_visib, 1);
  }
}     
// ============================================================================
// GLObjectWindow::callMe()                                                    
// Useless...                                                                  
void GLObjectWindow::callMe() {
  std::cerr << "called.....\n";
}
// ============================================================================
// GLObjectWindow::initStuff()                                                 
// Initialyze variables                                                        
void GLObjectWindow::initStuff()
{
  static bool first=true;
  x_mouse=0;
  y_mouse=0;
  z_mouse=0;
  tx_mouse=0;
  ty_mouse=0;
  tz_mouse=0;
  is_pressed_left_button=FALSE;
  is_pressed_right_button=FALSE;
  is_pressed_middle_button=FALSE;
  is_translation=FALSE;
  if (first) {  
    play_animation=FALSE;
    first=false;
  }
}
// ============================================================================
// GLObjectWindow::deleteFirstWidget()                                         
void GLObjectWindow::deleteFirstWidget()
{
}
// ============================================================================
// GLObjectWindow::connectToHostname()                                         
// network connect to the selected item                                        
void GLObjectWindow::connectToHostname(QListBoxItem * item,QString edit_port, const bool _vel)
{
  bool ok;
  const int port = edit_port.toInt(&ok,10);
  if (ok) {
    connectToHostname(item->text(),port, _vel);
  }
}
// ============================================================================
// GLObjectWindow::connectToHostname()                                         
// network connect to the selected host                                        
void GLObjectWindow::connectToHostname(QString  host,const int _port, const bool _vel)
{
  int port=(_port==-1?store_options->port:_port);
  store_options->vel_req = _vel;
  if (host != "" ) { // a hostname has been entered
    PRINT_D std::cerr <<" Host selected : [" << host << "]\n";          
    //
    // Instantiate a new NetworkData object
    //
    NetworkData * new_virtual_data= new NetworkData(host,port);
    PRINT_D std::cerr << "new_virtual_data = " << new_virtual_data << "\n";
    //LINE;
    
    // check connexion
    if (!new_virtual_data->isConnected()) { // not connected ?
      QString message="["+host+
                      "] is not a running simulation server\non port ["+
                      QString("%1").arg(port)+"]\n";
      QMessageBox::information( this,"Warning",message,"Ok");      
      delete new_virtual_data;
    } 
    else { // successfull connexion

      // Get NBODY
      int nbody   = new_virtual_data->getNbody();
      PRINT_D std::cerr << "I got nbody : " << nbody << "\n";
      
      //
      // Launch select particles dialog box
      //
      CSelectNbodyForm * select_part = new CSelectNbodyForm(this);
      select_part->setData("Gyrfalcon Network Data",host,nbody);
      
      if (select_part->exec()) {  // click OK
                  
        psv.clear();   // clear particles range vectors
        VirtualParticlesSelect::nb_select = 0;

        // Set particle range
        new_virtual_data->setSelectedRange(select_part->getSelectedRange());
        
        // establish Signal connexion
        connect(new_virtual_data,SIGNAL(loadedData(const ParticlesData *,
                ParticlesSelectVector * )),
                glbox,SLOT(getData(const ParticlesData *, 
                ParticlesSelectVector *)));
             // animation_engine 
       //connect(new_virtual_data,SIGNAL(newTime(const float )),anim_engine->record,SLOT(updateFrameTime(const float )));
   
        // load positions
        if ( ! new_virtual_data->loadPos(&psv,store_options->vel_req)) {
            QString message="error during Gyrfalcon loading";
            QMessageBox::information( this,"Warning",message,"Ok");
            std::cerr << "error during Gyrfalcon loading....\n";   
        } 
        else {
          pthread_mutex_lock(&mutex_timer);
          if (play_animation) {
            killTimer( play_timer );
            if (! anim_engine->record->isActivated()) {
              glbox->setHudActivate(GLHudObject::Loading, FALSE);
            }
            play_animation = FALSE;
            while ( !ad_thread->wait(500)) {
              PRINT_D std::cerr << 
              "AD_THREAD still running, waiting....\n";      
            }
          }
          if (virtual_data) {     // exist previous object?
            
            delete virtual_data;   // delete previous object
          }          
          virtual_data = new_virtual_data;// main object point to new object
          pthread_mutex_unlock(&mutex_timer);        
          
          part_data = virtual_data->getParticlesData();
          options_form->setVelBox((bool) (part_data->vel));
          glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
          glbox->setHud(GLHudObject::Time,*part_data->timu);
          glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
          glbox->setHud(GLHudObject::Getdata,
                        virtual_data->getDataType()+host);

          statusBar()->message("Network data loaded.");
        }
        optionsFitAllPartOnScreen();
      } // if (select_part->exec()) 
      else {  // cancel selected part form box
        if (new_virtual_data) {
          delete new_virtual_data;
        }
      }
    }
  }
  else {
    statusBar()->message( "Connexion aborted", 2000 );        
  }
}
// ============================================================================
// GLObjectWindow::selectConnexionOpen()                                       
// Launch select connexion box, and try to connect to the host selected        
void GLObjectWindow::selectConnexionOpen()
{
  static bool first=TRUE;
  //static HostnameSelectionForm *hsl;
  static HostnameSelectForm *hsl;
  int n; // happy Red Hat 9
  
  if (n); // remove compiler warning
  // first time, create Hostname selection box
  if (first) {
    first = FALSE;
    //hsl = new HostnameSelectionForm(this);
    hsl = new HostnameSelectForm(this);
  }
  //
  // Launch select host dialog box
  //
  if (hsl->exec()) { // if Ok pressed
    try {
      bool ok;
      const int port = (QString (hsl->edit_port->text())).toInt(&ok,10);
      if (ok) {
        connectToHostname(hsl->host_edit_list->currentText(),port, hsl->vel->isChecked());
      }   
      //connectToHostname(hsl->edit_hostname->text());
      
    } 
    catch (int n) {  // catch errors thrown by SEND and RECV
      switch(n) {
      case -1 : cerr << "Catch error during SEND\n";
        break;
      case -2 : cerr << "Catch error during RECV\n";
        break;
      case -3 : cerr << "Catch error during RECV, nbuffer=0!\n";
        break;
      default :
        assert(1);
      }      
    }
  }
}
// ============================================================================
// GLObjectWindow::optionsReloadSnapshot()                                     
// reload the current snapshot                                                 
void GLObjectWindow::optionsReloadSnapshot()
{
  if ( ! virtual_data ) { // no data loaded
    return;
  }
  pthread_mutex_lock(&mutex_timer);
  if (play_animation) {
    killTimer( play_timer );
    if (! anim_engine->record->isActivated()) {
      glbox->setHudActivate(GLHudObject::Loading, FALSE);
    }
    play_animation = FALSE;
    while ( !ad_thread->wait(500)) {
      PRINT_D std::cerr << "AD_THREAD still running, waiting....\n";      
    }
  }          
  pthread_mutex_unlock(&mutex_timer);
  int nb_select_copy=VirtualParticlesSelect::nb_select;
  ParticlesSelectVector psv2 = psv; 
  psv.clear();   // clear particles range vectors
  VirtualParticlesSelect::nb_select = 0;
  // load positions
  if ( virtual_data->reload(&psv, store_options->vel_req) <= 0) {
      QString message="Unable to reload snapshot!";
      //QMessageBox::information( this,"Warning",message,"Ok");
      

  } else {

      part_data = virtual_data->getParticlesData();
      options_form->setVelBox((bool) (part_data->vel));
      glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
      glbox->setHud(GLHudObject::Time,*part_data->timu);
      glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
      glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
      if (strcmp(virtual_data->getDataType(),"File list") ) { // do not preserve range for File List
        //!!!psv = psv2;
        // copy back color and visibility
        for (unsigned int i=0; i<psv2.size();i++) {
          if (i<psv.size()) {
            psv[i].vps->col = psv2[i].vps->col;
            psv[i].vps->is_visible = psv2[i].vps->is_visible;
          }
        }
        VirtualParticlesSelect::nb_select = nb_select_copy;
      }
      glbox->getData( part_data,&psv);   // upload the data with the right color
      statusBar()->message("Snapshot reloaded.");
  }  
}
// ============================================================================
// GLObjectWindow::selectFileOpen()                                            
// launch the select file box, and try to load the selected file               
void GLObjectWindow::selectFileOpen()
{
  bool ok = false;
  VirtualData * new_virtual_data;
  
  static QString fn;
  
  fn = QFileDialog::getOpenFileName( QString::null, QString::null,
                                            this);
  if ( !fn.isEmpty() ) {
    // try SnapshotData
    new_virtual_data=
        new SnapshotData(fn,"all","all",store_options->vel_req,&mutex_data);;
    
    if ( ! new_virtual_data->isValidData() ) {
      if (new_virtual_data) {
        delete new_virtual_data;
      }
      // try filelistdata
      new_virtual_data = new FileListData(fn,"all","all",store_options->vel_req,&mutex_data);
      if (! new_virtual_data->isValidData()) {
      QString message="file ["+fn+"] is neither a NEMO snapshot OR a FileList file";
      QMessageBox::information( this,"Warning",message,"Ok");
      if (new_virtual_data) delete new_virtual_data;
      } else ok = true;
    } else ok = true;
      
    if (ok) {
      pthread_mutex_lock(&mutex_timer);
      if (play_animation) {
        killTimer( play_timer );
        if (! anim_engine->record->isActivated()) {
          glbox->setHudActivate(GLHudObject::Loading, FALSE);
        }
        play_animation = FALSE;
        while ( !ad_thread->wait(500)) {
              PRINT_D std::cerr << "AD_THREAD still running, waiting....\n";      
        }
      }          
      if (virtual_data) {             // exist previous object?
        delete virtual_data;          // delete previous object
      }
      virtual_data = new_virtual_data;   // main object point to new object
      
      pthread_mutex_unlock(&mutex_timer);
      
      // establish Signal connexion
      connect(virtual_data,SIGNAL(loadedData(const ParticlesData *,
      ParticlesSelectVector * )),
      glbox,SLOT(getData(const ParticlesData *, ParticlesSelectVector *)));
      // animation_engine 
      //connect(virtual_data,SIGNAL(newTime(const float )),anim_engine->record,SLOT(updateFrameTime(const float )));

      psv.clear();   // clear particles range vectors
      VirtualParticlesSelect::nb_select = 0;

      // load positions
      if ( ! virtual_data->loadPos(&psv,store_options->vel_req)) {
          QString message="End of snapshot Reached !";
          QMessageBox::information( this,"Warning",message,"Ok");
          std::cerr << "error nemo loading....\n";

      } else {
          part_data = virtual_data->getParticlesData();
          options_form->setVelBox((bool) (part_data->vel));
          glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
          glbox->setHud(GLHudObject::Time,*part_data->timu);
          glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
          glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
          optionsFitAllPartOnScreen();
          statusBar()->message("Snapshot loaded.");
      }
    }
  }
  else
      statusBar()->message( "Loading aborted", 2000 );
}
// ============================================================================
// GLObjectWindow::loadNextFrame()                                             
// load next virtual_data frame                                                
void GLObjectWindow::loadNextFrame()
{
  // load positions
  if ( ! virtual_data->loadPos(&psv, store_options->vel_req)) {
      QString message="End of snapshot Reached !";
      QMessageBox::information( this,"Warning",message,"Ok");
      std::cerr << "error nemo loading....\n";

  } else {
      part_data = virtual_data->getParticlesData();
      options_form->setVelBox((bool) (part_data->vel));
      glbox->setHud(GLHudObject::Nbody,*part_data->nbody);
      glbox->setHud(GLHudObject::Time,*part_data->timu);
      glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
      glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
      //optionsFitAllPartOnScreen();
      statusBar()->message(QString("Snapshot loaded [%1]").arg(*part_data->timu));
      

  }

}
// ============================================================================
// GLObjectWindow::messageWarning()                                            
// Display a warning message                                                   
void GLObjectWindow::messageWarning(QString * message, int status)
{
  play_animation = FALSE;
  killTimer( play_timer );
  PRINT_D cerr << "IN GLObjectWindow::messageLoad\n";
  if (0) QMessageBox::information( this,"Warning",*message,"Ok");
  if (status); // do nothing (remove compiler warning)
}
// ============================================================================
// GLObjectWindow::infoMessage()                                               
// Display an info message                                                     
void GLObjectWindow::infoMessage(std::string message)
{
  std::cerr << "GOT info message ["
                    << message << "]\n";
  QMessageBox::information( this,QString(""),QString(message.c_str()));
}
// ============================================================================
// GLObjectWindow::optionsToggleTranslation()                                  
// Toggle translation                                                          
void GLObjectWindow::optionsToggleTranslation()
{
  is_translation=!is_translation;
  if (is_translation) {
    //QCursor q(SizeAllCursor);
    glbox->getPixelTranslation(&tx_mouse,&ty_mouse,&tz_mouse);
  }
}
// ============================================================================
// GLObjectWindow::optionsToggleFullScreen()                                   
// Toggle fullscreen mode                                                      
void GLObjectWindow::optionsToggleFullScreen()
{
  static bool full_screen=TRUE;

  if (full_screen) {
    showFullScreen();
  }
  else {
    setIcon( QPixmap( glnemo35_xpm ) );
    showNormal();
  }
  full_screen=!full_screen;
}
// ============================================================================
// GLObjectWindow::optionsToggleGrid()                                         
// Toggle grid view                                                            
void GLObjectWindow::optionsToggleGrid()
{
  glbox->toggleGrid();
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::optionsTogglePlay()                                         
// Toggle Animation                                                            
void GLObjectWindow::optionsTogglePlay()
{
  if (! virtual_data) {
    QString message="No Data";
    QMessageBox::information( this,"Warning",message,"Ok");
  }
  else {
    play_animation = !play_animation;
    
    if (play_animation) {
      if ( virtual_data->is_end_of_data ) { // no more data        
        play_animation = FALSE;
        QString message=virtual_data->endOfDataMessage();
        QMessageBox::information( this,"Warning",message,"Ok");
        if (! anim_engine->record->isActivated()) {
          glbox->setHudActivate(GLHudObject::Loading, FALSE);
        }
      } 
      else {
        play_timer = startTimer( 220 );  // 100 ms timer events
        //play_timer = startTimer( 40 ); 
      }
    } 
    else {
      killTimer( play_timer );
      if (! anim_engine->record->isActivated()) {
        glbox->setHudActivate(GLHudObject::Loading, FALSE);
      }          

    }
  }
}

// ============================================================================
// GLObjectWindow::timerEvent()                                                
// Each timer event (1/2 second) this routine is called, then a thread is      
// started to load the new data. Once loaded, positions are  uploaded to GLBox 
void GLObjectWindow::timerEvent( QTimerEvent *e )
{
  VirtualData * ptr_vd;
  PRINT_D std::cerr << "Timer event = " << e->timerId() << "\n";
  if ( virtual_data->is_end_of_data ) { // no more data
    killTimer( play_timer );
    play_animation = FALSE;
    QString message=virtual_data->endOfDataMessage();
    QMessageBox::information( this,"Warning",message,"Ok");   
    if (! anim_engine->record->isActivated()) {
      glbox->setHudActivate(GLHudObject::Loading, FALSE);
    }
  }
  else {

    if (! anim_engine->record->isActivated()) {
      glbox->setHudToggle(GLHudObject::Loading);
    }
    pthread_mutex_lock(&mutex_timer);
    if ( play_animation            && // playing ON
        !is_key_pressed            && // no interactive user request 
        !is_mouse_pressed          && // (mouse, keyboard)    
        e->timerId() == play_timer && // timer is ok
        ! virtual_data->is_loading_thread ) {   // no thread running
        
      virtual_data->is_loading_thread = TRUE; //
      if ( ad_thread ) { 
        delete ad_thread; // delete previous thread
      }
      
      // Start the new THREAD to acquiring data 
       
      ad_thread = new AcquireDataThread(virtual_data,&psv, store_options->vel_req);
      ptr_vd = virtual_data;
      
      
      //ad_thread->start();  
      
      // Upload the new data to GLBox                         
      // >> Protect this area with mutex in case the user load
      //    a new snapshot or establish a new connexion          
      
      if (ptr_vd == virtual_data) {
      //if (ad_thread->is_loaded) {
        emit newTime(virtual_data->getTime());
        PRINT_D std::cerr << "I am here THREAD............\n";
        glbox->setHud(GLHudObject::Nbody,virtual_data->getNbody());
        glbox->setHud(GLHudObject::Time,virtual_data->getTime());
        glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
        PRINT_D std::cerr << "New TIme = ["<<virtual_data->getTime()<<"]\n";
        // protect area
        pthread_mutex_lock(&mutex_data);
        part_data = virtual_data->getParticlesData();
        options_form->setVelBox((bool) (part_data->vel));
        virtual_data->uploadGlData(&psv);
        pthread_mutex_unlock(&mutex_data);
        // free area
        emit allowRecord();       // allow recording if activated
        if (anim_engine->record->isActivated()){
          glbox->updateGL();
#if 0          
          // options added for cosmological movies
          for (int i=0; i<4; i++) {
            optionsRotateAroundX(1,1);
            optionsRotateAroundY(1,1);
            optionsRotateAroundZ(1,1);
          }
#endif          
        } 
        //virtual_data->is_loading_thread = FALSE;
      } 
      else {
        std::cerr << "PTR_VD != VIRTUAL_DATA\n";
      }
      // protect area

      ad_thread->start();//QThread::LowestPriority);

      // free area
   } 
   pthread_mutex_unlock(&mutex_timer);
  }        
}
// ============================================================================
// GLObjectWindow::optionsCamera()                                             
// take a screenshot                                                           
void GLObjectWindow::optionsCamera()
{
  takeScreenshot("");
}
// ============================================================================
// GLObjectWindow::optionsToggleGazParticles()                                 
// toggle gaz particles effect                                                 
void GLObjectWindow::optionsToggleGazParticles()
{ 
  store_options->show_poly = !store_options->show_poly;
  options_form->downloadOptions(store_options);
  glbox->updateGL();
}
// ============================================================================
// GLObjectWindow::optionsRecord()                                             
// toggle movie record                                                         
void GLObjectWindow::optionsRecord()
{
#if 0
  static MovieThread * movie=NULL;
  static bool play=false;
  if ( ! movie ) {
    movie = new MovieThread(15,glbox);
  }
  play = !play;
  if (play) {
    movie->start(QThread::LowestPriority);
  } else {
    movie->stop();
    delete movie;
    movie=NULL;
  }
#endif
}
// ============================================================================
// GLObjectWindow::takeScreenshot()                                            
// Take a screenshot                                                           
void GLObjectWindow::takeScreenshot(QString savefilename)
{
  glbox->updateGL();
  QImage img=glbox->grabFrameBuffer();
  
  // Open a Save File Dialog Box
  if  (savefilename == "" ) {
     savefilename = QFileDialog::getSaveFileName(QString::null, QString::null,
                              this, "Save ScreenShot");
  }
  if ( !savefilename.isEmpty() )  // valid filename ?
    if  ( !img.save( savefilename, "PNG" ) ) { // save the image
      QMessageBox::warning( this, "Save failed", "Error saving file" );
    }    
}
// ============================================================================
// GLObjectWindow::optionsLookForNetworkServer()                               
// scan local network to find out running simulation server                    
void GLObjectWindow::optionsLookForNetworkServer()
{
  static bool first=TRUE;
  //static RunningServerForm * rs_form;
  static CListRunningServerForm * rs_form;
  if (first) {
    //rs_form = new RunningServerForm(this);
    rs_form = new CListRunningServerForm(this);
    first=FALSE;
  }
  if (rs_form->exec()) {  
     PRINT_D std::cerr << ">>> Host selected = " << rs_form->hostname << "\n";
     connectToHostname(rs_form->hostname, -1, false);
  } 
}
// ============================================================================
// GLObjectWindow::optionsParticlesRange()                                     
// Launch particles range and color box                                        
void GLObjectWindow::optionsParticlesRange()
{
  static bool first=TRUE;
  static ParticlesSelectForm * psf;
#if 0  
  if (first) {
    setParticlesRangeForm = new SetParticlesRangeForm(glbox,&psv,nbody,pos,this);
    first = FALSE;
  }
 setParticlesRangeForm->updateData(&psv,nbody,pos);
 setParticlesRangeForm->show(); 
#else
  if (first) {
    psf = new ParticlesSelectForm(this);
      // connection to update the data
    connect(psf,
            SIGNAL(applyData(const ParticlesData *, 
               		     ParticlesSelectVector * )),
	    this,
            SLOT(getData(const ParticlesData *, 
                     ParticlesSelectVector *)));
    first = FALSE;
  }
 psf->updateData(part_data,&psv,&mutex_data);
 psf->show(); 

#endif 
}
// ============================================================================
//
void GLObjectWindow::getData(const ParticlesData * _p_data,
                     ParticlesSelectVector  * _psv)
{
   if (_p_data) ;  // remove compiler warning
   pthread_mutex_lock(&mutex_data);
   virtual_data->uploadGlData(_psv);
   pthread_mutex_unlock(&mutex_data);
}

// ============================================================================
// GLObjectWindow::optionsFitAllPartOnScreen()                                 
// best fit all particles on screen according to the projection selected       
void GLObjectWindow::optionsFitAllPartOnScreen()
{
  if (store_options->perspective) {
    optionsFitAllPartOnScreenPersp();
  } else {
    optionsFitAllPartOnScreenOrtho(0.0);
  }
}
// ============================================================================
// GLObjectWindow::optionsFitAllPartOnScreenOrtho()                            
// best fit all particles on screen from orthographic projection               
void GLObjectWindow::optionsFitAllPartOnScreenOrtho(float range)
{
  store_options->zoomo = -1;
  glbox->setOrthoProjection(range);
}
// ============================================================================
// GLObjectWindow::optionsFitAllPartOnScreenPers()                             
// fit all the particles on the screen from perspective view                   
// For each point v(x,y,z), we have to compute its screen coordinates win(x,y) 
// according to its Projection and  ModelView matrix                           
// vp = MProj x MModelView x v                                                 
// winx = viewport[0] + (1 + vpx) * viewport[2] / 2                            
// winy = viewport[1] + (1 + vpy) * viewport[3] / 2                            
//                                                                             
// We have then to resolve the previous equation to figure out the best value  
// for zoom                                                                    
void GLObjectWindow::optionsFitAllPartOnScreenPersp()
{
  const double * mProj    = glbox->getProjMatrix();
        double * mMod     = const_cast<double*> 
                            (glbox->getModelMatrix());
  const int    * viewport = glbox->getViewPort();
#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mMod[col*4+row]
//  draw_box->setGeometry(viewport[0],viewport[1],viewport[2],viewport[3]);
//  draw_box->show();
//  QPainter paint(draw_box);
  pthread_mutex_lock(&mutex_data);
  if (psv.size()) {
    // force ZOOM to fit all particles                           
    // Zoom is located in ModelView matrix at coordinates MM(2,3)
    float best_zoom;
    if (store_options->perspective) {
      MM(2,3) = -20000000.0;
      best_zoom=1000;
    } 
    else {
      MM(2,3) =  20000000.0;
      best_zoom=-1000;
      return;
    }
    
    float mid_screenx = (viewport[2]-viewport[0])/2.;
    float mid_screeny = (viewport[3]-viewport[1])/2.;
    float coo[3];  
    // loop on all object
    for (int i=0; i< (int ) psv.size(); i++ ) {
      // loop on all visible object selected particles  
      if (psv[i].vps->is_visible ) {
        for (int j  = 0; 
                j  <  psv[i].vps->ni_index;
                j ++) {
            int jndex= psv[i].vps->index_tab[j];      
            float 
            x=virtual_data->part_data->pos[jndex*3  ]+store_options->xtrans,
            y=virtual_data->part_data->pos[jndex*3+1]+store_options->ytrans,
            z=virtual_data->part_data->pos[jndex*3+2]+store_options->ztrans,
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
	  //std::cerr << px << " " << py << " " << pz << "\n";
#if 1	  
          px /= pw;
          py /= pw;
          pz /= pw;
#endif          
	  //std::cerr << px << " " << py << " " << pz << "\n";
          // compute screen coordinates
          float winx=viewport[0] + (1 + px) * viewport[2] / 2;
          float winy=viewport[1] + (1 + py) * viewport[3] / 2;
          // paint particles
          //paint.setPen(red);
	  //paint.drawPoint((int) (winx), (int) (winy));
          //paint.drawPoint(90,10);
	  //std::cerr << winx << " " << winy << "\n";
          //std::cerr << "max winx = " << winx << " max_out_winx="<<max_out_winx<<"\n";
          //std::cerr << "max winy = " << winy << " max_out_winy="<<max_out_winy<<"\n";
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
            //std::cerr << "winx = "<<winx<<" new zoom x = " << new_zoomx << "\n";
            if (new_zoomx < best_zoom) {
              best_zoom = new_zoomx;
              coo[0] = x; coo[1] = y; coo[2] = z;
              //std::cerr << "best zoom x = " << best_zoom << "\n";
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
            //std::cerr << "new zoom y = " << new_zoomy << "\n";
            //std::cerr << "winy = "<<winy<<" new zoom y = " << new_zoomy << "\n";
            if (new_zoomy < best_zoom) {
              best_zoom = new_zoomy;
              coo[0] = x; coo[1] = y; coo[2] = z;
              //std::cerr << "best zoom y= " << best_zoom << "\n";
            }  
          }        
        }
      }
      else { // object not visible
      }
    }
    PRINT_D std::cerr << "[" << best_zoom << "] Cordinates for best zoom = " 
                      << coo[0] <<" " << coo[1] << " " << coo[2] << "\n";
    glbox->setZoom(best_zoom);
  }
  pthread_mutex_unlock(&mutex_data);
}
// ============================================================================
// GLObjectWindow::optionsReset()                                              
// reset rotation and translation to 0,0,0 coordinates                         
void GLObjectWindow::optionsReset()
{
  initStuff(); // reset Variables
  glbox->setRotation(0,0,0);
  glbox->setTranslation(0,0,0);
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::optionsBlending()                                           
// launch options box                                                          
void GLObjectWindow::optionsBlending()
{
  options_form->show();
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::optionsPrintBuffer()                                        
// print opengl buffer                                                         
void GLObjectWindow::optionsPrintBuffer()
{
  QPrinter m_printer;
  if ( m_printer.setup() ) {
    glbox->updateGL();
    QImage img=glbox->grabFrameBuffer();
    QPainter painter( &m_printer );
    painter.drawImage(0,0,img);
  }
}
// ============================================================================
// GLObjectWindow::wheelEvent()                                                
// manage zoom according to wheel event                                        
void GLObjectWindow::wheelEvent(QWheelEvent * e)
{
  glbox->setZoom(e->delta());
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::resizeEvent()                                               
// updaet OpenGL display according to resizeEvent                              
void GLObjectWindow::resizeEvent( QResizeEvent *e )
{

  if ( e ) ; // do nothing...just to remove the warning :p
  glbox->setWH(glframe->width(),glframe->height());
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::mousePressEvent()                                           
// manage rotation/translation according to mousePresssEvent                   
void GLObjectWindow::mousePressEvent( QMouseEvent *e )
{
  if ( e->button() == QMouseEvent::LeftButton ) {  // left button pressed
    is_mouse_pressed = TRUE;
    is_pressed_left_button = TRUE;
    setMouseTracking(TRUE);
    last_posx = e->x();
    last_posy = e->y();
    if ( is_translation) {
      statusBar()->message("Translating X/Y"); 
    }
    else {
      statusBar()->message("Rotating X/Y");

    }
  }
  if ( e->button() == QMouseEvent::RightButton ) { // right button pressed
    is_mouse_pressed = TRUE;
    is_pressed_right_button = TRUE;
    setMouseTracking(TRUE);
    last_posz = e->x();
    if (is_translation){
      statusBar()->message("Translating Z"); 
    }
    else {
      statusBar()->message("Rotating Z"); 
    }
  }
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::mouseReleaseEvent()                                         
// manage mouseReleaseEvent                                                    
void GLObjectWindow::mouseReleaseEvent( QMouseEvent *e )
{
  if (e) ;  // do nothing... just to remove the warning :p
  is_pressed_left_button = FALSE;
  is_pressed_right_button = FALSE;
  is_mouse_pressed = FALSE;
  setMouseTracking(FALSE);
  statusBar()->message("Ready");
  options_form->downloadOptions(store_options);
#if DRAWBOX
  draw_box->draw( glbox,part_data, &psv, store_options,"");
#endif
  //draw_box->show();
}

// ============================================================================
// GLObjectWindow::mouseMoveEvent()                                            
// manage mouseMoveEvent                                                       
void GLObjectWindow::mouseMoveEvent( QMouseEvent *e )
{
  int dx=0,dy=0,dz=0;

  if (is_pressed_left_button) {
    
    // offset displcacement
    dx = e->x()-last_posx;
    dy = e->y()-last_posy;
    // save last position
    last_posx = e->x();
    last_posy = e->y();
    
    if (is_translation) {
      // total rotation
      tx_mouse+=dx;
      ty_mouse+=dy;
    }
    else {
      // total rotation
      x_mouse+=dx;
      y_mouse+=dy;
    }
  }
  if ( is_pressed_right_button) {
    // offset displcacement
    dz = e->x()-last_posz;
    // save last position
    last_posz = e->x();
    
    if (is_translation) {
      tz_mouse-=dz; // total rotation
    }
    else {
      z_mouse-=dz;  // total translation
    }
  }
  if (is_translation) {
  glbox->setTranslation(tx_mouse,ty_mouse,tz_mouse);
  }
  else {
    glbox->setRotation(y_mouse,x_mouse,z_mouse);
  }
  options_form->downloadOptions(store_options);

}
// ============================================================================
// GLObjectWindow::keyPressEvent()                                             
// manage keyboard press events                                                
void GLObjectWindow::keyPressEvent(QKeyEvent * k)
{
  
  if (k->key() == Qt::Key_Control ) {  
    is_key_pressed = TRUE;  
    is_translation = TRUE;
    glbox->getPixelTranslation(&tx_mouse,&ty_mouse,&tz_mouse);
  }
  if (k->key() == Qt::Key_A) {
    is_key_pressed = TRUE;  
    glbox->toggleLineAliased();
  }
  if (k->key() == Qt::Key_Plus) {
    is_key_pressed = TRUE;  
    glbox->setZoom(-1);
    statusBar()->message("Zoom IN");
  }
  if (k->key() == Qt::Key_Minus) {
    is_key_pressed = TRUE;    
    glbox->setZoom(+1);
    statusBar()->message("Zoom OUT");
  }
  options_form->downloadOptions(store_options); 
}
// ============================================================================
// GLObjectWindow::keyReleaseEvent()                                           
// manage keyboard release events                                              
void GLObjectWindow::keyReleaseEvent(QKeyEvent * k)
{
  
  if (k->key() == Qt::Key_Control ) {
    is_translation = FALSE;    
  }
  is_key_pressed = FALSE;
  options_form->downloadOptions(store_options);
}
// ============================================================================
// GLObjectWindow::setProjection()                                             
// set projection                                                              
void GLObjectWindow::setProjection()
{
  setProjection(0.0);
}
// ============================================================================
// GLObjectWindow::setProjection()                                             
// set projection                                                              
void GLObjectWindow::setProjection(float range_ortho)
{
  glbox->hudProjection();
  if (store_options->perspective) {
    PRINT_D std::cerr << "Perspective projection\n"; 
  } 
  else {
    PRINT_D std::cerr << "Orthographic projection\n"; 
    glbox->setOrthoProjection(range_ortho);
  }
 glbox->updateGL();
}
// ============================================================================
//                             A N I M A T I O N                               
// ============================================================================

// ============================================================================
// GLObjectWindow::optionsAnimation()                                          
void GLObjectWindow::optionsAnimation(bool display_form)
{
  // create AnimationForm box
  static AnimationForm * anim_form;
  static bool first=true;
  
  if (first) {
    first=false;
    anim_form = new AnimationForm(this);
    // ---------- Signal from ANIMFORM --------------
    connect(anim_form,SIGNAL(start_play_signal()) , anim_engine, SLOT(playSlot()));
    connect(anim_form,SIGNAL(start_record_signal()), anim_engine, SLOT(recordSlot()));
    connect(anim_form,SIGNAL(stop_play_signal())  , anim_engine, SLOT(stopSlot()));
    connect(anim_form,SIGNAL(reset_signal()),anim_engine,SLOT(resetPrSlot()));
    // move the slider, display one frame
    connect(anim_form,SIGNAL(display_frame_index(int)),anim_engine,SLOT(displayFrameIndexSlot(int )));
    connect(anim_form,SIGNAL(pr_options(bool)),anim_engine->play,SLOT(selectOptions(bool )));
    connect(anim_form,SIGNAL(pr_options(bool)),anim_engine->render,SLOT(selectOptions(bool )));
    // send options to rendering engine
    connect(anim_form,
            SIGNAL(rendering_options(const QString&,const QString&,const int,const bool)),
            anim_engine->render,
            SLOT(initOptionsSlot(const QString&, const QString&, const int, const bool )));
     // load/save frame data from file
    connect(anim_form,
            SIGNAL(load_frame_data(const QString&)),anim_engine,
            SLOT(loadFrameData(const QString&)));
    connect(anim_form,
            SIGNAL(save_frame_data(const QString&)),anim_engine,
            SLOT(saveFrameData(const QString&)));
    // start bench
    connect(anim_form,
            SIGNAL(sig_start_bench()),anim_engine,
            SLOT(benchSlot()));
   // -------- PLAY/RECORD ------------
    connect(anim_engine->record,SIGNAL(infoRecord(int, int )),anim_form,SLOT(infoRecord(int, int )));   
    // update global options
    connect(anim_engine->play,SIGNAL(uploadToGL(GlobalOptions*, const bool )),
            glbox,SLOT(updateOptions(GlobalOptions*, const bool )));
    connect(anim_engine->play,SIGNAL(infoPlay(int, int, int )),
            anim_form,SLOT(infoPlay(int, int, int)));        
    // ask for loading snapshot
    connect(anim_engine->play,SIGNAL(loadNextFrame()),this,SLOT(loadNextFrame()));
    // notify status
    connect(anim_engine->play,SIGNAL(infoStatus(QString)),anim_form,SLOT(statusSlot(QString)));
    connect(anim_engine->record,SIGNAL(infoStatus(QString)),anim_form,SLOT(statusSlot(QString)));
    // notify animation script info
    connect(anim_engine,SIGNAL(sig_info_from_file( const int&, const QString& )),
            anim_form,SLOT(info_from_file( const int &, const QString & )));
    
    // -------- RENDER ------------
    // start render signal
    connect(anim_form,SIGNAL(start_render_signal()) , anim_engine, SLOT(renderSlot()));
    // reset render signale
    connect(anim_form,SIGNAL(reset_render_signal()) , anim_engine->render, SLOT(reset()));
    // update progress bar
    connect(anim_engine->render,SIGNAL(infoRender(int, int )),
            anim_form,SLOT(progressRender(int, int )));    
    // update global options
    connect(anim_engine->render,SIGNAL(uploadToGL(GlobalOptions*, const bool )),
            glbox,SLOT(updateOptions(GlobalOptions*, const bool )));
    // screenshot
    connect(anim_engine->render,SIGNAL(takeScreensgot(QImage& )),
            glbox,SLOT(takeScreenshot(QImage& )));
    // ask for loading snapshot
    connect(anim_engine->render,SIGNAL(loadNextFrame()),this,SLOT(loadNextFrame()));    
    // animform render button
    connect(anim_engine->render,SIGNAL(renderButtonTextSignal(QString )),anim_form,SLOT(renderButtonTextSlot( QString)));
    connect(anim_engine->render,SIGNAL(renderDrawBox(GlobalOptions *,const bool, const int,const QString ,const QString)),
	    this,SLOT(renderDrawBox(GlobalOptions *,const bool, const int,const QString ,const QString)));
  }
  if (display_form) {
    anim_form->show();
  }
}
// ============================================================================
void GLObjectWindow::renderDrawBox(GlobalOptions * go,const bool only_transform,
				   const int current_frame_render,const QString dirname,
                                   const QString framename)
{
  QString indx = QString("%1").arg(current_frame_render);
  indx = indx.rightJustify(6,'0');
  QString shot_name=QString("%1/soft%2.%3.png").arg(dirname).arg(framename).arg(indx);
  
  if (only_transform) {                // only                    
    store_options->copyTransform(*go); // copy transformation data
  } else {                             // else                    
    *store_options = *go;              // copy all                
    //treeUpdate(false);
  }
#if DRAWBOX
  draw_box->draw( glbox,part_data, &psv, store_options, shot_name);
#endif
}
// ============================================================================
void GLObjectWindow::optionsRotateAroundY(float ang, int step)
{
  for (int i=0;i<ang;i+=step) {
          store_options->yrot+=step;
          glbox->setRotation((int) store_options->xrot,(int) store_options->yrot,(int) store_options->zrot);
  }
}
// ============================================================================
void GLObjectWindow::optionsRotateAroundX(float ang, int step)
{
  step=1;
  for (int i=0;i<ang;i+=step) {
          store_options->xrot+=step;
          store_options->yrot+=step;
          store_options->zrot+=step;
          glbox->setRotation((int) store_options->xrot,(int) store_options->yrot,(int) store_options->zrot);
  }
}
// ============================================================================
void GLObjectWindow::optionsRotateAroundZ(float ang, int step)
{
  for (int i=0;i<ang;i+=step) {
          store_options->zrot+=step;
          glbox->setRotation((int) store_options->xrot,(int) store_options->yrot,(int) store_options->zrot);
  }
}
// ============================================================================
