// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
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
//                                                                             
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
#include <qcursor.h>
#include <qlineedit.h>
#include <qlistbox.h>

#include <iomanip>
#include <assert.h>
// Nemo stuffs
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "globjwin.h"
#include "glbox.h"

// My class
//#include "dialogform.h"
//#include "snapshot_data.h"
#include "particles_range.h"
#include "set_particles_range_form.h"
//#include "blending_options_form.h"
#include "optionsform.h"
#include "hostname_selection_form.h"
#include "network_data.h"
#include "select_particles_form.h"
#include "movie_thread.h"

// Include XPM images
#include "images/galaxy1.xpm"
#include "images/arrowsb.xpm"
#include "images/target1.xpm"
#include "images/screen2.xpm"
#include "images/grid.xpm"
#include "images/options_particles_range.xpm"
//#include "images/options_blending.xpm"
#include "images/tools.xpm"
#include "images/arw01rt.xpm"
#include "images/camera5.xpm"
#include "images/comm.xpm"
#include "images/allfit.xpm"
#include "images/broadcast.xpm"
#include "images/reload.xpm"
#include "images/recorder.xpm"

using namespace std;
#define LINE cerr << "Line: " << __LINE__ << "\n";

#define LOCAL_DEBUG 1
#include "print_debug.h"

#define NETWORK_BROWSE 0
#define RECORD_FRAME 0
// ----------------------------------------------------------------------------
// Constructor
GLObjectWindow::GLObjectWindow( QWidget* parent, const char* name)
  : QMainWindow( parent, name)
    //: QMainWindow( 0, 0, WDestructiveClose )
    // QWidget( parent, name )
{
  //                         Initialyze NEMO engine
  //initparam(argv,defv);
  ::string in       = getparam("in");
  QString  server   = getparam("server");
  ::string select   = getparam("select");
  ::string s_time   = getparam("times");
  bool     blending = getbparam("blending");
  bool     dbuffer  = getbparam("dbuffer");
  bool     show_grid= getbparam("grid");
  float    psize    = getiparam("psize");
  int      port     = getiparam("port");

  
  if (port) ; // do nothing (remove compiler warning)
  
  PRINT_D cerr << ">> select = [" << select << "]\n";
  //                         finish NEMO

  // no bodies yet ;)
  nbody = 0;
  is_mouse_pressed = FALSE;
  is_key_pressed   = FALSE;
  new_virtual_object = FALSE;
  
  // Initialyse mutex
  int status = pthread_mutex_init(&mutex_timer,NULL);
  if (status != 0 ) {
    std::cerr <<"error during [pthread_mutex_init], aborted...\n";
    std::exit(1);
  }
  // Initialyze Variable
  initStuff();

  // icon main window
  setIcon( QPixmap( galaxy1 ) );

  // Create Dialogs
  //DialogForm * dialogform = new DialogForm(this);
  setParticlesRangeForm = NULL;
  
  // create menu fileOpen
  QAction* fileOpenAction;
  fileOpenAction = new QAction( this, "fileOpenAction" );
  fileOpenAction->setText( tr( "Open File" ) );
  connect( fileOpenAction, SIGNAL( activated() ) , this, SLOT( selectFileOpen() ) );
  fileOpenAction->setIconSet( QIconSet( galaxy1 ) );
  
   // create menu connectionOpen
  QAction* connexionOpenAction;
  connexionOpenAction = new QAction( this, "connexionOpenAction" );
  connexionOpenAction->setText( tr( "Open Network Connexion" ) );
  connect( connexionOpenAction, SIGNAL( activated() ) , this, SLOT( selectConnexionOpen() ) );
  connexionOpenAction->setIconSet( QIconSet( comm ) );
   
  // Create a menu
  QPopupMenu *file = new QPopupMenu( this );
  fileOpenAction->addTo( file );
  connexionOpenAction->addTo( file );
  file->insertItem( "Exit",  qApp, SLOT(quit()), CTRL+Key_Q );
  //file->insertItem( "Options",  dialogform, SLOT(show()), CTRL+Key_O );

  // Create a menu bar
  QMenuBar *m = new QMenuBar( this );
  m->setSeparator( QMenuBar::InWindowsStyle );
  m->insertItem("&File", file );


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
#if  RECORD_FRAME
  optionsRecordAction->addTo(optionsTools);
#endif  
#if  NETWORK_BROWSE 
  optionsLookForNetworkServerAction->addTo(optionsTools);
#endif  
  
  // Create a Vertical Box to store the Display
  glframe =  new QHBox( this );
  glframe->setFrameStyle( QFrame::StyledPanel | QFrame::Sunken );
  setCentralWidget( glframe );

  // Create an OpenGL widget inside the vertical box
  glbox = new GLBox( glframe, "glbox",0,blending,dbuffer,show_grid,psize );
  //resizeEvent(NULL);
  
  // acquire_data_thread
  ad_thread = NULL;
  // initialyze status bar
  statusBar()->message("Ready");
  //
  // TRY TO LOAD A NEMO SNAPSHOT
  //
  if (hasvalue("in")) {
    virtual_data = new SnapshotData(in,select,s_time);
    if (! virtual_data->isValidData()) {
      cerr << "File [" << in << "] is not a NEMO snapshot, aborting...\n";
      exit(1);
    }
    connect(virtual_data,SIGNAL(loadedData(const int *, const float *, 
                                        const ParticlesRangeVector * )),
            glbox,SLOT(getData(const int *, const float *, const ParticlesRangeVector *)));
            
    //connect (virtual_data,SIGNAL(messageLoad(QString* )),
    //      this,SLOT(messageLoad(QString* )));
    
    if ( ! virtual_data->loadPos(&prv)) {
      std::cerr << "error nemo loading....\n";
    } else {
      nbody   = virtual_data->getNbody();
      timu    = virtual_data->getTime();
      pos     = virtual_data->getPos();
      coo_max = virtual_data->getCooMax();
      i_max   = virtual_data->getCooIndexMax();
      glbox->setHud(GLHudObject::Nbody,nbody);
      glbox->setHud(GLHudObject::Time,timu);
      glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
      glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
      optionsFitAllPartOnScreen();
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
      virtual_data=  new NetworkData(server);
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
        nbody   = virtual_data->getNbody();
        PRINT_D std::cerr << "I got nbody : " << nbody << "\n";
                    
        prv.clear();   // clear particles range vectors
        ParticlesRange::nb_select = 0;

        // Set particle range
        virtual_data->setSelectedRange(select);
        
        // establish Signal connexion
        connect(virtual_data,SIGNAL(loadedData(const int *, const float *,
        const ParticlesRangeVector * )),
        glbox,SLOT(getData(const int *, const float *, 
                    const ParticlesRangeVector *)));
              
        // load positions
        if ( ! virtual_data->loadPos(&prv)) {
            QString message="error during Gyrfalcon loading";
            QMessageBox::information( this,"Warning",message,"Ok");
            std::cerr << "error during Gyrfalcon loading....\n";
            delete virtual_data;
            virtual_data=NULL;
        } 
        else { //>> true loadPos
                  
          nbody   = virtual_data->getNbody();
          pos     = virtual_data->getPos();
          timu    = virtual_data->getTime();  
          coo_max = virtual_data->getCooMax();
          i_max   = virtual_data->getCooIndexMax();
          glbox->setHud(GLHudObject::Nbody,nbody);
          glbox->setHud(GLHudObject::Time,timu);
          glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
          glbox->setHud(GLHudObject::Getdata,
                        virtual_data->getDataType()+server);
          optionsFitAllPartOnScreen();
          statusBar()->message("Network data loaded.");
          
        } //<< true load pos
      } //<< successfull connexion
    } //<< if simulation_server...
    else { // no snapshot and no server in command line
      virtual_data = NULL;
    }
  }
  PRINT_D cerr << "End of globwin\n";
 
  //finiparam();  // garbage collecting for nemo
  
  //glbox->printViewport();  
  //glbox->ScreenDump();
}

// ----------------------------------------------------------------------------
// Destructor
GLObjectWindow::~GLObjectWindow()
{
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::initStuff()
{
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
  play_animation=FALSE;
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::deleteFirstWidget()
{

}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::connectToHostname(QListBoxItem * item)
{
  connectToHostname(item->text());
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::connectToHostname(QString  host)
{
  if (host != "" ) { // a hostname has been entered
    cerr <<" Host selected : [" << host << "]\n";
          
    //
    // Instantiate a new NetworkData object
    //
    NetworkData * new_virtual_data=
    new NetworkData(host);
    PRINT_D std::cerr << "new_virtual_data = " << new_virtual_data << "\n";
    LINE;
    
    // check connexion
    if (!new_virtual_data->isConnected()) { // not connected ?
      QString message="["+host+
                      "] is not a running simulation server\n";
      QMessageBox::information( this,"Warning",message,"Ok");      
      delete new_virtual_data;
    } 
    else { // successfull connexion

      // Get NBODY
      nbody   = new_virtual_data->getNbody();
      PRINT_D std::cerr << "I got nbody : " << nbody << "\n";
      //
      // Launch select particles dialog box
      //
      SelectParticlesForm * select_part = 
      new SelectParticlesForm("Gyrfalcon Network Data",
      nbody,this);
      
      if (select_part->exec()) {  // click OK
                  
        prv.clear();   // clear particles range vectors
        ParticlesRange::nb_select = 0;

        // Set particle range
        new_virtual_data->setSelectedRange(select_part->getSelectedRange());
        
        // establish Signal connexion
        connect(new_virtual_data,SIGNAL(loadedData(const int *, const float *,
        const ParticlesRangeVector * )),
        glbox,SLOT(getData(const int *, const float *, 
                    const ParticlesRangeVector *)));
            
        // load positions
        if ( ! new_virtual_data->loadPos(&prv)) {
            QString message="error during Gyrfalcon loading";
            QMessageBox::information( this,"Warning",message,"Ok");
            std::cerr << "error during Gyrfalcon loading....\n";   
        } 
        else {
          pthread_mutex_lock(&mutex_timer);
          if (play_animation) {
            killTimer( play_timer );
            glbox->setHudActivate(GLHudObject::Loading, FALSE);
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
          
          nbody   = virtual_data->getNbody();
          pos     = virtual_data->getPos();
          timu    = virtual_data->getTime();  
          coo_max = virtual_data->getCooMax();
          i_max   = virtual_data->getCooIndexMax();
          glbox->setHud(GLHudObject::Nbody,nbody);
          glbox->setHud(GLHudObject::Time,timu);
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
// ----------------------------------------------------------------------------
//
void GLObjectWindow::selectConnexionOpen()
{
  static bool first=TRUE;
  static HostnameSelectionForm *hsl;
  int n; // happy Red Hat 9
  
  // first time, create Hostname selection box
  if (first) {
    first = FALSE;
    hsl = new HostnameSelectionForm(this);
  }
  //
  // Launch select host dialog box
  //
  if (hsl->exec()) { // if Ok pressed
    try {
      connectToHostname(hsl->edit_hostname->text());
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
// ----------------------------------------------------------------------------
//
void GLObjectWindow::optionsReloadSnapshot()
{
  pthread_mutex_lock(&mutex_timer);
  if (play_animation) {
    killTimer( play_timer );
    glbox->setHudActivate(GLHudObject::Loading, FALSE);
    play_animation = FALSE;
    while ( !ad_thread->wait(500)) {
      PRINT_D std::cerr << "AD_THREAD still running, waiting....\n";      
    }
  }          
  pthread_mutex_unlock(&mutex_timer);
    
  //connect (virtual_data,SIGNAL(messageLoad(QString* )),
  //this,SLOT(messageLoad(QString* )));

  prv.clear();   // clear particles range vectors
  ParticlesRange::nb_select = 0;
  // load positions
  if ( virtual_data->reload(&prv) <= 0) {
      QString message="Unable to restart snapshot!";
      QMessageBox::information( this,"Warning",message,"Ok");

  } else {
      nbody   = virtual_data->getNbody();
      pos     = virtual_data->getPos();
      timu    = virtual_data->getTime();  
      coo_max = virtual_data->getCooMax();
      i_max   = virtual_data->getCooIndexMax();
      glbox->setHud(GLHudObject::Nbody,nbody);
      glbox->setHud(GLHudObject::Time,timu);
      glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
      glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
      statusBar()->message("Snapshot reloaded.");
  }  
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::selectFileOpen()
{
  if (setParticlesRangeForm) {
    setParticlesRangeForm->close();
  }
  
    QString fn = QFileDialog::getOpenFileName( QString::null, QString::null,
                                              this);
    if ( !fn.isEmpty() ) {
      SnapshotData * new_virtual_data=new SnapshotData(fn,"all","all");;
      
        if ( ! new_virtual_data->isValidData() ) {
          QString message="file ["+fn+"] is not a NEMO snapshot";
          QMessageBox::information( this,"Warning",message,"Ok");
          if (new_virtual_data) {
            delete new_virtual_data;
          }
        } 
        else {
          pthread_mutex_lock(&mutex_timer);
          if (play_animation) {
            killTimer( play_timer );
            glbox->setHudActivate(GLHudObject::Loading, FALSE);
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
          connect(virtual_data,SIGNAL(loadedData(const int *, const float *,
          const ParticlesRangeVector * )),
          glbox,SLOT(getData(const int *, const float *, const ParticlesRangeVector *)));
      
          //connect (virtual_data,SIGNAL(messageLoad(QString* )),
          //this,SLOT(messageLoad(QString* )));
      
          prv.clear();   // clear particles range vectors
          ParticlesRange::nb_select = 0;

          // load positions
          if ( ! virtual_data->loadPos(&prv)) {
             QString message="End of snapshot Reached !";
             QMessageBox::information( this,"Warning",message,"Ok");
             std::cerr << "error nemo loading....\n";

          } else {
              nbody   = virtual_data->getNbody();
              pos     = virtual_data->getPos();
              timu    = virtual_data->getTime();  
              coo_max = virtual_data->getCooMax();
              i_max   = virtual_data->getCooIndexMax();
              glbox->setHud(GLHudObject::Nbody,nbody);
              glbox->setHud(GLHudObject::Time,timu);
              glbox->setHud(GLHudObject::Title,virtual_data->getDataName());
              glbox->setHud(GLHudObject::Getdata,virtual_data->getDataType());
              statusBar()->message("Snapshot loaded.");
          }
        }
    }
    else
        statusBar()->message( "Loading aborted", 2000 );
}
// ----------------------------------------------------------------------------
// Toggle translation
void GLObjectWindow::messageWarning(QString * message, int status)
{
  play_animation = FALSE;
  killTimer( play_timer );
  PRINT_D cerr << "IN GLObjectWindow::messageLoad\n";
  if (0) QMessageBox::information( this,"Warning",*message,"Ok");
  if (status); // do nothing (remove compiler warning)
}
// ----------------------------------------------------------------------------
// Toggle translation
void GLObjectWindow::optionsToggleTranslation()
{
  is_translation=!is_translation;
  if (is_translation) {
    //QCursor q(SizeAllCursor);
    glbox->getPixelTranslation(&tx_mouse,&ty_mouse,&tz_mouse);
  }
}
// ----------------------------------------------------------------------------
// Toggle fullscreen mode
void GLObjectWindow::optionsToggleFullScreen()
{
  static bool full_screen=TRUE;

  if (full_screen) {
    showFullScreen();
  }
  else {
    setIcon( QPixmap( galaxy1 ) );
    showNormal();
  }
  full_screen=!full_screen;

}
// ----------------------------------------------------------------------------
// Toggle grid view
void GLObjectWindow::optionsToggleGrid()
{
  glbox->toggleGrid();
}
// ----------------------------------------------------------------------------
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
        glbox->setHudActivate(GLHudObject::Loading, FALSE);
      } 
      else {
        play_timer = startTimer( 1500 );  // 1/2 second timer events
        //play_timer = startTimer( 40 ); 
      }
    } 
    else {
      killTimer( play_timer );
      glbox->setHudActivate(GLHudObject::Loading, FALSE);
    }

  }
}

// ----------------------------------------------------------------------------
// Each timer event ( 1/2 second) this routine is called,
// then a thread is started to load the new data.
// Once loaded, the data are uploaded to GLBox
void GLObjectWindow::timerEvent( QTimerEvent *e )
{
  VirtualData * ptr_vd;

  if ( virtual_data->is_end_of_data ) { // no more data
    killTimer( play_timer );
    play_animation = FALSE;
    QString message=virtual_data->endOfDataMessage();
    QMessageBox::information( this,"Warning",message,"Ok");   
    glbox->setHudActivate(GLHudObject::Loading, FALSE);
  }
  else {
#if 0  
    glbox->setHudToggle(GLHudObject::Loading);
    if (!is_key_pressed            && // no interactive user request 
        !is_mouse_pressed          && // (mouse, keyboard)    
        e->timerId() == play_timer && // timer is ok
        ! virtual_data->is_loading_thread ) {   // no thread running
        
      virtual_data->is_loading_thread = TRUE; //
      if ( ad_thread ) { 
        delete ad_thread; // delete previous thread
      }
      
      // Start the new THREAD to acquiring data 
      pthread_mutex_lock(&mutex_timer); 
      ad_thread = new AcquireDataThread(virtual_data,&prv);
      ptr_vd = virtual_data;
      pthread_mutex_unlock(&mutex_timer);
      
      //ad_thread->start();  
      
      // Upload the new data to GLBox                         
      // >> Protect this area with mutex in case the user load
      //    a new snapshot or establish a new connexion          
      pthread_mutex_lock(&mutex_timer);
      if (ptr_vd == virtual_data) {
      //if (ad_thread->is_loaded) {
        std::cerr << "I am here THREAD............\n";
        glbox->setHud(GLHudObject::Nbody,virtual_data->getNbody());
        glbox->setHud(GLHudObject::Time,virtual_data->getTime());
        PRINT_D std::cerr << "New TIme = ["<<virtual_data->getTime()<<"]\n";
        virtual_data->uploadGlData(&prv);
        //virtual_data->is_loading_thread = FALSE;
      } 
      else {
        std::cerr << "PTR_VD != VIRTUAL_DATA\n";
      }
      ad_thread->start();
      pthread_mutex_unlock(&mutex_timer);
      
      
      
    }
#else 
  #if 0  
  static int no=0;
  
  char shot[100];
  sprintf(shot,"frame.%05d.png",no);
  takeScreenshot(shot);
  std::cerr <<  ">> shot = " << shot << "\n";
  no++;
  #else  
    glbox->setHudToggle(GLHudObject::Loading);
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
       
      ad_thread = new AcquireDataThread(virtual_data,&prv);
      ptr_vd = virtual_data;
      
      
      //ad_thread->start();  
      
      // Upload the new data to GLBox                         
      // >> Protect this area with mutex in case the user load
      //    a new snapshot or establish a new connexion          
      
      if (ptr_vd == virtual_data) {
      //if (ad_thread->is_loaded) {
        std::cerr << "I am here THREAD............\n";
        glbox->setHud(GLHudObject::Nbody,virtual_data->getNbody());
        glbox->setHud(GLHudObject::Time,virtual_data->getTime());
        PRINT_D std::cerr << "New TIme = ["<<virtual_data->getTime()<<"]\n";
        virtual_data->uploadGlData(&prv);        
        //virtual_data->is_loading_thread = FALSE;
      } 
      else {
        std::cerr << "PTR_VD != VIRTUAL_DATA\n";
      }
      ad_thread->start();//QThread::LowestPriority);
      
   } 
   pthread_mutex_unlock(&mutex_timer);
   #endif
#endif    
  }   
     
}
// ----------------------------------------------------------------------------
// options Camera
void GLObjectWindow::optionsCamera()
{
  takeScreenshot("");
}
// ----------------------------------------------------------------------------
// options record
void GLObjectWindow::optionsRecord()
{
  static MovieThread * movie=NULL;
  static bool play=false;
  if ( ! movie ) {
    movie = new MovieThread(15,glbox);
  }
  play = !play;
  if (play) {
    movie->start();
  } else {
    movie->stop();
    delete movie;
    movie=NULL;
  }
}
// ----------------------------------------------------------------------------
// Take a screenshot
void GLObjectWindow::takeScreenshot(QString savefilename)
{
  glbox->updateGL();
  QImage img=glbox->grabFrameBuffer();
  
  // Open a Save File Dialog Box
if (savefilename == "" ) {
   savefilename = QFileDialog::getSaveFileName(QString::null, QString::null,
                              this, "Save ScreenShot");
}
#if 0 
 else {
  QString savefilename = "test.png";                                        
}
#endif
#if 1
  if ( !savefilename.isEmpty() )  // valid filename ?
    if ( !img.save( savefilename, "PNG" ) ) { // save the image
      QMessageBox::warning( this, "Save failed", "Error saving file" );
    }
#endif    
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::optionsLookForNetworkServer()
{
  static bool first=TRUE;
  static RunningServerForm * rs_form;
  if (first) {
    rs_form = new RunningServerForm(this);
    first=FALSE;
  }
  //rs_form->fillList();
#if 1  
  if (rs_form->exec()) {
    std::cerr << ">>> Host selected = " << rs_form->hostname << "\n";
    connectToHostname(rs_form->hostname);
  }
#endif  
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::optionsParticlesRange()
{
#if 1
  //static SetParticlesRangeForm * setParticlesRangeForm;
  static bool first=TRUE;
  
  if (first) {
    setParticlesRangeForm = new SetParticlesRangeForm(glbox,&prv,nbody,pos,this);
    first = FALSE;
  }
 // if (setParticlesRangeForm->exec()) {
 //   glbox->getData(&nbody, pos, &prv);
 // }
 setParticlesRangeForm->updateData(&prv,nbody,pos);
 setParticlesRangeForm->show();
#else
  SetParticlesRangeForm * setParticlesRangeForm =
    new SetParticlesRangeForm(&prv,nbody,this,FALSE);
  if (setParticlesRangeForm->exec()) {

    glbox->getData(&nbody, pos, &prv);
  }
  delete setParticlesRangeForm;


#endif
  
}
// ----------------------------------------------------------------------------
// optionsFitAllPartOnScreen                                                   
// fit all the particles on the screen                                         
// For each point v(x,y,z), we have to compute its screen coordinates win(x,y) 
// according to its Projection and  ModelView matrix                           
// vp = MProj x MModelView x v                                                 
// winx = viewport[0] + (1 + vpx) * viewport[2] / 2                            
// winy = viewport[1] + (1 + vpy) * viewport[3] / 2                            
//                                                                             
// We have then to resolve the previous equation to figure out the best value  
// for zoom                                                                    
void GLObjectWindow::optionsFitAllPartOnScreen()
{
  const double * mProj    = glbox->getProjMatrix();
        double * mMod     = const_cast<double*> 
                            (glbox->getModelMatrix());
  const int    * viewport = glbox->getViewPort();
#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mMod[col*4+row]

// force ZOOM to fit all particles                           
// Zoom is located in ModelView matrix at coordinates MM(2,3)
MM(2,3) = -200000.0;
//MM(0,3) = glbox->getXtrans();
//MM(1,3) = glbox->getYtrans();

//std::cerr << "Viewport :" << viewport[0] << " "
//          << viewport[1]  << " "  << viewport[2] << " " << viewport[3] << "\n";
#if 0
  std::cerr << "ModelView Matrix:\n";
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {
       std::cerr << setw(10) << MM(i,j) <<"\t";
    }
    std::cerr << "\n";
  }  
  std::cerr << "Projection Matrix:\n";
  for (int i=0;i<4;i++) {
    for (int j=0;j<4;j++) {

       std::cerr << setw(10) << MP(i,j) <<"\t";
    }
    std::cerr << "\n";
  } 
#endif
  
  //
  float best_zoom=1000;
  float mid_screenx = (viewport[2]-viewport[0])/2.;
  float mid_screeny = (viewport[3]-viewport[1])/2.;
    
  // loop on all object
  for (int i=0; i< (int ) prv.size(); i++ ) {
    // loop on all visible object selected particles  
    if (prv[i].is_visible ) {
      for (int j  = prv[i].first_part; 
              j  <= prv[i].last_part; 
              j += prv[i].step_part) {
        float 
          x=pos[j*3  ],
          y=pos[j*3+1],
          z=pos[j*3+2],
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
        // compute screen coordinates
        float winx=viewport[0] + (1 + px) * viewport[2] / 2;
        float winy=viewport[1] + (1 + py) * viewport[3] / 2;
        //std::cerr << "max winx = " << winx << " max_out_winx="<<max_out_winx<<"\n";
        //std::cerr << "max winy = " << winy << " max_out_winy="<<max_out_winy<<"\n";
        // check farest particle
        bool guess_out_zoomx=false;
        bool guess_out_zoomy=false;
        float screen_coo;
#if 0        
        // proceed left side of the screen
        if (winx <= viewport[0] && (fabs(winx-viewport[0]) > max_out_winx)) {
          max_out_winx = fabs(winx-viewport[0]);
          screen_coo=viewport[0];
          guess_out_zoomx=true;
        }        
        else {
          // proceed right side of the screen    
          if (winx >= viewport[2] && (fabs(winx-viewport[2]) > max_out_winx)) {
            max_out_winx = fabs(winx-viewport[2]);
            screen_coo=viewport[2];
            guess_out_zoomx=true;
          }
        }
#endif
#if 1
        // proceed from left to mid-side of the screen
        if (winx >= viewport[0] && winx < mid_screenx) {
          screen_coo=viewport[0];
          guess_out_zoomx=true;
        } 
        else {
          // proceed from right to mid-side of the screen
          if (winx> mid_screenx && winx < viewport[2]) {
            screen_coo=viewport[2];
            guess_out_zoomx=true;
          }
        }
#endif                   
        if (guess_out_zoomx) {
          float A=(2.*(screen_coo-viewport[0])/viewport[2])-1.;
          float new_zoomx=-Mmz + (-Ppx+A*Ppw)/(MP(0,2)-A*MP(3,2));
          //std::cerr << "winx = "<<winx<<" new zoom x = " << new_zoomx << "\n";
          if (new_zoomx < best_zoom) {
            best_zoom = new_zoomx;
            //std::cerr << "best zoom x = " << best_zoom << "\n";
          }  
        }        
#if 0        
        // proceed bottom side of the screen
        if (winy >= viewport[3] && (fabs(winy-viewport[3]) > max_out_winy)) {
          max_out_winy = fabs(winy-viewport[3]);
          screen_coo=viewport[3];
          guess_out_zoomy=true;
        }        
        else {
          // proceed top side of the screen    
          if (winy <= viewport[1] && (fabs(winy-viewport[1]) > max_out_winy)) {
            max_out_winy = fabs(winy-viewport[1]);
            screen_coo=viewport[1];
            guess_out_zoomy=true;
          }
        }
#endif     
#if 1
        // proceed from bottom to mid-side of the screen           
        if (winy < viewport[3] && winy > mid_screeny) {
          screen_coo=viewport[3];
          guess_out_zoomy=true;
        }
        else {
          // proceed from top to mid-side of the screen   
          if (winy > viewport[1] && winy < mid_screeny) {
            screen_coo=viewport[1];
            guess_out_zoomy=true;
          }
        }
#endif        
        if (guess_out_zoomy) {
          float A=(2*(screen_coo-viewport[1])/viewport[3])-1;
          float new_zoomy=-Mmz + (-Ppy+A*Ppw)/(MP(1,2)-A*MP(3,2));
          //std::cerr << "new zoom y = " << new_zoomy << "\n";
          //std::cerr << "winy = "<<winy<<" new zoom y = " << new_zoomy << "\n";
          if (new_zoomy < best_zoom) {
            best_zoom = new_zoomy;
            //std::cerr << "best zoom y= " << best_zoom << "\n";
          }  
        }        
      }
    }
    else { // object not visible
    }
  }
  glbox->setZoom(best_zoom);
}
// ----------------------------------------------------------------------------
// optionsReset:                                                               
// reset rotation and translation to 0,0,0 coordinates                         
void GLObjectWindow::optionsReset()
{
  initStuff(); // reset Variables
  glbox->setRotation(0,0,0);
  glbox->setTranslation(0,0,0);
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::optionsBlending()
{
#if 0
  static BlendingOptionsForm * blending_options;
  static bool first=TRUE;

  if (first) {
    blending_options = new BlendingOptionsForm(glbox,this);
    first = FALSE;
  }

  blending_options->show();
#else
  static OptionsForm * options;
  static bool first=TRUE;

  if (first) {
    options = new OptionsForm(glbox,this);
    first = FALSE;
  }
  options->show();
#endif


}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::wheelEvent(QWheelEvent * e)
{
  glbox->setZoom(e->delta());
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::resizeEvent( QResizeEvent *e )
{

  if ( e ) ; // do nothing...just to remove the warning :p
  glbox->setWH(glframe->width(),glframe->height());
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::mousePressEvent( QMouseEvent *e )
{
  if ( e->button() == QMouseEvent::LeftButton ) {  // left button pressed
    is_mouse_pressed = TRUE;
    is_pressed_left_button = TRUE;
    setMouseTracking(TRUE);
    last_posx = e->x();
    last_posy = e->y();
    if ( is_translation) {
      statusBar()->message("translating X/Y"); 
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
    //last_posy = e->y();
    if (is_translation){
      statusBar()->message("Translating Z"); 
    }
    else {
      statusBar()->message("Rotating Z"); 
    }
  }
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::mouseReleaseEvent( QMouseEvent *e )
{
  if (e) ;  // do nothing... just to remove the warning :p
  is_pressed_left_button = FALSE;
  is_pressed_right_button = FALSE;
  is_mouse_pressed = FALSE;
  setMouseTracking(FALSE);
  statusBar()->message("Ready");
}

// ----------------------------------------------------------------------------
//
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
    
    //cerr << "dx = " << dx << "\n";
    //cerr << "dy = " << dy << "\n";
    if (is_translation) {
      // total rotation
      tx_mouse+=dx;
      ty_mouse+=dy;
      //cerr << "tx = " << tx_mouse << "\n";
      //cerr << "ty = " << ty_mouse << "\n";
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
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::keyPressEvent(QKeyEvent * k)
{
  if (k->key() == Qt::Key_Control ) {
    is_key_pressed = TRUE;
    is_translation = TRUE;
    glbox->getPixelTranslation(&tx_mouse,&ty_mouse,&tz_mouse);
  }
  if (k->key() == Qt::Key_A) {
    glbox->toggleLineAliased();
  }
  if (k->key() == Qt::Key_Plus) {
    glbox->setZoom(-1);
    statusBar()->message("Zoom IN");
  }
  if (k->key() == Qt::Key_Minus) {
    glbox->setZoom(+1);
    statusBar()->message("Zoom OUT");
  } 
}
// ----------------------------------------------------------------------------
//
void GLObjectWindow::keyReleaseEvent(QKeyEvent * k)
{
  if (k->key() == Qt::Key_Control ) {
    is_key_pressed = FALSE;
    is_translation = FALSE;    
    //cerr << "control released\n";
  }
}
// ----------------------------------------------------------------------------
