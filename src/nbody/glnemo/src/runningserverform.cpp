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
// RunningServerForm class implementation                                      
//                                                                             
// ============================================================================
#include "runningserverform.h"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qlistbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

    char * host[] = { 
          "node1"  ,"node2"  ,"node3"   , "node4",
          "node5"  ,"node6"  ,"node7"   , "node8", 
          "paxi"   ,"ouzo"   ,"beaver"  , "batis",
          "teddy"  ,"koala"  ,"amos"    , "anixi",
          "loutre" ,"slinky" ,"pentium" , "dunk" ,
          "oniro"  ,"panda"  ,"ouki"    , "heron",
          "meli"   ,"ionio"  ,"kerkira" , "pilos",
          "pirgos" ,"raccoon","tripa"   ,
          NULL };
          
/*
 *  Constructs a RunningServerForm as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
RunningServerForm::RunningServerForm(GLObjectWindow* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "RunningServerForm" );

    cancel_button = new QPushButton( this, "cancel_button" );
    cancel_button->setGeometry( QRect( 320, 280, 111, 41 ) );

    server_title = new QLabel( this, "server_title" );
    server_title->setGeometry( QRect( 20, 10, 450, 30 ) );

    server_list = new QListBox( this, "server_list" );
    server_list->setGeometry( QRect( 10, 50, 460, 220 ) );
    server_list->setColumnMode( QListBox::FixedNumber );
    server_list->setVariableWidth( TRUE );

    connect_button = new QPushButton( this, "connect_button" );
    connect_button->setGeometry( QRect( 40, 280, 111, 41 ) );
    languageChange();
    resize( QSize(480, 332).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );

    // signals and slots connections
    connect( connect_button, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancel_button, SIGNAL( clicked() ), this, SLOT( reject() ) );
    //connect( server_list, SIGNAL( clicked(QListBoxItem*) ), this, SLOT( acceptHostSelected(QListBoxItem* ) ));
    connect( server_list, SIGNAL( selected(QListBoxItem*) ), this, SLOT( hostSelected(QListBoxItem* ) ));
    connect( server_list, SIGNAL( clicked(QListBoxItem*) ), this,SLOT(hostSelected(QListBoxItem* ) ));
    connect( server_list, SIGNAL( doubleClicked(QListBoxItem*) ), this,SLOT(acceptHostSelected(QListBoxItem* ) ));
    
    fillList();
    clearWState( WState_Polished ); 
}


/*
 *  Destroys the object and frees any allocated resources
 */
RunningServerForm::~RunningServerForm()
{
    // no need to delete child widgets, Qt does it all for us
}
int RunningServerForm::findHost(QListBoxItem * item)
{
  char * p= host[0];
  int i=0;
  while (p) {
    if (item->text().find(QString(p),0) != -1 ) {
      return i;
    }
    i++;
    p=host[i];
  }
  return -1;
}
void RunningServerForm::hostSelected(QListBoxItem * item)
{
  int i = findHost(item);
  if ( i != -1 ) {
    hostname = QString(host[i]);
  } else {
    hostname ="";
  }
}
void RunningServerForm::acceptHostSelected(QListBoxItem * item)
{
  hostSelected(item);
  accept();
}
/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void RunningServerForm::languageChange()
{
    setCaption( tr( "Running Server Box" ) );
    cancel_button->setText( tr( "Cancel" ) );
    server_title->setText( tr( "<p align=\"center\">List of running simulation server</p><font></font>" ) );
    server_list->clear();
    //server_list->insertItem( tr( "New Item" ) );
    QToolTip::add( server_list, tr( "Double click to connect" ) );
    connect_button->setText( tr( "Connect" ) );
    QToolTip::add( connect_button, tr( "connect to the highlighted server" ) );
}

void RunningServerForm::fillList()
{

  int i=0;        
  char * p = host[i];  
#if 1 
  while (p) {
    
    //
    // Instantiate a new NetworkData object
    //
    NetworkData * new_virtual_data=
    new NetworkData(p);
    //PRINT_D std::cerr << "new_virtual_data = " << new_virtual_data << "\n";
    //LINE;
    
    // check connexion
    if (!new_virtual_data->isConnected()) { // not connected ?
      //QString message="["+hsl->edit_hostname->text()+
      //                "] is not a running simulation server\n";
      //QMessageBox::information( this,"Warning",message,"Ok");      
      //delete new_virtual_data;
    } 
    else { // successfull connexion
  
      // Get NBODY
      int nbody   = new_virtual_data->getNbody();
      
      QString host = QString( "%1%2" ).arg(QString(p),-30)
                                      .arg(QString( "%1" ).arg(nbody,9),20);
      server_list->insertItem(host);
    }
    delete new_virtual_data;
    i++;
    p = host[i]; 
  }
#endif  
}
