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

#include <string.h>
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>

#include "select_particles_form.h"

#include <qvariant.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <qmessagebox.h>


/*
 *  Constructs a SelectParticlesForm as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
SelectParticlesForm::SelectParticlesForm( const char * title, const int p_nbody,QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "SelectParticlesForm" );
    setSizeGripEnabled( TRUE );

    //
    nbody = p_nbody;
    
    text_nbody = new QLabel( this, "text_nbody" );
    text_nbody->setGeometry( QRect( 10, 26, 54, 20 ) );

    text_select_part = new QLabel( this, "text_select_part" );
    text_select_part->setGeometry( QRect( 10, 60, 210, 21 ) );

    text_nbody_value = new QLabel( this, "text_nbody_value" );
    text_nbody_value->setGeometry( QRect( 70, 30, 141, 16 ) );
    text_nbody_value->setText(QString("%1").arg(nbody));
    
    line_edit_select_part = new QLineEdit( this, "line_edit_select_part" );
    line_edit_select_part->setGeometry( QRect( 10, 80, 370, 21 ) );

    buttonOk = new QPushButton( this, "buttonOk" );
    buttonOk->setGeometry( QRect( 40, 110, 100, 30 ) );
    buttonOk->setAutoDefault( TRUE );
    buttonOk->setDefault( TRUE );

    buttonCancel = new QPushButton( this, "buttonCancel" );
    buttonCancel->setGeometry( QRect( 250, 110, 100, 30 ) );
    buttonCancel->setAutoDefault( TRUE );

    text_title = new QLabel( this, title );
    text_title->setGeometry( QRect( 1, 0, 350, 20 ) );
    text_title->setText(title);
    
    languageChange();
    resize( QSize(389, 152).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );

    // signals and slots connections
    //connect( line_edit_select_part, SIGNAL(),this, SLOT());
    connect( buttonOk, SIGNAL( clicked() ), this, SLOT(isValidParticlesRange()) );
    connect( buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );

    // tab order
    setTabOrder( buttonOk, buttonCancel );
}

/*
 *  Destroys the object and frees any allocated resources
 */
SelectParticlesForm::~SelectParticlesForm()
{
    // no need to delete child widgets, Qt does it all for us
}

bool SelectParticlesForm::isValidParticlesRange()
{
  bool status=true;
  QString select_string=line_edit_select_part->text();
  
  if (  select_string!="all") {
    int * int_array = new int[nbody];
    const char * s = select_string;
    int npart = nemoinpi(const_cast<char*> (s) , int_array, nbody);
    if (npart <=0 ) {
      status = FALSE;
      cerr << "nemoinpi = [" << select_string << "] npart = "<<npart
      <<"and nbody["<<nbody<<"\n";
      //exit(1);
    }

    delete int_array;  // useless anymore
  }
  
  if (status) {
    accept();
  } 
  else {
    QString message="Particles selection string misformated";
    QMessageBox::information( this,"Warning",message,"Ok");
    //std::cerr << "error nemo loading....\n";    
  }
  return true;
}
/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void SelectParticlesForm::languageChange()
{
    setCaption( tr( "Select Particles Form" ) );
    text_nbody->setText( tr( "<b>Nbody :</b>" ) );
    text_select_part->setText( tr( "Enter a range of particles to select" ) );
    //text_nbody_value->setText( tr( "textLabel3" ) );
    buttonOk->setText( tr( "&OK" ) );
    buttonOk->setAccel( QKeySequence( QString::null ) );
    buttonCancel->setText( tr( "&Cancel" ) );
    buttonCancel->setAccel( QKeySequence( QString::null ) );
    //text_title->setText( tr( "Title<p align=\"center\"></p>" ) );
}

