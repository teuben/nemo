// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#include "hostname_selection_form.h"

#include <qvariant.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

/*
 *  Constructs a HostnameSelectionForm as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
HostnameSelectionForm::HostnameSelectionForm( QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "HostnameSelectionForm" );
    //setSizeGripEnabled( FALSE );
    //setModal( TRUE );
    setWFlags(Qt::WType_Dialog);
    QWidget* privateLayoutWidget = new QWidget( this, "layout4" );
    privateLayoutWidget->setGeometry( QRect( 10, 10, 341, 70 ) );
    layout4 = new QGridLayout( privateLayoutWidget, 1, 1, 11, 6, "layout4"); 

    layout3 = new QHBoxLayout( 0, 0, 6, "layout3"); 

    buttonOk = new QPushButton( privateLayoutWidget, "buttonOk" );
    buttonOk->setAutoDefault( TRUE );
    buttonOk->setDefault( TRUE );
    layout3->addWidget( buttonOk );
    QSpacerItem* spacer = new QSpacerItem( 161, 20, QSizePolicy::Expanding, QSizePolicy::Minimum );
    layout3->addItem( spacer );

    buttonCancel = new QPushButton( privateLayoutWidget, "buttonCancel" );
    buttonCancel->setAutoDefault( TRUE );
    layout3->addWidget( buttonCancel );

    layout4->addMultiCellLayout( layout3, 1, 1, 0, 1 );

    text_hostname = new QLabel( privateLayoutWidget, "text_hostname" );

    layout4->addWidget( text_hostname, 0, 0 );

    edit_hostname = new QLineEdit( privateLayoutWidget, "edit_hostname" );

    layout4->addWidget( edit_hostname, 0, 1 );
    languageChange();
    resize( QSize(356, 74).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );

    // signals and slots connections
    connect( buttonOk, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( buttonCancel, SIGNAL( clicked() ), this, SLOT( reject() ) );
}

/*
 *  Destroys the object and frees any allocated resources
 */
HostnameSelectionForm::~HostnameSelectionForm()
{
    // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void HostnameSelectionForm::languageChange()
{
    setCaption( tr( "Network host selection form" ) );
    buttonOk->setText( tr( "&OK" ) );
    buttonOk->setAccel( QKeySequence( QString::null ) );
    buttonCancel->setText( tr( "&Cancel" ) );
    buttonCancel->setAccel( QKeySequence( QString::null ) );
    text_hostname->setText( tr( "Enter Hostname" ) );
}

