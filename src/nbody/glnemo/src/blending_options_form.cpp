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


#include "blending_options_form.h"

#include <qvariant.h>
#include <qlabel.h>
#include <qslider.h>
#include <qcheckbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

/* 
 *  Constructs a BlendingOptionsForm as a child of 'parent', with the 
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
BlendingOptionsForm::BlendingOptionsForm( GLBox * glbox,QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "BlendingOptionsForm" );

    text_particle_size = new QLabel( this, "text_particle_size" );
    text_particle_size->setGeometry( QRect( 20, 10, 92, 21 ) );

    slider_part_size = new QSlider( this, "slider_part_size" );
    slider_part_size->setGeometry( QRect( 130, 10, 141, 21 ) );
    slider_part_size->setOrientation( QSlider::Horizontal );
    slider_part_size->setMinValue(0);
    slider_part_size->setMaxValue(100);
    int value=(int) ((glbox->getParticlesSize()-1.0)*100./(glbox->MAX_PARTICLES_SIZE));
    slider_part_size->setValue(value);
    connect(slider_part_size,SIGNAL(valueChanged(int)),glbox,SLOT(setParticlesSize(int)));

    checkBox_blending = new QCheckBox( this, "checkBox_blending" );
    checkBox_blending->setGeometry( QRect( 20, 40, 80, 20 ) );
    checkBox_blending->setChecked(glbox->statusBlending());
    connect(checkBox_blending,SIGNAL(toggled(bool)),glbox,SLOT(toggleBlending()));

    checkBox_depthbuffer = new QCheckBox( this, "checkBox_depthbuffer" );
    checkBox_depthbuffer->setGeometry( QRect( 170, 40, 100, 20 ) );
    checkBox_depthbuffer->setChecked(glbox->statusDepthBuffer());
    connect(checkBox_depthbuffer,SIGNAL(toggled(bool)),glbox,SLOT(toggleDepthBuffer()));

    languageChange();
    resize( QSize(280, 71).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );
}

/*
 *  Destroys the object and frees any allocated resources
 */
BlendingOptionsForm::~BlendingOptionsForm()
{
    // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void BlendingOptionsForm::languageChange()
{
    setCaption( tr( "Blending Options" ) );
    text_particle_size->setText( tr( "Particles's size" ) );
    QToolTip::add( slider_part_size, tr( "Increase particle size" ) );
    checkBox_blending->setText( tr( "Blending" ) );
    QToolTip::add( checkBox_blending, tr( "Enable blending" ) );
    checkBox_depthbuffer->setText( tr( "Depth Buffer" ) );
    QToolTip::add( checkBox_depthbuffer, tr( "Enable depth buffer" ) );
}

