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

#include "optionsform.h"

#include <qvariant.h>
#include <qslider.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qcheckbox.h>
#include <qlayout.h>
#include <qtooltip.h>
#include <qwhatsthis.h>

/*
 *  Constructs a OptionsForm as a child of 'parent', with the
 *  name 'name' and widget flags set to 'f'.
 *
 *  The dialog will by default be modeless, unless you set 'modal' to
 *  TRUE to construct a modal dialog.
 */
OptionsForm::OptionsForm( GLBox * glbox,QWidget* parent, const char* name, bool modal, WFlags fl )
    : QDialog( parent, name, modal, fl )
{
    if ( !name )
	setName( "OptionsForm" );

    slider_part_size = new QSlider( this, "slider_part_size" );
    slider_part_size->setGeometry( QRect( 130, 10, 141, 21 ) );
    slider_part_size->setOrientation( QSlider::Horizontal );    
    // Implement
    slider_part_size->setMinValue(0);
    slider_part_size->setMaxValue(100);
    int value=(int) ((glbox->getParticlesSize()-1.0)*100./(glbox->MAX_PARTICLES_SIZE));
    slider_part_size->setValue(value);
    connect(slider_part_size,SIGNAL(valueChanged(int)),glbox,SLOT(setParticlesSize(int)));
    
    text_particle_size = new QLabel( this, "text_particle_size" );
    text_particle_size->setGeometry( QRect( 20, 10, 92, 21 ) );

    blendingBox = new QGroupBox( this, "blendingBox" );
    blendingBox->setGeometry( QRect( 20, 40, 260, 50 ) );
    
    checkBox_depthbuffer = new QCheckBox( blendingBox, "checkBox_depthbuffer" );
    checkBox_depthbuffer->setGeometry( QRect( 120, 20, 100, 20 ) );
    checkBox_depthbuffer->setChecked(glbox->statusDepthBuffer());
    connect(checkBox_depthbuffer,SIGNAL(toggled(bool)),glbox,SLOT(toggleDepthBuffer()));
    
    checkBox_blending = new QCheckBox( blendingBox, "checkBox_blending" );
    checkBox_blending->setGeometry( QRect( 20, 20, 80, 20 ) );
    checkBox_blending->setChecked(glbox->statusBlending());
    connect(checkBox_blending,SIGNAL(toggled(bool)),glbox,SLOT(toggleBlending()));
    
    gridsBox = new QGroupBox( this, "gridsBox" );
    gridsBox->setGeometry( QRect( 20, 100, 260, 50 ) );

    checkBox_gridX = new QCheckBox( gridsBox, "checkBox_gridX" );
    checkBox_gridX->setGeometry( QRect( 20, 20, 40, 21 ) );
    checkBox_gridX->setChecked(glbox->statusGridX());
    connect(checkBox_gridX,SIGNAL(toggled(bool)),glbox,SLOT(toggleGridX()));
    
    checkBox_gridZ = new QCheckBox( gridsBox, "checkBox_gridZ" );
    checkBox_gridZ->setGeometry( QRect( 190, 20, 40, 21 ) );
    checkBox_gridZ->setChecked(glbox->statusGridZ());
    connect(checkBox_gridZ,SIGNAL(toggled(bool)),glbox,SLOT(toggleGridZ()));

    checkBox_gridY = new QCheckBox( gridsBox, "checkBox_gridY" );
    checkBox_gridY->setGeometry( QRect( 110, 20, 40, 21 ) );
    checkBox_gridY->setChecked(glbox->statusGridY());
    connect(checkBox_gridY,SIGNAL(toggled(bool)),glbox,SLOT(toggleGridY()));
    
    polyGroup = new QGroupBox( this, "polyGroup" );
    polyGroup->setGeometry( QRect( 20, 160, 260, 50 ) );
    polyCheck = new QCheckBox( polyGroup, "" );
    polyCheck->setGeometry( QRect( 20, 20, 40, 21 ) );
    polyCheck->setChecked(glbox->statusPoly());
    connect(polyCheck,SIGNAL(toggled(bool)),glbox,SLOT(togglePoly()));
    hudBox = new QGroupBox( this, "hudBox" );
    //hudBox->setGeometry( QRect( 20, 160, 260, 50 ) );
    hudBox->setGeometry( QRect( 20, 210, 260, 50 ) );
    
    hud_display = new QCheckBox( hudBox, "hud box display" );
    hud_display->setGeometry( QRect( 20, 20, 200, 21 ) );
    hud_display->setChecked(glbox->statusHUD());
    connect(hud_display,SIGNAL(toggled(bool)),glbox,SLOT(toggleHUD()));
    
    
    languageChange();
    
    // main Box size
    //resize( QSize(289, 164).expandedTo(minimumSizeHint()) );
    resize( QSize(289, 285).expandedTo(minimumSizeHint()) );
    clearWState( WState_Polished );
}

/*
 *  Destroys the object and frees any allocated resources
 */
OptionsForm::~OptionsForm()
{
    // no need to delete child widgets, Qt does it all for us
}

/*
 *  Sets the strings of the subwidgets using the current
 *  language.
 */
void OptionsForm::languageChange()
{
    setCaption( tr( "Options" ) );
    QToolTip::add( slider_part_size, tr( "Increase particle size" ) );
    text_particle_size->setText( tr( "Particles's size" ) );
    blendingBox->setTitle( tr( "Blending and DepthBuffer" ) );
    checkBox_depthbuffer->setText( tr( "Depth Buffer" ) );
    QToolTip::add( checkBox_depthbuffer, tr( "Enable depth buffer" ) );
    checkBox_blending->setText( tr( "Blending" ) );
    QToolTip::add( checkBox_blending, tr( "Enable blending" ) );
    gridsBox->setTitle( tr( "Grids" ) );
    checkBox_gridX->setText( tr( "X" ) );
    checkBox_gridZ->setText( tr( "Z" ) );
    checkBox_gridY->setText( tr( "Y" ) );
    polyGroup->setTitle(tr("Polygones"));
    hudBox->setTitle( tr( "Head Up Display"));
    hud_display->setText( tr( "HUD" ) );
}

