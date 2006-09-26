// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2006                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include <qapplication.h>
#include <qstylefactory.h>
#include <iostream>
#include "images/glnemo35.xpm"
void OptionsForm::init()
{
    // Form's icon
    setIcon( QPixmap( glnemo35_xpm ) );
    
    first = true;
    clearWState( WState_Polished ); 
}
void OptionsForm::destroy()
{
}
//============================================================================
// Upload Global Options
void OptionsForm::uploadOptions()
{
    if (first) return; // do not enter while not downloaded yet
    
    // match all the OptionsForm variable with their corresponding
    // from OpenGL TAB
    store_options->show_part=show_part->isChecked();
    store_options->show_vel=show_vel->isChecked();
    store_options->psize=1.+(psize->value()*(store_options->MAX_PARTICLES_SIZE-1)/psize->maxValue());
    store_options->vel_vector_size=(vel_size->value()*(store_options->MAX_VEL_VECTOR_SIZE)/vel_size->maxValue());
    store_options->blending=blending->isChecked();
    store_options->dbuffer=dbuffer->isChecked();
    store_options->particles_alpha=alpha_slider->value();
    store_options->perspective=perspective->isChecked();
    store_options->orthographic=!perspective->isChecked();
    // from Grids TAB
    //show_grid->setChecked(store_options->show_grid);
    //store_options->mesh_length=square_size->value();
    store_options->nb_meshs=nb_meshs->value();
    store_options->xy_grid=xy_grid->isChecked();
    store_options->yz_grid=yz_grid->isChecked();
    store_options->xz_grid=xz_grid->isChecked();
    store_options->show_grid=display_grid->isChecked();
    store_options->show_cube=cube_display->isChecked();
    
    // from HUD TAB
    store_options->hud=hud->isChecked();
    store_options->hud_title=hud_title->isChecked();
    store_options->hud_time=hud_time->isChecked();
    store_options->hud_zoom=hud_zoom->isChecked();
    store_options->hud_rot=hud_rot->isChecked();
    store_options->hud_trans=hud_trans->isChecked();
    store_options->hud_data_type=hud_data_type->isChecked();
    store_options->hud_nbody=hud_nbody->isChecked();
    store_options->hud_projection=hud_projection->isChecked();
    // from experimental TAB
    store_options->show_poly=show_poly->isChecked();     
    //    Octree
    store_options->octree_enable=enable_tree->isChecked();
    store_options->octree_display=display_tree->isChecked();
    store_options->octree_level=level_tree->value();
    // launch signal to glbox
    emit updateAlphaSignal(store_options->particles_alpha);
    emit updateGL();
}


//============================================================================
// download Global options
void OptionsForm::downloadOptions(GlobalOptions * options)
{
    if (! isShown()) return; // not necessary to update otions if not shown
    // position pointer on GlobalOptions object (store_options)
    store_options = options;
 
    // match all the OptionsForm variable with their corresponding
    // from OpenGL TAB
    show_part->setChecked(store_options->show_part);
    show_vel->setChecked(store_options->show_vel);
    psize->setValue((int) (psize->maxValue()*(store_options->psize-1.)/
      (store_options->MAX_PARTICLES_SIZE-1)));
    vel_size->setValue((int) (vel_size->maxValue()*(store_options->vel_vector_size)/
      (store_options->MAX_VEL_VECTOR_SIZE)));
    blending->setChecked(store_options->blending);
    dbuffer->setChecked(store_options->dbuffer);
    alpha_slider->setValue(store_options->particles_alpha);
    perspective->setChecked(store_options->perspective);
    orthographic->setChecked(store_options->orthographic);
    // from Scene Orientation TAB
    zoom->setText(QString("%1").arg(store_options->zoom,0,'f',3));
    xrot->setText(QString("%1").arg(store_options->xrot,0,'f',3));
    yrot->setText(QString("%1").arg(store_options->yrot,0,'f',3));
    zrot->setText(QString("%1").arg(store_options->zrot,0,'f',3));
    xtrans->setText(QString("%1").arg(store_options->xtrans,0,'f',3));
    ytrans->setText(QString("%1").arg(store_options->ytrans,0,'f',3));
    ztrans->setText(QString("%1").arg(store_options->ztrans,0,'f',3));
    // from Grids TAB
    //show_grid->setChecked(store_options->show_grid);
    //changeGridSize();
    square_size->setText(QString("%1").arg(store_options->mesh_length,0,'f',3));
    nb_meshs->setValue(store_options->nb_meshs);
    xy_grid->setChecked(store_options->xy_grid);
    yz_grid->setChecked(store_options->yz_grid);
    xz_grid->setChecked(store_options->xz_grid);
    display_grid->setChecked(store_options->show_grid);
    cube_display->setChecked(store_options->show_cube);
    grid_xy_button->setPalette(QPalette(store_options->col_x_grid));
    grid_yz_button->setPalette(QPalette(store_options->col_y_grid)); 
    grid_xz_button->setPalette(QPalette(store_options->col_z_grid));
    cube_button->setPalette(QPalette(store_options->col_cube));
    // from HUD TAB
    hud->setChecked(store_options->hud);
    hud_title->setChecked(store_options->hud_title);
    hud_time->setChecked(store_options->hud_time);
    hud_zoom->setChecked(store_options->hud_zoom);
    hud_rot->setChecked(store_options->hud_rot);
    hud_trans->setChecked(store_options->hud_trans);
    hud_data_type->setChecked(store_options->hud_data_type);
    hud_nbody->setChecked(store_options->hud_nbody);
    hud_projection->setChecked(store_options->hud_projection);
    bg_col_button->setPalette(QPalette(store_options->background_color)); 
    text_col_button->setPalette(QPalette(store_options->hud_color)); 	
    // from experimental TAB
    show_poly->setChecked(store_options->show_poly);     
    texture_size->setValue((int) (texture_size->maxValue()*(store_options->texture_size)/ (store_options->MAX_TEXTURE_SIZE)));    
    texture_alpha_color->setValue(store_options->texture_alpha_color);
    //      Octree
    enable_tree->setChecked(store_options->octree_enable);
    display_tree->setChecked(store_options->octree_display);
    level_tree->setValue(store_options->octree_level);
    first = false; 
}

//============================================================================
// Change gris size
void OptionsForm::changeGridSize()
{
    bool ok1=toFloat(square_size->text(), &(store_options->mesh_length));
   // std::cerr << "Mesh length ="<< store_options->mesh_length <<"\n";
    bool ok2=toInt(nb_meshs->text(), &(store_options->nb_meshs));
   // std::cerr << "nb_meshs ="<< store_options->nb_meshs <<"\n";
    if (ok1 || ok2) {
	emit resizeGrid(store_options->mesh_length,
		         store_options->nb_meshs);
            emit updateGL(); 
    }
}

//============================================================================
// Convert String to Float value
bool OptionsForm::toFloat( QString  s, float * fvalue )
{
    bool ok;
    float data=s.toFloat(&ok);
    //std::cerr << "float data =" << data<< " str [" << s <<" ] ok="<<ok<<"\n";
    if (ok) {
 if (data != *fvalue) {
     *fvalue = data; 
 } 
 else {
     ok = false;
 } 
    }
    return ok;
}


//============================================================================
// Convert String to Int value
bool OptionsForm::toInt( QString s, int * ivalue )
{
    bool ok;
    int data=s.toInt(&ok);
    //std::cerr << "Int data =" << data<< " str [" << s <<" ] ok="<<ok<<"\n";
    if (ok) {
 if (data != *ivalue) {
     *ivalue = data;
 } 
 else {
     ok = false;
 }      
    }
    return ok;
}

//============================================================================
// commit change
void OptionsForm::commitChange()
{
    changeGridSize();
    changeTransformations();
}

//============================================================================
// manga Ok button
void OptionsForm::okChange()
{
    commitChange();
    accept();
}

//============================================================================
// update transformation coordinates
void OptionsForm::changeTransformations()
{
    bool ok1=toFloat(zoom->text(), &(store_options->zoom));
    bool ok2=toFloat(xrot->text(), &(store_options->xrot));
    bool ok3=toFloat(yrot->text(), &(store_options->yrot));
    bool ok4=toFloat(zrot->text(), &(store_options->zrot));      
    bool ok5=toFloat(xtrans->text(), &(store_options->xtrans));
    bool ok6=toFloat(ytrans->text(), &(store_options->ytrans));
    bool ok7=toFloat(ztrans->text(), &(store_options->ztrans));
     if (ok1 || ok2 || ok3 || ok4 || ok5 || ok6 || ok7) {
 emit updateGL();
    }   
}

//============================================================================
// Change XY grid plan color.
void OptionsForm::setColorGridXY()
{
  const QColor color=setColor(&(grid_xy_button->paletteBackgroundColor()));
  grid_xy_button->setPalette(QPalette(color));
  store_options->col_x_grid=color;
  emit changeColorGridXY(color);
}

//============================================================================
// Change YZ grid plan color.
void OptionsForm::setColorGridYZ()
{
  const QColor color=setColor(&(grid_yz_button->paletteBackgroundColor()));
  grid_yz_button->setPalette(QPalette(color));
  store_options->col_y_grid=color;
  emit changeColorGridYZ(color);  
}

//============================================================================
// Change XZ grid plan color.
void OptionsForm::setColorGridXZ()
{
  const QColor color=setColor(&(grid_xz_button->paletteBackgroundColor()));
  grid_xz_button->setPalette(QPalette(color));
  store_options->col_z_grid=color;
  emit changeColorGridXZ(color); 
}
//===========================================================================
void OptionsForm::setColorCube()
{
  const QColor color=setColor(&(cube_button->paletteBackgroundColor()));
  cube_button->setPalette(QPalette(color));
  store_options->col_cube=color;
  emit changeColorCube(color); 
}
//============================================================================
// set color for background display
void OptionsForm::setColorBackground()
{
  const QColor color=setColor(&(bg_col_button->paletteBackgroundColor()));
  bg_col_button->setPalette(QPalette(color));
  store_options->background_color=color;
  emit changeColorBackground();     
}
//============================================================================
// set color for HUD tex 
void OptionsForm::setColorHUD()
{
  const QColor color=setColor(&(text_col_button->paletteBackgroundColor()));
  text_col_button->setPalette(QPalette(color));
  store_options->hud_color=color;
  emit changeColorHUD(color);     
}
//============================================================================
// Select color for color dialog box
QColor OptionsForm::setColor( const QColor * _color)
{
    QColor color=QColorDialog::getColor(
             QColor(*_color),
              this, "color dialog" );
    return color;
}
//============================================================================// Toggle Grids and cube
void OptionsForm::selectGrid()
{
    //store_options->show_grid = display_grid->isChecked();
    emit toggleGrid();
}
//===========================================================================
// Toggl XY grid .
void OptionsForm::selectGridX()
{
    store_options->xy_grid=xy_grid->isChecked();
    emit toggleGridX();
}
//============================================================================
// Toggl YZ grid .
void OptionsForm::selectGridY()
{
    store_options->yz_grid=yz_grid->isChecked();    
    emit toggleGridY();
}
//============================================================================
// Toggl XZ grid .
void OptionsForm::selectGridZ()
{
    store_options->xz_grid=xz_grid->isChecked(); 
    emit toggleGridZ();
}
//===========================================================================
void OptionsForm::selectCube()
{
    store_options->show_cube=cube_display->isChecked(); 
    emit toggleCube();
}
//===========================================================================
// Change texture size (gaz like particles)
void OptionsForm::changeTextureSize(int)
{
    store_options->texture_size=0.0001+(texture_size->value()*
        (store_options->MAX_TEXTURE_SIZE)/
        texture_size->maxValue()); 
    //std::cerr<< "texture size =" << store_options->texture_size <<"\n";
    emit setTextureSize(store_options->texture_size);
}

//============================================================================
// checkProjection
void OptionsForm::checkProjection()
{
    store_options->perspective=perspective->isChecked();
    store_options->orthographic=!perspective->isChecked();    
    emit setProjection();
}

//============================================================================
// set Texture Alpha Color
void OptionsForm::setTextureAlphaColor( int )
{
    store_options->texture_alpha_color = texture_alpha_color->value();
    emit changeTextureAlphaColor(store_options->texture_alpha_color);
}

//============================================================================
// set windows stype
void OptionsForm::setStyle()
{
  static int s=0;
  QStringList styles = QStyleFactory::keys();
  s++;
  s = (s)%styles.count();
  qApp->setStyle( styles[ s] );
  //WidgetView::button1Clicked();
}

//============================================================================
// 
void OptionsForm::toggleHud( bool)
{
    bool b=false;
    // activate/desactivate widgets according b status;
    if (hud->isChecked()) {
	b=true;
    }
    hud_title->setEnabled(b);
    hud_time->setEnabled(b);  
    hud_zoom->setEnabled(b);
    hud_rot->setEnabled(b);
    hud_trans->setEnabled(b);
    hud_data_type->setEnabled(b);
    hud_nbody->setEnabled(b);
    hud_projection->setEnabled(b);
    // update Options
    uploadOptions();
    // emit signal to glbox
    emit setHudActivate();
}

//============================================================================
void OptionsForm::setupOctreeSlot()
{
    if (first) return; // do not enter while not downloaded yet
    // get octree info
    store_options->octree_enable=enable_tree->isChecked();
    store_options->octree_display=display_tree->isChecked();
    store_options->octree_level=level_tree->value();
    emit sigUpdateTree();
}


void OptionsForm::updateVelVectorFactor()
{
  store_options->MAX_VEL_VECTOR_SIZE = max_vel_vec->value(); 
  store_options->vel_vector_size=(vel_size->value()*(store_options->MAX_VEL_VECTOR_SIZE)/vel_size->maxValue());
  emit sigVelVectorFactor();
}


void OptionsForm::setVelBox( bool b)
{
    vel_box->setEnabled(b);
}
