// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "formoptions.h"
#include <QFileDialog>
#include <QDir>
namespace glnemo {
  int FormOptions::windows_size[13][2] = {
    {2560,1600},{1920,1200},{1920,1080},{1600,1200},
    {1680,1050},{1440,900},{1280,1024},{1280,800},{1280,720},
    {1024,768},{800,600},{640,480},{320,200}
  };
// ============================================================================
// Constructor                                                                 
FormOptions::FormOptions(GlobalOptions * _go, QWidget *parent):QDialog(parent)
{
  if (parent) {;}  // remove compiler warning
  form.setupUi(this);
  start=false;
  limited_timer = new QTimer(this);
  connect(limited_timer, SIGNAL(timeout()), this, SLOT(stop_bench()));
  
  // density tab
  go = _go;
  // default screen resolution for offscreen rendering set to 1280x720
  form.screen_size->setCurrentIndex(8);
  form.frame_name_text->setText(QString(go->base_frame_name));
  // activate the first TAB by default
  form.options_dialog->setCurrentIndex(0);
  form.com->setChecked(go->auto_com);
  

  QString css;
  // ---------- Grid tab
  // Color cube button
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->col_cube.red()).
        arg(go->col_cube.green()).
        arg(go->col_cube.blue());  
  form.cube_color->setStyleSheet(css);
  // Color XY button
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->col_x_grid.red()).
        arg(go->col_x_grid.green()).
        arg(go->col_x_grid.blue());  
  form.xy_grid_color->setStyleSheet(css);
  // Color YZ button
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->col_y_grid.red()).
        arg(go->col_y_grid.green()).
        arg(go->col_y_grid.blue());  
  form.yz_grid_color->setStyleSheet(css);
  // Color XZ button
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->col_z_grid.red()).
        arg(go->col_z_grid.green()).
        arg(go->col_z_grid.blue());  
  form.xz_grid_color->setStyleSheet(css);
  
  // ------------- OSD tab
  // font color
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->osd_color.red()).
        arg(go->osd_color.green()).
        arg(go->osd_color.blue());  
  form.font_color->setStyleSheet(css);
  // background color
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->background_color.red()).
        arg(go->background_color.green()).
        arg(go->background_color.blue());  
  form.background_color->setStyleSheet(css);
  
  // ------------- Colorbar tab
  // font color
  css=QString("background:rgb(%1,%2,%3)").
        arg(go->gcb_color.red()).
        arg(go->gcb_color.green()).
        arg(go->gcb_color.blue());  
  form.gcb_font_color->setStyleSheet(css);
  
  update();
}
// ============================================================================
// Destructor                                                                  
FormOptions::~FormOptions()
{
}
// ============================================================================
// update                                                                 
void FormOptions::update()
{
  // Grid tabs
  form.show_grid_checkb->setChecked(go->show_grid);
  form.xy_checkb->setChecked(go->xy_grid);
  form.xz_checkb->setChecked(go->xz_grid);
  form.yz_checkb->setChecked(go->yz_grid);
  
  form.mesh_length_spin->setValue(go->mesh_length);
  form.mesh_nb_spin->setValue(go->nb_meshs);
  form.cube_checkb->setChecked(go->show_cube);
  
  // OSD tabs
  form.show_osd_checkb->setChecked(go->show_osd);
  form.osd_datatype->setChecked(go->osd_data_type);
  form.osd_title->setChecked(go->osd_title);
  form.osd_time->setChecked(go->osd_time);
  form.osd_nbody->setChecked(go->osd_nbody);
  form.osd_trans->setChecked(go->osd_trans);
  form.osd_zoom->setChecked(go->osd_zoom);
  form.osd_rot->setChecked(go->osd_rot);
  form.spin_font_size->setValue(go->osd_font_size);
  form.title_name->setText(go->osd_title_name);
  
  // ColorBar tab
  form.gcb_enable->setChecked(go->gcb_enable);
  form.gcb_height->setValue(go->gcb_pheight*100.);
  form.gcb_width->setValue(go->gcb_pwidth*100.);
  form.gcb_log->setChecked(go->gcb_logmode);
  if (go->gcb_orientation==0) form.gcb_radio_north->setChecked(true);
  if (go->gcb_orientation==1) form.gcb_radio_est->setChecked(true);
  if (go->gcb_orientation==2) form.gcb_radio_south->setChecked(true);
  if (go->gcb_orientation==3) form.gcb_radio_west->setChecked(true);
  form.gcb_spin_digit->setValue(go->gcb_ndigits);
  form.gcb_spin_font_size->setValue(go->gcb_font_size);
  form.gcb_spin_offset->setValue(go->gcb_offset);
  
}
// ============================================================================
// updateFrame                                                                 
void FormOptions::updateFrame(const int frame, const int tot)
{
  form.fps_lcd->display(1000*frame/time.elapsed());
  form.total_lcd->display(tot);
  time.restart();
}
// ============================================================================
// updateFrame                                                                 
void FormOptions::on_bench_button_pressed()
{
  if ( ! start ) {
   if (form.limit_radio->isChecked())            // ask for limited time
     limited_timer->start((int)form.time_box->value()*1000);// limited timer start
  form.time_choice->setDisabled(true);
  form.fps_lcd->display(0);
  form.total_lcd->display(0);
  form.status_label->setText("Running");
  form.bench_button->setText("Stop");
  time.restart();
  emit start_bench(true);
  start = true;
 }
 else {
  stop_bench();
  start = false;
 }
  
}
// ============================================================================
// stop_bench                                                                  
void FormOptions::stop_bench()
{
  form.time_choice->setDisabled(false);
  form.status_label->setText("Not Running");
  form.bench_button->setText("Start");
  start=false;
  limited_timer->stop();
  emit start_bench(false);
}
// ============================================================================
// TAB INTERACTIVE SELECT OPTIONS                                              
// ============================================================================

// ============================================================================
// updateParticlesSelect                                                       
void FormOptions::updateParticlesSelect(const int nbody)
{
  form.nsel_edit->setText(QString("%1").arg(nbody));
}
// ============================================================================
// TAB CAMERA OPTIONS                                                          
// ============================================================================
void FormOptions::on_cam_play_pressed()
{
  static bool playing=false;
  if (playing) {
    form.cam_play->setText("Start");
  } 
  else {
    form.cam_play->setText("Stop");
  }
  emit startStopPlay();
  playing = !playing;
}
// ============================================================================
// TAB PLAY OPTIONS                                                            
// ============================================================================
void FormOptions::on_frame_name_pressed()
{
  QString fileName = QFileDialog::getSaveFileName(this,tr("Select Frame directory"),go->base_frame_name);
  if (!fileName.isEmpty()) {
    go->base_frame_name = fileName;
    form.frame_name_text->setText(QString(fileName));
  }
}
// ============================================================================
void FormOptions::on_play_pressed()
{
  static bool play=false;
  play = ! play;
  if (play) {
    form.play->setText("STOP");
  } else {
    form.play->setText("PLAY");
  }
  emit playPressed();  
}
// ============================================================================
void FormOptions::on_screen_size_activated(int index)
{
  go->frame_width  = windows_size[index][0];
  go->frame_height = windows_size[index][1];
}
// ============================================================================
// TAB GRIDS and CUBE OPTIONS                                                            
// ============================================================================
}
