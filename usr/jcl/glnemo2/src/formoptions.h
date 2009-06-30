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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
 */
#ifndef GLNEMOFORMOPTIONS_H
#define GLNEMOFORMOPTIONS_H
#include <QTime>
#include <QTimer>
#include <iostream>
#include "ui_formoptions.h"
#include "globaloptions.h"
namespace glnemo {

class FormOptions: public QDialog {
  Q_OBJECT
  public:
    FormOptions(GlobalOptions * ,QWidget *parent = 0);
    ~FormOptions();
  public slots:
    void updateFrame(const int, const int);
    void updateParticlesSelect(const int);
  private:
    Ui::FormOptions form;
    GlobalOptions * go;
    bool start;
    QTime time;
    QTimer * limited_timer;
    static int windows_size[][2];
  private slots:
    void leaveEvent ( QEvent * event ) {
      if (event) {;}
      emit leaveEvent();
    }
    //                       
    // bench tab             
    void on_bench_button_pressed();
    void stop_bench();
    // interactive select tab
    void on_save_button_pressed()    { emit save_selected(); }
    void on_obj_button_pressed()     { emit create_obj_selected(); }
    void on_sel_zoom_radio_pressed() { emit select_and_zoom(true);}
    void on_sel_only_radio_pressed() { emit select_and_zoom(false);}
    //                     
    // camera selection tab
    void on_cam_pts_display_clicked() { 
      emit setCamDisplay(form.cam_pts_display->isChecked(),
                         form.cam_path_display->isChecked());
    }
    void on_cam_path_display_clicked() { 
      emit setCamDisplay(form.cam_pts_display->isChecked(),
                         form.cam_path_display->isChecked());
    }
    void on_spline_points_valueChanged(int v) {
      emit setSplineParam(v,form.spline_scale->value());
    }
    void on_spline_scale_valueChanged(double v) {
      emit setSplineParam(form.spline_points->value(),v);
    }
    void on_cam_play_pressed();
    //                   
    // play selection tab
    void on_play_pressed();
    void on_com_clicked() { go->auto_com = form.com->isChecked();}
    //                   
    // auto-screenshot selection tab
    void on_frame_name_pressed();
    void on_reset_pressed() { go->frame_index = 0; }
    void on_auto_save_clicked() { go->auto_play_screenshot = form.auto_save->isChecked();}
    void on_auto_save_gl_clicked() { go->auto_gl_screenshot = form.auto_save_gl->isChecked();}
    void on_screen_size_activated(int index);
    
    
  signals:
    void start_bench(const bool);
    void select_and_zoom(const bool);
    void save_selected();
    void create_obj_selected();
    void leaveEvent();
    //           
    // camera tab
    void setCamDisplay(const bool, const bool);
    void setSplineParam(const int, const double);
    void startStopPlay();
    //         
    // play tab
    void playPressed();
};

}

#endif
