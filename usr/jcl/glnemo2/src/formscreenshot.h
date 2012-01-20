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
/**
	@author Jean-Charles Lambert <jean-charles.lambert@oamp.fr>
 */
#ifndef GLNEMOFORMSCREENSHOT_H
#define GLNEMOFORMSCREENSHOT_H
#include "ui_formscreenshot.h"
namespace glnemo {

class FormScreenshot: public QDialog {
  Q_OBJECT
  public:
    FormScreenshot(QWidget *parent = 0);
  
    ~FormScreenshot();
  private slots:
    // screen
    void on_method_screen_clicked()   { method=0; width=-1  ; height=-1  ;}
    void on_method_standard_clicked() { method=1; }
    void on_method_custom_clicked()   { method=2; }
    // standard
    void on_stand_320_clicked()     { width=320 ; height=240 ;}
    void on_stand_640_clicked()     { width=640 ; height=480 ;}
    void on_stand_800_clicked()     { width=800 ; height=600 ;}
    void on_stand_1024_clicked()    { width=1024; height=768 ;}
    void on_stand_1152_clicked()    { width=1152; height=864 ;}
    void on_stand_1280_clicked()    { width=1280; height=1024;}
    // shot button
    void on_shot_button_clicked();
  private:
    Ui::FormScreenshot form;
    int width,height;
    int method;
  signals:
    void screenshot(const int, const int);
};

}

#endif
