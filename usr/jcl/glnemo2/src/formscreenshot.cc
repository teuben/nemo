// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include "formscreenshot.h"


namespace glnemo {

// ============================================================================
// constructor
FormScreenshot::FormScreenshot(QWidget *parent)
{
  if (parent) {;}  // remove compiler warning
  form.setupUi(this);
  form.method_standard->setChecked(true);
  form.custom_box->setDisabled(true);
  width=1024;
  height=768;
  method =1;
}
// ============================================================================
// destructor
FormScreenshot::~FormScreenshot()
{
}
// ============================================================================
void FormScreenshot::on_shot_button_clicked()
{
  close();
  if (method == 2) { // custom selection
    bool ok;
    width  = (form.custom_w->text()).toInt(&ok,10);
    height = (form.custom_h->text()).toInt(&ok,10);
  }
  emit screenshot(width,height);
}

}
// ============================================================================
