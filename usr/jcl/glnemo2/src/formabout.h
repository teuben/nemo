// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014                                  
// e-mail:   Jean-Charles.Lambert@lam.fr                                      
// address:  Centre de donneeS Astrophysique de Marseille (CeSAM)              
//           Laboratoire d'Astrophysique de Marseille                          
//           Pôle de l'Etoile, site de Château-Gombert                         
//           38, rue Frédéric Joliot-Curie                                     
//           13388 Marseille cedex 13 France                                   
//           CNRS U.M.R 7326                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
/**
	@author Jean-Charles Lambert <jean-charles.lambert@lam.fr>
*/
#ifndef GLNEMOFORMABOUT_H
#define GLNEMOFORMABOUT_H
#include "ui_formabout.h"

namespace glnemo {


class FormAbout: public QDialog {
Q_OBJECT
public:
    FormAbout(QWidget *parent = 0);
    void setVersion(const QString & _ver) { form.version->setText(_ver);}
    ~FormAbout();
private:
    Ui::FormAbout form;
    QGraphicsScene scene;  
};

}

#endif
