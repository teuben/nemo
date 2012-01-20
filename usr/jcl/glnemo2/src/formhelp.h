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
#ifndef FORMHELP_H
#define FORMHELP_H

#include <QDialog>

namespace Ui {
    class FormHelp;
}

namespace glnemo {

class FormHelp : public QDialog
{
    Q_OBJECT

public:
    explicit FormHelp(QWidget *parent = 0);
    ~FormHelp();

private:
    Ui::FormHelp *ui;
};

}
#endif // FORMHELP_H
