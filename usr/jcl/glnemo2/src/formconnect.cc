// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
//           Yannick Dalbin
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
#include <QMessageBox>
#include "formconnect.h"
#include "ui_formconnect.h"

namespace glnemo {

	FormConnect::FormConnect(QWidget *parent) {
          if (parent) {;} // remove compiler warning
		ui.setupUi(this);
	}

	FormConnect::~FormConnect() {

	}

	void FormConnect::on_buttonBox_accepted() {
		if(ui.txtIP->text().length() == 0 || ui.txtPort->text().length() == 0) {
			QMessageBox::warning(this, tr("Missing"), tr("Please enter an IP address and a communication port"));
			show();
		}
		else {
            bool vel, dens;
            vel = (ui.chkVel->checkState() == 2);
            dens = (ui.chkDens->checkState() == 2);
            emit newConnect(ui.txtIP->text().toStdString(), ui.txtPort->text().toInt(), vel, dens, true);
		}
	}

}
