// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014      
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
#ifndef GLNEMOFORMCONNECT_H
#define GLNEMOFORMCONNECT_H

#include "ui_formconnect.h"

namespace glnemo {

	class FormConnect : public QDialog {

		Q_OBJECT

		public:
			FormConnect(QWidget *parent = 0);
			~FormConnect();

		private:
			Ui::FormConnect ui;

		private slots:
			void on_buttonBox_accepted();

		signals:
            void newConnect(std::string adresseIP, int Port, bool velocities, bool densities, bool);

};

}

#endif // GLNEMOFORMCONNECT_H
