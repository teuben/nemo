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
#include "formabout.h"
#include <QPainter>
#include <iostream>
#include "globaloptions.h"

namespace glnemo {

FormAbout::FormAbout(QWidget *parent):QDialog(parent)
{
  if (parent) {;}  // remove compiler warning
  //this->setPalette(parent->palette());
  form.setupUi(this);
  //QFont font("Times",14);
  //QWidget::setFont(font);
  //"<b>Glnemo2:</b> an interactive 3D visualisation program for nbody simulation data"
  
  QString info(tr(
       "<center>Copyright (c) <b>Jean-Charles LAMBERT</b> 2007-2012"
      "<br><a href=\"mailto:Jean-Charles.Lambert@oamp.fr\">Jean-Charles.Lambert@oamp.fr</a>"
       "<br><br>Dynamique des Galaxies"
       "<br>Centre de donnéeS Astrophysique de Marseille (CeSAM)"
       "<br><br>Aix Marseille Université, CNRS"
       "<br>Laboratoire d'Astrophysique de Marseille (LAM)"
       "<br>UMR 7326"
       "<br>38, rue Frédéric Joliot-Curie"
       "<br>13388 Marseille cedex 13"
       "<br>FRANCE"
                 "<br><br><a href=http://projets.oamp.fr/projects/glnemo2>GLNEMO2 home page</a>"
       "</center>"));
  form.text_info->setHtml(info);
  QRect geo = form.view_picture->geometry();       // initial geometry defined with designer
  QPixmap pix =  QPixmap(GlobalOptions::RESPATH+"/images/glnemo2.png");  // picture to display
  scene.setSceneRect(pix.rect());                  // size of the scene
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);  
  scene.addPixmap(pix);                            // add picture
  form.view_picture->setGeometry(geo.x(),geo.y(),pix.width(),pix.height()); // rebuild geometry
  form.view_picture->updateGeometry(); // mandatory because FIXED size policy
  form.view_picture->setScene(&scene);
  form.tabWidget->setCurrentIndex(0); // set position to first tab
}


FormAbout::~FormAbout()
{
}

}
