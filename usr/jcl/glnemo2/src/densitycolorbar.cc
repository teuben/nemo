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
#include "densitycolorbar.h"
#include <QPainter>
#include <QList>
#include <QGraphicsItem>
#include <cmath>
#include <iostream>
namespace glnemo {

// ============================================================================
// Constructor                                                                 
  DensityColorBar::DensityColorBar(GlobalOptions * _go,QWidget *_parent):parent(_parent)
{
  go = _go;
}
// ============================================================================
// Destructor                                                                  
DensityColorBar::~DensityColorBar()
{
    clearScene();
}
// ============================================================================
// Draw                                                                        
void DensityColorBar::draw()
{
  if (go) {
    int ncolors=go->R->size();
    clearScene();
    int R,G,B;
    for (int i=0; i<parent->width();i++) {
      //QPainterPath path;
      int index;
      if (!go->reverse_cmap) 
        index=i*(ncolors-1)/(parent->width()-1);
      else
        index=(parent->width()-1-i)*(ncolors-1)/(parent->width()-1);
      
      //path.moveTo(i,0);
      //path.lineTo(i,parent->height());
      R=pow((*go->R)[index],go->powercolor)*255;
      G=pow((*go->G)[index],go->powercolor)*255;
      B=pow((*go->B)[index],go->powercolor)*255;
      addLine(i,0,i,parent->height(),QPen(QColor(R,G,B)));
      //addPath(path,QPen(QColor(R,G,B)));
    }
  }
}
// ============================================================================
// Draw                                                                        
void DensityColorBar::draw(const int min, const int max)
{
  if (go && !go->dynamic_cmap) {
    draw();  // constant colormap
  } 
  else {     // Dynamic colormap 
    clearScene();
    if (go) {
        clearScene();
        int ncolors=go->R->size();
        setSceneRect(0,0,parent->width(),parent->height());
        //QPainterPath path;
        int R,G,B;
        
        // loop from min to max range
        for (int i=min*parent->width()/100.; i<max*parent->width()/100.;i++) {

          //int index=(i-min*parent->width()/100.)*ncolors/((max-min)*parent->width()/100.);
          int index=(i-min*parent->width()/100.)*ncolors/((max-min)*parent->width()/100.);
          if (go->reverse_cmap) {
             index = ncolors-1-index;
          }
          if (index>=0 && index<ncolors) {
            R=pow((*go->R)[index],go->powercolor)*255;
            G=pow((*go->G)[index],go->powercolor)*255;
            B=pow((*go->B)[index],go->powercolor)*255;
            addLine(i,0,i,parent->height(),QPen(QColor(R,G,B)));
          }
        }
        //addPath(path,QPen(QColor(R,G,B)));
   }
 }
}
// ============================================================================
// clearScene
// Clear the scene by removing all the items. I used to used clear() method, but
// since qt 4.5 it seems there is a bug in it...
void DensityColorBar::clearScene()
{
   QList<QGraphicsItem *> list=items(); // get list of items
   if (list.size() > 0) {
     for (QList<QGraphicsItem *>::iterator it=list.begin(); it!=list.end();it++) {
         removeItem(*it);
         delete (*it);
     }
     list.clear();
     //clear();
   }
}
// ============================================================================
// resizeEvent
void DensityColorBar::resizeEvent ( QResizeEvent * event )
{
  if (event) {;}
  setSceneRect(0,0,parent->width(),parent->height());
  //std::cerr << "DensityHisto::resizeEvent resize parent->width =" << parent->width() << " my width="<<width()<<"\n";
}
}
