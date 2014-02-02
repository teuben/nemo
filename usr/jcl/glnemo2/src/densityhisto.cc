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
#include "densityhisto.h"
#include <QPainter>
#include <QRadialGradient>
#include <QGraphicsItem>
#include <iostream>
#include <cmath>
#include <cstring>
namespace glnemo {

// ============================================================================
// Constructor                                                                 
DensityHisto::DensityHisto(QWidget *_parent):parent(_parent),border(0),nhisto(300)
{
  //drawGrid();

 // a white semi-transparent foreground
//setBackgroundBrush(QColor(Qt::darkYellow));
  
  // reserve vector
  density_histo.reserve(nhisto);
 // a grid foreground
 setBackgroundBrush(QBrush(Qt::lightGray, Qt::CrossPattern));
}

// ============================================================================
// Constructor                                                                 
DensityHisto::~DensityHisto()
{
    clearScene();
}
void DensityHisto::resizeEvent ( QResizeEvent * event )
{
  if (event) {;}
  //std::cerr << "DensityHisto::resizeEvent resize parent->width =" << parent->width() << " my width="<<width()<<"\n";
}
// ============================================================================
// drawGrid                                                                    
void DensityHisto::drawGrid()
{
  QPainterPath path;

  
  // vertical lines
  for (int i=0; i<=10; i++) {
    int offset=(parent->width()-border)/10;
    path.moveTo(0+border/2+i*offset, 0);
    path.lineTo(0+border/2+i*offset, parent->height());
  }
  // horizontal lines
  for (int i=0; i<=10; i++) {
    int offset=(parent->width()-border)/10;
    path.moveTo(0,0+border/2+i*offset);
    path.lineTo(parent->width(),0+border/2+i*offset);
  }

  addPath(path);
}
// ============================================================================
// drawDensity  
// density_histo stores the number of particles foreach percentage
// of density from 0 to 99 %
void DensityHisto::drawDensity(const std::vector <int> _density_histo) //_density_histo[100])
{
  //memcpy(density_histo,_density_histo,sizeof(int)*100);
  density_histo = _density_histo;
  // compute maxhisto
  for (int i=0; i<nhisto; i++) {
    maxhisto=std::max(maxhisto,density_histo[i]);
  }
  // draw density curve
  drawDensity();
}
// ============================================================================
// drawDensity                                                                 
void DensityHisto::drawDensity(int _min, int _max)
{
  //clearScene();
  clear();
  // draw density curve
  QPainterPath path1;
  for (int i=0; i<nhisto; i++) {
    int y;
    if (density_histo[i] > 0) {
      y=parent->height()-log(density_histo[i])*(parent->height()-border)/log(maxhisto);
    } else {
      y = parent->height()-border;
    }
    int x=i*parent->width()/nhisto;
    if (i==0) {
      path1.moveTo(x+border/2,y);
    } else {
      path1.lineTo(x+border/2,y);
    }
  }
  addPath(path1,QPen(Qt::blue));
  QPainterPath path2;
  int x,y;

  // draw histogram between min and max
  int lastx=-1;
  for (int i=_min*nhisto/100; i<_max*nhisto/100; i++) {
    if (density_histo[i] > 0) {
      y=parent->height()+1-log(density_histo[i])*(parent->height()-border)/log(maxhisto);
    } else {
      y = parent->height()-border;
    }
    x=i*parent->width()/(float) (nhisto);
    if (i ==_min*nhisto/100) {
      lastx=x;
    }
    for (int xx=lastx+1; xx<=x; xx++) {
      path2.moveTo(xx+border/2,parent->height()-border);
      path2.lineTo(xx+border/2,y);
    }
    //std::cerr << "x="<<x<<" float="<<i*parent->width()/100.<<"\n";
    lastx=x;
  }

  QPen pen(QColor(255,0,0,127));
  //pen.setWidth(2);
  addPath(path2,pen);
}
// ============================================================================
// clearScene
// Clear the scene by removing all the items. I used to used clear() method, but
// since qt 4.5 it seems there is a bug in it...
void DensityHisto::clearScene()
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
}
