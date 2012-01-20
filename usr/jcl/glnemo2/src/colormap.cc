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
#include "colormap.h"
#include <sstream>
#include <QFile>
#include <QTextStream>

namespace glnemo {
// ============================================================================
// Constructor                                                                 
Colormap::Colormap(GlobalOptions * _go)
{
  go = _go;
  load(go->colormap);
}
// ============================================================================
// Destructor                                                                  
Colormap::~Colormap()
{
}
// ============================================================================
// load                                                                        
int Colormap::load(const int _cmap)
{
  cmap = _cmap;
  load();
  return cmap;
}
// ============================================================================
// load                                                                        
int Colormap::load()
{
  if (cmap > 125) cmap=101;
  else  if (cmap < 101) cmap=125;

  std::ostringstream fortmap;
  fortmap << GlobalOptions::RESPATH.toStdString()+"/colormaps/fort." << cmap;
  //QFile file(":/colormaps/fort.108");  //fort.106" 108 117
  QFile file((fortmap.str()).c_str());  //fort.106" 108 117
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    return 0;
  QTextStream in(&file);
  QString line;
  R.clear();
  G.clear();
  B.clear();

  int cpt=1;
  do {
    line = in.readLine();
    if (!line.isNull()) {
      //std::cerr << "line :" << line.toStdString() <<"\n";
      std::istringstream ss(line.toStdString());
      float r,g,b;
      ss >> r;
      ss >> g;
      ss >> b;
      R.push_back(r);
      G.push_back(g);
      B.push_back(b);
      cpt++;
      //std::cerr << "R :" << r <<"G : " << g << "B : "<< b <<"\n";
    }
  } while (!line.isNull());
  file.close();
  go->R = &R;
  go->G = &G;
  go->B = &B;
  emit newColorMap();
  return cpt;
}
// ============================================================================
// next                                                                        
int Colormap::next()
{
  load(cmap+1);
  return (cmap+1);
}
// ============================================================================
// prev                                                                        
int Colormap::prev()
{
  load(cmap-1);
  return (cmap-1);
}
}
