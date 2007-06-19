// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2007                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef DRAWBOX_H
#define DRAWBOX_H

#include <qwidget.h>
#include <qpainter.h>
#include <qimage.h>
#include <qcolor.h>
#include "glbox.h"
#include "particles_data.h"
#include "particles_select.h"
#include "global_options.h"

/**
@author Jean-Charles Lambert
*/
// ============================================================================
class ColorRGBA {
public:
  ColorRGBA(float red=0, float green=0, float blue=0, float alpha=1 );
  ~ColorRGBA();
  float red,green,blue;
  float alpha;
  bool enable;
  void clear(float red=0, float green=0, float blue=0, float alpha=1 );
  void blend(float , float , float, float  );
  void blend(QColor c, float  );
  void blendFilter(QColor c, float, float  );
};
// ============================================================================
class DrawBox;
class TexAlpha {
public:
  ~TexAlpha();
  void init(const QImage *,const  int);
  void draw(DrawBox * dbox, const ParticlesSelectVector *,const int, const float, const float, const float); // dbox, obj_index, wx, wy
  int size;           // texture's width and height size             
  int half_size;
  float * alpha_array; // linear 2D array to store texture alpha value
  void print();
  float alpha2[520*520];
private:
  
  void loadAlpha(const QImage*);   // scanline texture to load alpha              
  void resize(float );
};
// ============================================================================
class Texture {
public:
  Texture(QString, int);
  ~Texture();
  void print();
  void draw(DrawBox * dbox, const ParticlesSelectVector *, const int, const float, const float, const int, const float); // dbox, obj_index, wx, wy, sample(size)
private:
  QImage * texture;     // texture's image 
  int sample;           // #textures       
  TexAlpha * tex_alpha; // array of texture
};
// ============================================================================
class DrawBox: public QWidget {
public:
    DrawBox(QWidget *parent=0, const char *name=0 );

    ~DrawBox();
    
    void draw(GLBox * glbox, const ParticlesData * part_data, 
	      const ParticlesSelectVector * psv, const GlobalOptions *, const QString shot_name);
    
    ColorRGBA color_rgba[1024][1024];
    int width,    // windows's width 
    height;       // windows's height
    const GlobalOptions * store_options;
    const ParticlesSelectVector * psv;
    const int * viewport;
protected:
    //void	paintEvent( QPaintEvent * );
private:
 
    //QImage     * texture;
    Texture * texture;
    int tw,th;
    void filterPoint(float x, float y, QColor, float alpha);
    void clearBuffer();
    void paintBuffer(QString);
};
// ============================================================================
#endif
