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
#include "glcolorbar.h"

namespace glnemo {

// ============================================================================
// constructor    
GLColorbar::GLColorbar(const GlobalOptions  * _go,bool _enable ):GLObject()
{
  go     = _go;
  is_activated = _enable;
  font = NULL;
  legend = new GLTextObject(); // new object for text display
  updateFont();
}

// ============================================================================
// destructor                                                                  
GLColorbar::~GLColorbar()
{
  if (font) delete font;
  delete  legend;
}

// ============================================================================
// void updateFont
void GLColorbar::updateFont()
{
  fntRenderer text;
  if (font) delete font;
  font = new fntTexFont(go->gcb_font_name.toStdString().c_str());  
  text.setFont(font);
  text.setPointSize(go->gcb_font_size);  
  legend->setFont(text);
  legend->setColor(go->gcb_color);
}

// ============================================================================
// void update
void GLColorbar::update(GLObjectParticlesVector * _gpv, PhysicalData * _phys_select,
                         GlobalOptions   * _go, QMutex * _mutex)
{
  // update variables
  gpv           = _gpv;
  go            = _go;
  mutex_data    = _mutex;
  phys_select   = _phys_select;
  
}
// ============================================================================
// void GLSelection::display
void GLColorbar::display(const int _width, const int _height)
{
  height = _height;
  width  = _width;
  
  if (go && go->gcb_enable && phys_select && phys_select->isValid()) {
    glDisable( GL_DEPTH_TEST );
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluOrtho2D(0.,width,0.,height);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
        
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);// No Alpha bending accumulation
    glEnable(GL_BLEND);
    // draw box
    drawBox();
    // draw text
    drawLegend();
    glDisable(GL_BLEND);
    // draw color
    drawColor();

    // go back to normal mode
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    glEnable( GL_DEPTH_TEST );
  }
}

// ============================================================================
// void GLColorbar::drawBox
// ***************************************************
// Est and West
// 
//       x=x[0][0]           x=x[1][0]
//       y=x[0][1]           y=x[1][1]
//                +-------+
//                |       |
//                |       |
//                |       |
//                |       |
//                |       |
//                |       |
//                |       |
//                |       |
//                +-------+
//       x=x[3][0]           x=x[2][0]
//       y=x[3][1]           y=x[2][1]
// ***************************************************
// North and South
// 
//       x=x[1][0]                      x=x[2][0]
//       y=x[1][1]                      y=x[2][1]
//                +----------------------+
//                |                      |
//                |                      |
//                |                      |
//                +----------------------+
//       x=x[0][0]                      x=x[3][0]
//       y=x[0][1]                      y=x[3][1]
//
// ***************************************************
//        (0,0) OpenGL   (screen)           (width,height) glText
//          +-------------------------------------+
//          |                                     |
//          |                                     |
//          |                                     |
//          |                                     |
//          |                                     |
//          |                                     |
//          |                                     |
//          |                                     |
//          +-------------------------------------+
//        (0,0) glText                         (width,height) OpenGL
void GLColorbar::drawBox()
{
  switch (go->gcb_orientation) { 
  case 0:// NORTH
    x[0][0] = width/2-go->gcb_pheight*width/2;
    x[0][1] = height-go->gcb_offset-go->gcb_pwidth*height;
    
    x[1][0] = x[0][0];
    x[1][1] = height-go->gcb_offset;
    
    x[2][0] = width/2+go->gcb_pheight*width/2;
    x[2][1] = x[1][1];
    
    x[3][0] = x[2][0];
    x[3][1] = x[0][1];
    break;

  case 1:// EST
    x[0][0] = width-go->gcb_offset-go->gcb_pwidth*width;
    x[0][1] = height/2-go->gcb_pheight*height/2;
    
    x[1][0] = width-go->gcb_offset;
    x[1][1] = x[0][1];
    
    x[2][0] = x[1][0];
    x[2][1] = height/2+go->gcb_pheight*height/2;
    
    x[3][0] = x[0][0];
    x[3][1] = x[2][1];
    break;
  case 2:// SOUTH
    x[0][0] = width/2-go->gcb_pheight*width/2;
    x[0][1] = go->gcb_offset;
    
    x[1][0] = x[0][0];
    x[1][1] = go->gcb_offset+go->gcb_pwidth*height;
    
    x[2][0] = width/2+go->gcb_pheight*width/2;
    x[2][1] = x[1][1];
    
    x[3][0] = x[2][0];
    x[3][1] = x[0][1];
    break;

  case 3:// WEST
    x[0][0] = go->gcb_offset;
    x[0][1] = height/2-go->gcb_pheight*height/2;
    
    x[1][0] = go->gcb_offset+go->gcb_pwidth*width;
    x[1][1] = x[0][1];
    
    x[2][0] = x[1][0];
    x[2][1] = height/2+go->gcb_pheight*height/2;
    
    x[3][0] = x[0][0];
    x[3][1] = x[2][1];
    break;
  default: break;  
  }
  // draw box
  glColor4f( 1.0f, 0.f, 0.f,1.f );
  glBegin(GL_LINE_STRIP);
  glVertex2i(x[0][0],x[0][1]);
  glVertex2i(x[1][0],x[1][1]);
  glVertex2i(x[2][0],x[2][1]);
  glVertex2i(x[3][0],x[3][1]);
  glVertex2i(x[0][0],x[0][1]);
  glEnd();
}
// ============================================================================
// void GLColorbar::drawColor
void GLColorbar::drawColor()
{
  if (go && phys_select && phys_select->isValid()) {
    int large_box,long_box;
    if (go->gcb_orientation==1 || go->gcb_orientation==3) {  // Est or West
      long_box=x[3][1]-x[0][1];// -2 pixels
      large_box=x[2][0]-x[3][0];
    } else {                     // South or North
      long_box=x[3][0]-x[0][0];// -2 pixels
      large_box=x[0][1]-x[1][1];

    }
    int ncolors=go->R->size();
    //ncolors = (go->gcb_max-percmin)*ncolors;
        
    int R,G,B;
    int cpt=0;
    //for (int i=0; i<=vbox;i++) {
    for (int i=0; i<long_box-2;i++) {
      int index;
      if (go->dynamic_cmap) { // dynamic cmap
        if (!go->reverse_cmap)//    normal cmap 
          index=i*ncolors/(long_box-2);
        else                  //    reverse cmap
          index=(long_box-i-2)*ncolors/(long_box-2);
      } else {                // constant cmap
        int ncolors2 = (go->gcb_max-go->gcb_min)/100.*ncolors;
        if (!go->reverse_cmap)//    normal cmap 
          index=(go->gcb_min*ncolors)/100.+i*ncolors2/(long_box-2);
        else                  //    reverse cmap
	  index = ncolors - (go->gcb_min*ncolors)/100.-i*ncolors2/(long_box-2);

        //std::cerr << index << " " << percmin << " " << go->gcb_max << " " << ncolors2 << " "<< ncolors<< "\n";
      }
      
      if (index>=0 && index<ncolors) {
        cpt++;
        R=pow((*go->R)[index],go->powercolor)*255;
        G=pow((*go->G)[index],go->powercolor)*255;
        B=pow((*go->B)[index],go->powercolor)*255;
        // draw color line
        glColor3ub(R,G,B);        
        glBegin(GL_LINES);     
	if (go->gcb_orientation==1 || go->gcb_orientation==3) {   // Est or West
	  glVertex2i(x[0][0]+1,x[0][1]+i+1);
	  glVertex2i(x[1][0]-1,x[1][1]+i+1);
	} else {                      // South or North
	  float fac=-1;
	  //if (place==0) fac=1;
	  glVertex2i(x[0][0]+i+1,x[0][1]-fac);
	  glVertex2i(x[1][0]+i+1,x[1][1]+fac);
	}
        glEnd();
      } else {
      }      
      
    }
    //std::cerr << "cpt =  " << cpt << " " << long_box << " +++ " <<  x[3][0] -  x[0][0] << "\n";
  }
}
// ============================================================================
// void GLColorbar::drawLegend
void GLColorbar::drawLegend()
{
  if (go && phys_select && phys_select->isValid()) {
    float value;
    float diff_rho=(log(phys_select->getMax())-log(phys_select->getMin()))/100.;
    //max
    value=log(phys_select->getMin())+go->gcb_max*diff_rho;
    drawText(value,0);
    //max - 1/3 (max-min)
    value=log(phys_select->getMin())+(go->gcb_max-(go->gcb_max-go->gcb_min)/3.)*diff_rho;
    drawText(value,1);
    // max - 2/3 (max-min)
    value=log(phys_select->getMin())+(go->gcb_max-2.*(go->gcb_max-go->gcb_min)/3.)*diff_rho;
    drawText(value,2);
    // max - 3/3 (max-min) 
    value = log(phys_select->getMin())+go->gcb_min*diff_rho; 
    drawText(value,3);
  }
}
// ============================================================================
// void GLColorbar::drawText
void GLColorbar::drawText(float value, int fac)
{
  QString text1,text0="";
  int xx=0,yy=0,tw=0,th=0;
  // max
  if (!go->gcb_logmode) value=exp(value);
  text1=QString("%1").arg(value,0,'E',go->gcb_ndigits);
  legend->setText(text0,text1);    
  switch (go->gcb_orientation) {
  case 0: // North
    tw=legend->getTextWidth();    
    th=legend->getHeight();    
    xx=x[3][0]-tw/2-fac*(x[3][0]-x[0][0])/3.;
    yy=height-(-5+x[0][1]-th);
    break;
  case 1: // Est
    tw=legend->getTextWidth();
    xx=x[0][0]-tw-5;
    yy=5+x[0][1]+fac*(x[3][1]-x[0][1])/3.;
    break;
  case 2: // South
    tw=legend->getTextWidth();    
    xx=x[2][0]-tw/2-fac*(x[2][0]-x[1][0])/3.;
    yy=height-(5+x[1][1]);
    break;
  case 3: // West
    xx=5+x[1][0];
    yy=5+x[1][1]+fac*(x[2][1]-x[1][1])/3.;
    break;
  }
  legend->setPos(xx,yy,xx);
  font->begin(); // mandatory !!
  legend->display(width,height);
  font->end();   // mandatory !!
}
} // namespace glnemo
