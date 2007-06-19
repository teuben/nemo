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
#include <qpen.h>
#include <qpicture.h>
#include <qstring.h>
#include <math.h>
#include <qcolor.h>
#include "drawbox.h"
#include "smoke.h"
#include <assert.h>

#define MAX_TEXTURES 512

// ============================================================================
DrawBox::DrawBox( QWidget *parent, const char *name)
  : QWidget( parent, name, WStaticContents )
{
  setBackgroundColor( black );		// white background

  std::cerr << "Generating texture....\n";
  texture = new Texture("smoke",MAX_TEXTURES);
  std::cerr << "End of Generating texture....\n";
  //texture->print();
}

// ============================================================================
DrawBox::~DrawBox()
{
  delete texture;
}
// ============================================================================
// DrawBox::draw()                                                             
// convert particles from 3D space to 2D screen coordinates                    
// compute blending values for particles and texture                           
void DrawBox::draw( GLBox * glbox, const ParticlesData * part_data, 
		    const ParticlesSelectVector * _psv, const GlobalOptions * _store_options,
		    QString shot_name)
{
  const double * mProj    = glbox->getProjMatrix();
  double * mMod     = const_cast<double*> 
	(glbox->getModelMatrix());
  // get data pointers
  viewport      = glbox->getViewPort();
  store_options = _store_options;
  psv           = _psv;
#define MP(row,col)  mProj[col*4+row]
#define MM(row,col)  mMod[col*4+row]

  // get windows size 
  width = viewport[2];
  height= viewport[3];
  // set windows size
  setGeometry(viewport[0],viewport[1],viewport[2],viewport[3]);
  
  //QPainter paint(this);
  //paint.eraseRect(viewport[0],viewport[1],viewport[2],viewport[3]);
  clearBuffer();
  
  if (psv->size()) {
    // Zoom is located in ModelView matrix at coordinates MM(2,3)
    // loop on all object
    for (int i=0; i< (int ) psv->size(); i++ ) {
      // loop on all visible object selected particles  
      //std::cerr << "----------\n";  
      if ((*psv)[i].vps->is_visible ) {
	for (int j  = 0; 
	     j  <  (*psv)[i].vps->ni_index;
	     j ++) {
	  int jndex= (*psv)[i].vps->index_tab[j];      
          float 
              x=part_data->pos[jndex*3  ],
	      y=part_data->pos[jndex*3+1],
	      z=part_data->pos[jndex*3+2],
	      w=1;   
          // do the product Mmodel X point = mxyzw
          float myzw0 = MM(0,1)*y + MM(0,2)*z + MM(0,3)*w;
          float mx = MM(0,0)*x  + myzw0;
          float mx1= mx + store_options->texture_size/2.; // move X from half a texture
          
          float myzw1 = MM(1,1)*y + MM(1,2)*z + MM(1,3)*w;
          float my = MM(1,0)*x  + myzw1;
          
          float MmzG = MM(2,1)*y + MM(2,2)*z;
          float Mmz = MM(2,0)*x  + MmzG;
          
          float mz = Mmz  + MM(2,3)*w;
          
          float mwG= MM(3,1)*y + MM(3,2)*z + MM(3,3)*w;
          float mw = MM(3,0)*x  + mwG;
          
          // do the product Mproj X mxyzw  = pxyzw
          float Ppx = MP(0,0)*mx + MP(0,1)*my + MP(0,3)*mw;
          float px  = Ppx + MP(0,2)*mz;
          float px1 = MP(0,0)*mx1 + MP(0,1)*my + MP(0,3)*mw + MP(0,2)*mz;
          
          float Ppy= MP(1,0)*mx + MP(1,1)*my + MP(1,3)*mw;
          float py = Ppy + MP(1,2)*mz;
          
          float Ppz= MP(2,0)*mx + MP(2,1)*my + MP(2,3)*mw;
          float pz = Ppz + MP(2,2)*mz;
          
          float Ppw= MP(3,0)*mx + MP(3,1)*my + MP(3,3)*mw;
          float pw = Ppw + MP(3,2)*mz;
          
          float pw1 = MP(3,0)*mx1 + MP(3,1)*my + MP(3,3)*mw + MP(3,2)*mz;
          // normalyze
	  //std::cerr << px << " " << py << " " << pz << "\n";
          px /= pw;
          px1/= pw1;
          py /= pw;
	  pz /= pw;
	  //std::cerr << px << " " << py << " " << pz << "\n";
          // compute screen coordinates
          float winx =viewport[0] + (1 + px)  * viewport[2] / 2;
          float winx1=viewport[0] + (1 + px1) * viewport[2] / 2;
          float winy =viewport[1] + (1 + py)  * viewport[3] / 2;
          // paint particles
	  //paint.setPen(red);
	  //paint.drawPoint((int) (winx), viewport[3]-(int) (winy));

          #define MMAX(A,B) ((A)>(B)?(A):(B))
          #define MMIN(A,B) ((A)<(B)?(A):(B))

          //if (winx >= 0. && winy >= 0 && winx < 1024 && winy < 1024) {
	  if (winx >= 0. && winy >= 0 && winx < width && winy < height && pz<1.0) {
	    float size_texture = (int) fabs(winx1-winx);
            //std::cerr << "x=" << winx << "  x1=" << winx1 << " size="<<size_texture<< "\n";
	    
	    //color_rgba[(int) (winx)][(int) (winy)].blend(psv[i].vps->col,(float) (alpha)/255.);
	    if (store_options->show_part) {
	      //color_rgba[(int) (winx)][viewport[3]-1-(int) (winy)].blend((*psv)[i].vps->col,
	      //(float) (store_options->particles_alpha)/255.);
	      filterPoint( winx,winy,(*psv)[i].vps->col,(float) (store_options->particles_alpha)/255.);
	    }
	  //paint.drawPoint((int) (winx), (int) (winy));
	    if (store_options->show_poly) {
	      float s4 = size_texture*4.;
	      float ratio_texture =  s4 - floor(s4);
	      float ceil_texture  = ceil(s4);
	      texture->draw(this,psv,i,winx,winy,(int)(ceil_texture), ratio_texture); 
	    }
          }
	}
      }
    }
  }
  
  paintBuffer(shot_name);
}
// ============================================================================
// DrawBox::filterPoint()                                                      
//                                                                             
void DrawBox::filterPoint( float x, float y, QColor c,float alpha)
{
  float u  = x-floor(x);
  float uo = 1. - u;
  float v  = y-floor(y);
  float vo = 1. - v;
  int a = (int) floor(x);
  int b = (int) floor(y);
  int off = viewport[3]-1;
  //assert((off-b) > 0);
  //assert((off-b-1) > 0);	
  if (a < 1024 && ((off-b) < 1024) && ((off-b) >= 0))
    color_rgba[a][off-b].blendFilter(c,uo*vo,alpha);
  if ((a+1) < 1024 && ((off-b) < 1024 )&& ((off-b) >= 0))
    color_rgba[a+1][off-b].blendFilter(c,u*vo,alpha);
  if (a < 1024 && (off-b-1) < 1024 && ((off-b-1) >= 0) )
    color_rgba[a][off-b-1].blendFilter(c,uo*v,alpha);
  if ((a+1) < 1024 && (off-b-1) < 1024 && ((off-b-1) >= 0))
    color_rgba[a+1][off-b-1].blendFilter(c,u*v,alpha);
 
 
}
// ============================================================================
// DrawBox::clearBuffer()                                                      
// clear RGBA buffer with black color                                          
void DrawBox::clearBuffer()
{ 
  for (int i=0; i<1024; i++) {
    for (int j=0; j<1024; j++) {
      color_rgba[i][j].clear();
    }
  }
}
// ============================================================================
// DrawBox::paintBuffer()                                                      
// paint in a pixmap the RGBA buffer                                           
void DrawBox::paintBuffer(QString shot_name)
{
  //const int    * viewport = (const int *) glbox->getViewPort();
  QPixmap  pic(width,height);
  //QPainter paint(this);
  QPainter paint(&pic);
  //paint.end();
  paint.eraseRect(viewport[0],viewport[1],viewport[2],viewport[3]);
  //paint.begin(&pic);
  for (int i=0; i<viewport[2]; i++) {
    for (int j=0; j<viewport[3]; j++) {
#if 1
      QColor c((int) color_rgba[i][j].red    ,
            (int) color_rgba[i][j].green  ,
            (int) color_rgba[i][j].blue  );
      if (color_rgba[i][j].red < 0) std::cerr << "red 255\n";
      if (color_rgba[i][j].green < 0) std::cerr << "green 255\n";
      if (color_rgba[i][j].blue < 0) std::cerr << "blue 255\n";
#endif    

        QPen pen;
        pen.setColor(c);
        paint.setPen(pen);
        //paint.drawPoint(i, viewport[3]-j-1);
	paint.drawPoint(i, j);

      if (color_rgba[i][j].enable) {
        //paint.drawImage(i-tw, viewport[3]-j-th,*texture,DiffuseAlphaDither|AlphaDither_Mask);
      }
    }
  }
  //paint.end();
  if (shot_name != "") {
    std::cerr << "Render soft: Saving frame [" << shot_name << "]\n";
    if (!pic.save(shot_name,"PNG")) {  // save image on disk
      std::cerr << "DrawBox::paintBuffer, UNABLE to save frame......\n";
    }
  }
  QPainter p(this);       
  p.drawPixmap(0,0,pic);       // draw image on screen
  
}
// ============================================================================
ColorRGBA::ColorRGBA(float _red, float _green, float _blue, float _alpha)
{
  clear(_red, _green, _blue, _alpha);
}
// ============================================================================
ColorRGBA::~ColorRGBA()
{

}
// ============================================================================
void ColorRGBA::clear(float _red, float _green, float _blue, float _alpha)
{
  red   = _red;
  green = _green;
  blue  = _blue;
  alpha = _alpha; 
  enable=false;
}
// ============================================================================
void ColorRGBA::blendFilter(QColor c, float filter,float _alpha)
{
  c.setRgb((int) (c.red()*filter),(int) (c.green()*filter),(int)(c.blue()*filter));
  blend(c,_alpha);
}
// ============================================================================
void ColorRGBA::blend(QColor c, float _alpha)
{
  blend((float) c.red(), (float) c.green(), (float) c.blue(), _alpha);
}
// ============================================================================
void ColorRGBA::blend(float _red, float _green, float _blue, float _alpha)
{
#if 0
  red   = _alpha * _red   + ( 1 - _alpha) * red;
  green = _alpha * _green + ( 1 - _alpha) * green;
  blue  = _alpha * _blue  + ( 1 - _alpha) * blue;
  //alpha = _alpha * _alpha + (1  - _alpha) * alpha; 
#endif
  #define LIMIT255(A) ((A>=255.)?255:A)
  red   = _red   * _alpha  + red;     //+ ( 1 - alpha) * red;
  red   = LIMIT255(red);
  
  green = _green * _alpha  + green;   //+ ( 1 - alpha) * green;
  green = LIMIT255(green);
  
  blue  = _blue  * _alpha  + blue;    //+ ( 1 - alpha) * blue;
  blue  = LIMIT255(blue);
  //alpha =_alpha;
#if 0
  std::cerr << red << " aared 255\n";
  std::cerr << green <<" aagreen 255\n";
  std::cerr << blue <<" aablue 255\n";
#endif  
enable=true;
}
// ============================================================================
//                  T   E   X    T    U    R     E     S                       
// ============================================================================
Texture::Texture(QString name, int _sample)
{
  // load texture from disk or embeded    
#if 0
  texture= new QImage(qembed_findImage(name));
#else  
  texture= new QImage("/home/jcl/text1.png");
#endif
  if (name) ; // remove compiler warning
  sample = _sample;                 // number of samples for the texture
  
  tex_alpha = new TexAlpha[sample]; // allocate array to store textures's alpha channel;
  // initialyse texture
  for (int i=0; i<sample; i++) {
    tex_alpha[i].init(texture,i+1);
  }    
  //tex_alpha[50].print();
}
// ============================================================================
Texture::~Texture()
{
  if (sample) {
/*    for (int i=0; i<sample; i++) {
      delete tex_alpha[i];
    }*/
    delete [] tex_alpha;
    delete texture;
  }
}
// ============================================================================
void Texture::print()
{
  for (int i=0; i<sample; i++) {
    tex_alpha[i].print(); 
  }
}

// ============================================================================
// Texture::draw()                                                             
// draw texture in the RGBA buffer                                             
void Texture::draw(DrawBox * dbox, const ParticlesSelectVector * psv,const int obj_index, 
		   const float wx, const float wy, const int _sample, const float ratio_texture)
{
  int size = _sample;
  if (size >= sample) size = MAX_TEXTURES-1;
  if (size == 0     ) size = 1;
  tex_alpha[size].draw(dbox, psv,obj_index, wx, wy, ratio_texture);
}
// ============================================================================
TexAlpha::~TexAlpha()
{
  delete [] alpha_array;
}
// ============================================================================
// TexAlpha::init()                                                            
// re-scale texture's size and load alpha channel                              
void TexAlpha::init(const QImage * texture,const int _size)
{  
  size = _size;                            // news texture's size 
  half_size = size >> 1;                   // half size           
  QImage new_texture;                      // temporary texture   
  new_texture = texture->scale(size,size); // rescale with size   
  alpha_array = new float[size*size];       // array to store alpha
  loadAlpha(&new_texture);
}
// ============================================================================
// TexAlpha::loadAlpha()                                                       
// load texture's alpha channel                                                
void TexAlpha::loadAlpha(const QImage * texture)
{
  int h=texture->height(),
  w=texture->width();
  assert(size*size==(h*w));
  for ( int y=0; y<size; y++ ) {                 // set image pixels
    QRgb * p = (QRgb *) texture->scanLine(y);
    for ( int x=0; x<size; x++ ) {
      assert( qAlpha((QRgb ) *p) >= 0);
      alpha_array[y*size + x] = qAlpha((QRgb ) *p)/255.;
    //std::cerr << qAlpha((QRgb ) *p)  << " ";
      if (alpha_array[y*size + x] < 0) {
	
	//alpha_array[y*size + x] = (char) fake;
	std::cerr << "negatif : size = " << size << " : " << (int) alpha_array[y*size + x] << " qAlpha =" << qAlpha((QRgb ) *p) <<"\n";
      }
      p++;
    }
    //std::cerr << "\n";
  }
}
// ============================================================================
// TexAlpha::resize()                                                            
void TexAlpha::resize( float _size)
{
#if 1

  float offset = (size-1.+_size)/(float) size;
  // clear texture
  for ( int y=0; y<size; y++ ) {
    for ( int x=0; x<size; x++ ) {
      TexAlpha::alpha2[y*size+x] = 0;
    }
  }
  for ( int y=0; y<size-1; y++ ) {
    for ( int x=0; x<size-1; x++ ) {
      float u  = (x+1)*offset-floor((x+1)*offset);
      float uo = 1. - u;
      float v  = u;
      float vo = uo;
      alpha2[y*size+x]      += alpha_array[y*size+x]      *uo *vo;
      alpha2[y*size+x+1]    += alpha_array[y*size+x+1]    *u  *vo;
      alpha2[(y+1)*size+x]  += alpha_array[(y+1)*size+x]  *uo *v ;
      alpha2[(y+1)*size+x+1]+= alpha_array[(y+1)*size+x+1]*u  *v ;
    }
  }
  if (size==1)
      alpha2[0] = 1.0;//alpha_array[0];
#endif
}
// ============================================================================
// TexAlpha::draw()                                                            
// draw texture in the RGBA buffer                                             
void TexAlpha::draw(DrawBox * dbox, const ParticlesSelectVector * psv,const int obj_index, 
		   const float wx, const float wy, const float ratio_texture)
{
#define T2D(A,B) (alpha_array[A*size+B]);
#define X1 ((wx-half_size)<0 ? 0 : wx-half_size)
#define X2 ((wx+half_size)>dbox->width ? dbox->width : wx+half_size)
#define Y1 ((dbox->viewport[3]-1-wy-half_size)<0 ? 0 : dbox->viewport[3]-1-wy-half_size)
#define Y2 ((dbox->viewport[3]-1-wy+half_size)>dbox->height ? dbox->height: dbox->viewport[3]-1-wy+half_size)

int x1,x2,y1,y2,tx1,ty1;
if ((wx-half_size)<0) {
  x1  = 0;
  tx1 =(int) fabs(wx-half_size);
} else {
  x1  = (int) (wx-half_size);
  tx1 = 0;
}

if ((wx+half_size)>dbox->width) {
  x2  = dbox->width;
} else {
  x2  = (int) (wx+half_size);
}

if (dbox->viewport[3]-1-wy-half_size<0) {
  y1  = 0;
  ty1 = (int) fabs(dbox->viewport[3]-1-wy-half_size);
} else {
  y1  = (int) (dbox->viewport[3]-1-wy-half_size);
  ty1 = 0;
}

if ((dbox->viewport[3]-1-wy+half_size)>dbox->height) {
  y2  = dbox->height;
} else {
  y2  = (int) (dbox->viewport[3]-1-wy+half_size);
}

  float u  = wx-floor(wx);
  float uo = 1. - u;
  float v  = wy-floor(wy);
  float vo = 1. - v;
  //int a = (int) floor(wx);
  //int b = (int) floor(wy);

  //std::cerr << size << " " << wx << " " << wy << " " << X1 << " " << X2 << " " << Y1 << " " << Y2 << "\n";
  float alpha_user = (float) (dbox->store_options->texture_alpha_color)/255.;
  QColor col = (*psv)[obj_index].vps->col;
#if 0  
  int * alpha_array2 = new int[size*size];
  alpha_array2 = alpha_array;
  // reset row 0 and column size-1
  for (int i=0; i< size; i++) {
    alpha_array2[i] = 0;
    alpha_array2[i*size+size-1] = 0;
  }
#endif

  resize( ratio_texture );
  
  for (int i=x1; i<x2; i++) {
    for (int j=y1; j<y2; j++) {
      //assert(((j-Y1)*size+(i-X1))<size*size);
      //float alpha_texture = (float)((int)(alpha_array[(j-Y1)*size+(i-X1)]))/255.;

      //float alpha_texture = (float)((int)(alpha_array[(j-y1+ty1)*size+(i-x1)+tx1]))/255.;

      float alpha_texture = alpha2[(j-y1+ty1)*size+(i-x1)+tx1];
      
#if 0
      dbox->color_rgba[i][j].blend(col,uo*vo*alpha_user*alpha_texture);
      dbox->color_rgba[i+1][j].blend(col,u*vo*alpha_user*alpha_texture);
      dbox->color_rgba[i][j+1].blend(col,uo*v*alpha_user*alpha_texture);
      dbox->color_rgba[i+1][j+1].blend(col,u*v*alpha_user*alpha_texture);
#else
      dbox->color_rgba[i  ][j  ].blendFilter(col,uo*vo,alpha_user*alpha_texture);
      dbox->color_rgba[i+1][j  ].blendFilter(col,u *vo,alpha_user*alpha_texture);
      dbox->color_rgba[i  ][j+1].blendFilter(col,uo*v ,alpha_user*alpha_texture);
      dbox->color_rgba[i+1][j+1].blendFilter(col,u *v ,alpha_user*alpha_texture);
#endif            
      if ((alpha_user*alpha_texture) > 1.) {
	std::cerr << alpha_user*alpha_texture << "\n";
      };
    }
  }

}
// ============================================================================
void TexAlpha::print()
{
  std::cerr << "TEXTUR #" << size << "\n";
  for ( int y=0; y<size; y++ ) {                 // set image pixels
    for ( int x=0; x<size; x++ ) {
      std::cerr << " " << (int) alpha_array[y*size + x];
    }
    std::cerr << "\n";
  }  
}

