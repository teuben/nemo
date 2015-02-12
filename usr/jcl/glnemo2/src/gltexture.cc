// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015                                  
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
#include <GL/glew.h>
#include "gltexture.h"
#include "globaloptions.h"
#include <QImage>
#include <math.h>
#include <iostream>
namespace glnemo {
QString GLTexture::TEXTURE[][2] = {
  {"gaussian"              , "gaussian"      },
  {"/textures/text1.png"   , "particles"     },
  {"/textures/smoke10.png" , "gas"           },
  {"/textures/text2.png"   , "you are here"  },
  {"/textures/identite.jpg", "jcl's identity"},
  {NULL                            , NULL            }
};
// ============================================================================
// Constructor                                                                 
GLTexture::GLTexture()
{
  status=false;
}

// ============================================================================
// Destructor                                                                  
GLTexture::~GLTexture()
{
  //glDeleteTextures(1,&texture);
}
// ============================================================================
// load()                                                                      
// load the given texture from embeded memory or from disk                     
bool GLTexture::load(QString _texture_name, QString _path,bool embeded)
{
  status=false;
  if (embeded) {;} // ! remove compiler warning
  if (_path != NULL && _path != "none")
    _texture_name = _path + _texture_name;

  if (_texture_name == texture_name) {
    std::cerr << "Texture name already exist :" << texture_name.toStdString() << "\n";
    status = true; // already loaded
  }
  else {
    if (status) { // texture already loaded
      std::cerr << "Get rid off previous texture\n";
      glDeleteTextures(1,&texture);
    }
    texture_name = _texture_name;
    if ( _texture_name == GlobalOptions::RESPATH+"gaussian") {
      createGaussian(1024);
      u_max = v_max = 1.;
      status = true;
    } else {
      QImage img(texture_name);
      if (img.isNull()) {
	status = false;
	std::cerr << "Unable to load " << texture_name.toStdString() <<
	",unsupported file format" << std::endl;
      }
      else {
	// 1E-3 needed. Just try with width=128 and see !
	int newWidth  = 1<<(int)(1+log(img.width() -1+1E-3) / log(2.0));
	int newHeight = 1<<(int)(1+log(img.height()-1+1E-3) / log(2.0));
	// compute uv max
	u_max = img.width()  / (float)newWidth;
	v_max = img.height() / (float)newHeight;
        u_max = v_max= 1.0;
	if ((img.width()!=newWidth) || (img.height()!=newHeight)) {
	  //PRINT_D std::cout << "Image size set to " << newWidth
	  //                  << "x" << newHeight << " pixels" << std::endl;
	  img = img.copy(0, 0, newWidth, newHeight);
	}
	//tratio = newWidth / (float)newHeight;
	QImage glImg = QGLWidget::convertToGLFormat(img);  // flipped 32bit RGBA
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1,&texture);             // Create The Texture
	//std::cerr << "TEXTURE Id =" << texture << "\n";
	// Typical Texture Generation Using Data From The Bitmap
	::glBindTexture(GL_TEXTURE_2D, texture);
	// Bind the img texture...
	
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);	

        glTexImage2D(GL_TEXTURE_2D, 0, 4, glImg.width(), glImg.height(), 0,
        GL_RGBA, GL_UNSIGNED_BYTE, glImg.bits());

	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glDisable(GL_TEXTURE_2D);
	status = true;
      }
    }
  }
  return status;
}
// ============================================================================
//
void GLTexture::createGaussian(int resolution)
{
  unsigned char* data = createGaussianMap(resolution);
  glGenTextures(1, &texture);
  ::glBindTexture(GL_TEXTURE_2D, texture);

 glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP_SGIS, GL_TRUE);
 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, resolution, resolution, 0, 
	       GL_RGBA, GL_UNSIGNED_BYTE, data);
}
// ============================================================================
//
unsigned char * GLTexture::createGaussianMap(int N)
{
  float *M = new float[2*N*N];
  unsigned char *B = new unsigned char[4*N*N];
  float X,Y,Y2,Dist;
  float Incr = 2.0f/N;
  int i=0;  
  int j = 0;
  Y = -1.0f;
    //float mmax = 0;
  for (int y=0; y<N; y++, Y+=Incr)
  {
    Y2=Y*Y;
    X = -1.0f;
    for (int x=0; x<N; x++, X+=Incr, i+=2, j+=4)
    {
      Dist = (float)sqrtf(X*X+Y2);
      if (Dist>1) Dist=1;
      M[i+1] = M[i] = evalHermite(1.0f,0,0,0,Dist);
      B[j+3] =  (unsigned char)(M[i] * 255);
      B[j+2] = B[j+1] = B[j] =255;
    }
  }
  delete [] M;
  return(B);
}
// ============================================================================
//
float GLTexture::evalHermite(float pA, float pB, float vA, float vB, float u)
{
  float u2=(u*u), u3=u2*u;
  float B0 = 2*u3 - 3*u2 + 1;
  float B1 = -2*u3 + 3*u2;
  float B2 = u3 - 2*u2 + u;
  float B3 = u3 - u;
  return( B0*pA + B1*pB + B2*vA + B3*vB );
}
// ============================================================================
// int glBindTexture()                                                         
int GLTexture::glBindTexture()
{
 if (status) {
  ::glBindTexture(GL_TEXTURE_2D, texture); // Select texture
 }
 return status;
}
// ============================================================================
// int loadTextureVector()                                                     
// Load embeded texture in a vector                                            
int GLTexture::loadTextureVector(GLTextureVector & gtv)
{
  // delete previous binded texture
  for (GLTextureVector::iterator iv=gtv.begin(); iv<gtv.end();iv++) {
    iv->deleteTexture();
  }
  gtv.clear(); // delete texture vector
  int i=0;
  // loop and load all embeded texture
  while (GLTexture::TEXTURE[i][0]!=NULL) {
    GLTexture * p = new GLTexture();
    p->load(QString(GlobalOptions::RESPATH+TEXTURE[i][0]),NULL);
    gtv.push_back(*p);
    delete p;
    i++;
  }
  return gtv.size();
}
}
