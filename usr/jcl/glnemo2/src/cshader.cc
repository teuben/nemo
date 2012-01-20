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
#include "cshader.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <glwindow.h>
#include <QFile>
#include <QString>
#include <QTextStream>
using namespace glnemo;
using namespace std;
// ============================================================================
// constructor        
CShader::CShader(std::string _vert_file, std::string _frag_file, bool _v)
{
  vert_file  = _vert_file;
  frag_file  = _frag_file;
  verbose = _v;
}

// ============================================================================
// init      
bool CShader::init()
{
  bool ret=false;
  if (processVertex()) {
    if(processPixel()) {
      if (createProgram()) {
        ret=true;
      }
    }
  }
  return ret;  
}
// ============================================================================
// start
void CShader::start() {
  glUseProgramObjectARB(m_program);
}
// ============================================================================
// stop
void CShader::stop() {
  glUseProgramObjectARB(0);
}
// ============================================================================
// sendUniformXfv
void CShader::sendUniformXfv(const char * s,const int _dim, const int _count, const float * _v)
{
  GLint loc = glGetUniformLocation(m_program, s);
  if (loc == -1) {
    std::cerr << "Error occured when sending \""<<s<<"\" to shader..\n";
    exit(1);
  }
  switch (_dim) {
  case 1: glUniform1fv(loc,_count,_v); break;// 
  case 2: glUniform2fv(loc,_count,_v); break;// 
  case 3: glUniform3fv(loc,_count,_v); break;// 
  case 4: glUniform4fv(loc,_count,_v); break;// 
  default: 
    std::cerr << "CShader::sendUniformXfv unknown dimension ["<<_dim<<"], abort\n";
    std::exit(1);
  }
}
// ============================================================================
// sendUniformf
void CShader::sendUniformf(const char * s,const float _v)
{
  GLint loc = glGetUniformLocation(m_program, s);
  if (loc == -1) {
    std::cerr << "CShader::sendUniformf Error occured when sending \""<<s<<"\" to shader..\n";
    exit(1);
  }
  glUniform1f(loc, _v);  
}
// ============================================================================
// sendUniformi
void CShader::sendUniformi(const char * s,const int _v)
{
  GLint loc = glGetUniformLocation(m_program, s);
  if (loc == -1) {
    std::cerr << "CShader::sendUniformi Error occured when sending \""<<s<<"\" to shader..\n";
    exit(1);
  }
  glUniform1i(loc, _v);  
}
// ============================================================================
// createProgram()      
bool CShader::createProgram()
{
  bool ret=true;
  std::cerr << "Creating program\n";
  m_program = glCreateProgramObjectARB();
  std::cerr << "Attaching vertex shader\n";
  glAttachObjectARB(m_program, m_vertexShader);
  std::cerr << "Attaching pixel shader\n";
  glAttachObjectARB(m_program, m_pixelShader);
  
  // bind attribute
  //glBindAttribLocation(m_program, 100, "a_sprite_size");
  std::cerr << "Linking  program\n";
  glLinkProgramARB(m_program);
  GLWindow::checkGLErrors("link Shader program");
  printLog(m_program,"Linking  program");
  int  link_status;
  glGetProgramiv(m_program, GL_LINK_STATUS, &link_status);
  if(link_status != GL_TRUE) {
    cerr << "Unable to LINK Shader program.....\n";
    exit(1);
  }
  glDeleteShader(m_vertexShader);
  glDeleteShader(m_pixelShader);
  std::cerr << "ending init shader\n";
  return ret;
}
// ============================================================================
// processVertex()      
bool CShader::processVertex()
{
  bool ret=true;
  // process VERTEX
  std::string vert_src = load(vert_file);
  if (vert_src.size()>0) {
    m_vertexShader = glCreateShaderObjectARB(GL_VERTEX_SHADER);
    if (!m_vertexShader) {
      cerr << "Unable to create VERTEX SHADER.....\n";
      exit(1);
    }
    const char * v =  vert_src.c_str();
    glShaderSourceARB(m_vertexShader, 1, &v , NULL);
    std::cerr << "Compiling vertex shader\n";
    GLint compile_status;
    glCompileShaderARB(m_vertexShader);
    GLWindow::checkGLErrors("compile Vertex Shader");
    printLog(m_vertexShader,"Compiling vertex shader");
    glGetShaderiv(m_vertexShader, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE) {
      cerr << "Unable to COMPILE VERTEX SHADER.....\n";
      exit(1);
    }
    
  }
  return ret;
}
// ============================================================================
// processPixel()      
bool CShader::processPixel()
{
  bool ret=true;
  // process PIXEL
  std::string  frag_src = load(frag_file);
  if (frag_src.size()>0) {
    m_pixelShader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER);
    if (!m_pixelShader) {
      cerr << "Unable to create PIXEL SHADER.....\n";
      exit(1);
    }
    const char * v =  frag_src.c_str();
    glShaderSourceARB(m_pixelShader, 1, &v , NULL);
    std::cerr << "Compiling pixel shader\n";
    GLint compile_status;
    glCompileShaderARB(m_pixelShader);
    GLWindow::checkGLErrors("compile Pixel Shader");
    printLog(m_pixelShader,"Compiling pixel shader");
    glGetShaderiv(m_pixelShader, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE) {
      cerr << "Unable to COMPILE PIXEL SHADER.....\n";
      exit(1);
    }
  }
  return ret;
}
// ============================================================================
// load      
std::string CShader::load(std::string filename)
{
  std::ifstream fi;
  std::string src="";
  
  // Open file
  QFile infile(filename.c_str());
  if (!infile.open(QIODevice::ReadOnly | QIODevice::Text)) {
    std::cerr << "CShader::load Unable to open file ["<<filename<<"] for reading, ...\n";
    std::exit(1);
  }
  else {  
    QTextStream in(&infile);
    QString line;
    do {
      line = in.readLine();
      if (!line.isNull()) {
        src = src + "\n" + line.toStdString();
      }
    } while (!line.isNull());
    infile.close();
  }
  return src;
}
// ============================================================================
// printLog      
void CShader::printLog(GLuint obj,std::string s)
{
  int infologLength = 0;
  int maxLength;
  
  if(glIsShader(obj))
    glGetShaderiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
  else
    glGetProgramiv(obj,GL_INFO_LOG_LENGTH,&maxLength);
  
  char infoLog[maxLength];
  
  if (glIsShader(obj))
    glGetShaderInfoLog(obj, maxLength, &infologLength, infoLog);
  else
    glGetProgramInfoLog(obj, maxLength, &infologLength, infoLog);
  
  if (infologLength > 0) {
    std::cerr << "CShader::printLog [" << s<< "]:"<<infoLog<<"\n";
  }
  
}
