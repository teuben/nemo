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
#include "frame_data.h"
#include <qmessagebox.h>
#include <string>
FrameData::FrameData()
{
}
// ============================================================================
// ============================================================================

FrameData::~FrameData()
{
}

// ============================================================================
// Constructor                                                                 
// ============================================================================
IOFrameData::IOFrameData(const QString& filename_, bool is_load_)
{
  filename = filename_;
  is_load  = is_load_;
  if (is_load) {
    in.open(filename,std::ios::in| std::ios::binary);
    if ( ! in.is_open()) {
      QString message = "Unable to open file\n"+filename+"\nfor reading";
      QMessageBox::information(NULL,"Warning",message,"Ok");
    }
  } else {
    out.open(filename,std::ios::out| std::ios::binary);
    if ( ! out.is_open()) {
      QString message = "Unable to open file\n"+filename+"\nfor writing";
      QMessageBox::information(NULL,"Warning",message,"Ok");
    }

  }
}

// ============================================================================
// ============================================================================
IOFrameData::~IOFrameData()
{
}
// ============================================================================
// ============================================================================
int IOFrameData::save( FrameDataVector*fdv)
{
  int status=0;
  if (is_load) return status;
  if ( fdv->size() > 0 ) {                          // is thre something to save ?
    status=1;
    std::string cookie(MAGICFDCOOKIE);
    out.write((char *) cookie.c_str(),cookie.length()); // magic number
    int nframe = fdv->size();                           // #frames
    out.write((char *) &nframe,sizeof(int));            // save #frames
    for (unsigned int i=0; i<fdv->size(); i++) {        // loop on frames
      FrameData * fd = &(*fdv)[i];                      // point on current frame
      out.write((char *) fd,sizeof(FrameData));         // save
    }
  }
  out.close();
  return status;
}
// ============================================================================
// ============================================================================
int IOFrameData::load( FrameDataVector*fdv)
{
  int status=0;
  if (!is_load) return status;
  // clear current frame data vector
  int nframe=0;
  std::string cookie(MAGICFDCOOKIE);
  char * read_cookie = new char[cookie.length()+1];   // space to read cookie
  in.read(read_cookie,cookie.length());              // read cookie
  read_cookie[cookie.length()]='\0';
  std::string ss(read_cookie);                        // char to string
  delete [] read_cookie;                              // free memory
  if ( ss != cookie) {                                // is valid cookie?
    // not a valid cookie
    QString message = "File:"+filename+
                      "\nis not a valid GLNEMO FRAME DATA file";
    QMessageBox::warning(NULL,"Warning",message,"Ok");
  } else {
    fdv->clear();
    status=1;
    in.read((char *) &nframe,sizeof(int));     // read #frames
    for (int i=0; i<nframe; i++) {
      FrameData fd;
      in.read((char *) &fd,sizeof(FrameData));   // read a new frame
      fdv->push_back(fd);                        // store new frame
    }
    //std::cerr << "[" << fdv->size() << "] frames read \n";
    
  }
  in.close();
  return status;
}
// ============================================================================
