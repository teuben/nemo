// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                  
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#ifndef ANIMATION_ENGINE_H
#define ANIMATION_ENGINE_H

#include <qobject.h>
#include <qdatetime.h>
#include <qthread.h>
#include <qimage.h>
#include <qtimer.h>

#include "frame_data.h"
// ----------------------------------------------------------------------------
// -------------------------------- B A S E -----------------------------------
// ----------------------------------------------------------------------------
class AnimationBase: public QObject{
Q_OBJECT
public:
 AnimationBase(FrameDataVector *);
 ~AnimationBase();
 
  bool setActivated(bool);
  bool isActivated() { return activated; };
  int wholeTime(); // return whole time in ms     


signals:
  void uploadToGL(GlobalOptions *, const bool);  // send to GLBOX to update GL  
  void loadNextFrame();                          // request next frame loading  
  void infoStatus(QString);
protected:  
  FrameDataVector * fdv;
  bool activated;
  bool select_options_gui;        // true from GUI, false from record file   
  QTimer my_timer;
protected slots:
  void selectOptions(bool);       // set options from GUI or from record file
  bool isNewTimeFrame(const int);  
};
// ----------------------------------------------------------------------------
// ---------------------------- R E C O R D -----------------------------------
// ----------------------------------------------------------------------------
class AnimationRecord: public AnimationBase {
Q_OBJECT
public:
 AnimationRecord(FrameDataVector *);
 ~AnimationRecord();
 
 int beginFrame();
 int endFrame(GlobalOptions *);
 void reset() { fdv->clear(); current_frame_time=-1; setActivated(false);} ;
 int start();
 int stop();
signals:
  void infoRecord(int,int);

private slots:
  void updateFrameTime(const float);
  void allowRecord() { allow_record = true; };
private:
  QTime t; 
  int elapsed;
  int restart() { return (t.restart());};
  bool allow_record;
  float current_frame_time;

};
// ----------------------------------------------------------------------------
// ---------------------------- R E N D E R -----------------------------------
// ----------------------------------------------------------------------------
class AnimationRender: public AnimationBase{
Q_OBJECT
public:
 AnimationRender(FrameDataVector *);
 ~AnimationRender();
  bool start();
  void stop();
signals:
  void infoRender(int, int);
  void takeScreensgot(QImage &);
  void renderButtonTextSignal(QString);
public slots:
  void runTimeout();
  void reset();
private slots:
  void initOptionsSlot(const QString &, const QString & , const int, const bool);
public:
  void pause();
  int rendering; // 0 - false / 1 - true / 2 - pause           
  

private:
  unsigned int current_frame_render;  
  QString dirname,framename;
  int fps;
  bool pngformat;
  
  // method
  bool saveScreenshot(const QImage);
  bool init();
};
// ----------------------------------------------------------------------------
// ----------------------------- P L A Y --------------------------------------
// ----------------------------------------------------------------------------
class AnimationPlay: public AnimationBase{
Q_OBJECT
public:
  AnimationPlay(FrameDataVector *);
  ~AnimationPlay();  
  bool start();
  bool stop();
  void pause();
signals:
  //void uploadToGL(GlobalOptions *, const bool);  // send to GLBOX to update GL
  void endOfRecord();                            // send to main appli        
  void infoPlay(int,int,int);
  void infoStatus(QString);
public slots:
int runTimeout();            // called from QTimer signal   
void displayFrameIndex(int); // called from slider move     

private slots:
//void slotPlay();


private:
  int display(int);        // display OpenGL Frame      
  
  int playing;             // 0 - false / 1 - true / 2 - pause           
  int current_frame_index; // easy                                       
  int cumul_elapsed;       // cumulated eplased time at the current frame
  bool reset;              // ...                                        
  bool first_frame;        // special check for the first frame          
  QTime elapsed_in_timer;  // time spent                                 
  int whole_time;          // whole time animation                       
  //QTimer my_timer;
};
// ----------------------------------------------------------------------------
// ---------------------------- E N G I N E -----------------------------------
// ----------------------------------------------------------------------------
class AnimationEngine: public QObject {
Q_OBJECT
public:
  AnimationEngine();

  ~AnimationEngine();
    
  AnimationRecord * record;
  AnimationRender * render;
  AnimationPlay   * play;    
  FrameDataVector fdv; 

private slots:
  void playSlot();
  void pauseSlot();
  void stopSlot();
  void recordSlot();
  void renderSlot();
  void resetPrSlot();
  void firstFrameSlot();
  void lastFrameSlot();
  void displayFrameIndexSlot(int);
private:
  enum STATUS {
    NONE, PLAY, RECORD, RENDER
  };
  int status;
  // frame data base 
    
  // timing variables
  QTime t;
  
};
// ============================================================================
#endif
