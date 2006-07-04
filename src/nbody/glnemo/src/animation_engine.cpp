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
#include <iostream>
#include "animation_engine.h"
#include <qmessagebox.h>
#include <iomanip>

// ----------------------------------------------------------------------------
// ---------------------------- E N G I N E -----------------------------------
// ----------------------------------------------------------------------------

// ============================================================================
// Constructor                                                                 
AnimationEngine::AnimationEngine()
{
  // allocate memory
  record = new AnimationRecord(&fdv);
  play   = new AnimationPlay(&fdv);
  render = new AnimationRender(&fdv);
  status = AnimationEngine::NONE;
}
// ============================================================================
// AnimationEngine::playSlot()
void AnimationEngine::playSlot()
{
  switch (status) {
  case  NONE :
  case  PLAY : 
    if (fdv.size()) {  // there is frame to be played
      status = PLAY;
      play->start();
    }
    break;
  case  RECORD       :
    record->stop();
    status = NONE;
    if (play->start()) {
      status = PLAY;
    }
    break;
  case  RENDER :
    if (render->rendering == 0) { // not rendering anymore
      if (fdv.size()) {  // there is frame to be played
        status = PLAY;
        play->start();
      }       
    }  
    //render->start();
    break;  
  }
}
// ============================================================================
// AnimationEngine::pauseSlot()
void AnimationEngine::pauseSlot()
{
}
// ============================================================================
// AnimationEngine::stopSlot()
void AnimationEngine::stopSlot()
{
  switch (status) {
  case  NONE :
  case  PLAY : 
    if (fdv.size()) {  // there is
      play->stop();
    }
    break;
  case  RECORD       :
    record->stop();
    break;
  case  RENDER :
    render->stop();
    break;  
  }
  status=NONE;
}
// ============================================================================
// AnimationEngine::recordSlot()
void AnimationEngine::recordSlot()
{
  switch (status) {
  case NONE:
  case RECORD :
      status = RECORD;
      record->start();
      break;
  case  PLAY  : 
      break;
  case  RENDER :
      if (render->rendering == 0) { // not rendering anymore
        status = RECORD;
        record->start();
      }
    break;  
  }
}
// ============================================================================
// AnimationEngine::displayFrameIndexSlot()
void AnimationEngine::displayFrameIndexSlot(int index)
{
  switch (status) {
  case NONE:
  case  PLAY  : 
    play->displayFrameIndex(index);
    break;
  case RECORD :
    break;   
  case  RENDER :
    break;  
  }
}
// ============================================================================
// AnimationEngine::renderSlot()
void AnimationEngine::renderSlot()
{
  switch (status) {
  case  NONE :
  case  PLAY : 
    play->stop();
    status=NONE;
    if (render->start()) {
      status=RENDER;
    }
    break;
  case  RECORD       :
    record->stop();
    status = NONE;
    if (render->start()) {
      status = RENDER;
    }
    break;
  case  RENDER :
      render->start();
    break;  
  }
}
// ============================================================================
// AnimationEngine::resetPrSlot()
void AnimationEngine::resetPrSlot()
{
  stopSlot();
  fdv.clear();  
}
// ============================================================================
// AnimationEngine::firstFrameSlot()
void AnimationEngine::firstFrameSlot()
{
}
// ============================================================================
// AnimationEngine::lastFrameSlot()
void AnimationEngine::lastFrameSlot()
{
}
// ============================================================================
// Destructor                                                                  
AnimationEngine::~AnimationEngine()
{
  fdv.clear();
  delete record;
  delete play;
  delete render;  
}
// ----------------------------------------------------------------------------
// -------------------------------- B A S E -----------------------------------
// ----------------------------------------------------------------------------

// ============================================================================
// Constructor                                                                 
AnimationBase::AnimationBase(FrameDataVector * _fdv)
{
  // get pointer on frame data vector
  fdv = _fdv;
  // initialyse variable             
  activated=false;
  //
  select_options_gui = true;
}

// ============================================================================
// Destructor                                                                  
AnimationBase::~AnimationBase()
{
}
// ============================================================================
// AnimationBase::selectOptions(bool _b)                                       
// set playing/rendering option                                                
void AnimationBase::selectOptions(bool _b)
{
  select_options_gui = _b; 
}

// ============================================================================
// AnimationBase::setActivated()                                               
bool AnimationBase::setActivated(bool b)
{
  return(activated = b); 
}
// ============================================================================
// AnimationBase::wholeTime()                                                  
// return elapsed time stored in FrameDataVector                               
int AnimationBase::wholeTime()
{
  int whole_elapsed=0;
  for (unsigned int i=0; i< fdv->size(); i++) {
    whole_elapsed += (*fdv)[i].elapsed;
  }
  return whole_elapsed;
}
// ============================================================================
// AnimationBase::isNewTimeFrame()                                             
// return true if the time of the current frame is newer                       
bool AnimationBase::isNewTimeFrame(const int current_index)
{
  bool ret=false;
  if (current_index > 0) {
    if ( (*fdv)[current_index].time >      // current time >
         (*fdv)[current_index-1].time ) {  // previous time 
      ret = true;
    } 
  }
  return ret;
}
// ----------------------------------------------------------------------------
// ---------------------------- R E C O R D -----------------------------------
// -------------------------------------------------------------------------

// ============================================================================
// Constructor                                                                 
AnimationRecord::AnimationRecord(FrameDataVector * _fdv):AnimationBase(_fdv)
{
  current_frame_time = -1;
}

// ============================================================================
// Destructor                                                                  
AnimationRecord::~AnimationRecord()
{
}
// ============================================================================
// AnimationRecord::beginFrame()                                               
// called from glbox::paintGL() at the beginning                               
int AnimationRecord::beginFrame()
{ if (activated && allow_record) {
    return(elapsed = t.elapsed()); 
  } else {
    return 0;
  }
};
// ============================================================================
// AnimationRecord::endFrame()                                                 
// Function called at the end of paintGL() function. Allow storing of the      
// current frame options.                                                      
int AnimationRecord::endFrame(GlobalOptions * go)
{
  if (activated && allow_record) {
    restart(); // restart time
    
    FrameData fd;
    fd.store_options = *go;
    fd.time          = current_frame_time;
    fdv->push_back(fd); // record current options
 
    int size = fdv->size();
    if (size > 1) { // from frame #2
      (*fdv)[size-2].elapsed=elapsed;
      if (size > 2 ) {   // #frame > 2
        (*fdv)[size-2].cumul_elapsed = elapsed + (*fdv)[size-3].cumul_elapsed;
      } else {            // #frame == 2
        (*fdv)[size-2].cumul_elapsed = elapsed;
      }
      //std::cerr << "elapsed cumuled = ["<<(*fdv)[size-2].cumul_elapsed<<"]\n";
    }
      
    //int whole_elapsed=wholeTime();
    int whole_elapsed=(*fdv)[size-2].cumul_elapsed;
    //std::cerr << "elapsed computed = ["<<whole_elapsed<<"\n";
    emit infoRecord(whole_elapsed,size-1);
  }
  return 1;
}
// ============================================================================
// AnimationRecord::start()                                                     
int AnimationRecord::start()
{ 
  int status=0;
  if (!activated) {
    allow_record=true; 
    setActivated(true);
    restart();
    status=1;
  }
  emit infoStatus("Recording");
  return(status);
};
// ============================================================================
// AnimationRecord::stop()                                                     
// Record the elapsed time for the last frame                                  
int AnimationRecord::stop()
{
  setActivated(false);
  int size = fdv->size();
  (*fdv)[size-1].elapsed = t.elapsed();
  if (size > 1) {
    (*fdv)[size-1].cumul_elapsed = (*fdv)[size-1].elapsed + (*fdv)[size-2].cumul_elapsed;
  } else {
    (*fdv)[size-1].cumul_elapsed = (*fdv)[size-1].elapsed ;
  }
  emit infoStatus("Stop recording");
  return 1;
}
// ============================================================================
// AnimationRecord::stop()                                                     
// slot called by newTime() signal to inform a new frame arrival               
void AnimationRecord::updateFrameTime(const float new_time)
{
  if (activated) { 
    allow_record = false;      // forbid recording after a new time step
    current_frame_time = new_time;  // update frame time
  }
}
// ----------------------------------------------------------------------------
// ----------------------------- P L A Y --------------------------------------
// ----------------------------------------------------------------------------

// ============================================================================
// Constructor                                                                 
AnimationPlay::AnimationPlay(FrameDataVector * _fdv):AnimationBase(_fdv)
{
  playing = 0;
  current_frame_index = 0;
  reset=true;
  first_frame=true;
  cumul_elapsed=0;
  connect(&my_timer,SIGNAL(timeout()),this,SLOT(runTimeout()));
}

// ============================================================================
// Destructor                                                                  
AnimationPlay::~AnimationPlay()
{
}
// ============================================================================
// AnimationPlay::runTimeout()                                                 
int AnimationPlay::runTimeout()
{
  
  switch (playing) {
    case 0: // stop re
      reset=true;
      first_frame=true;
      return 0;
      break;  
      
    case 1: // start
      if (reset) { // reset frame timer
        reset=false;
        elapsed_in_timer.restart();
      }
      //
      if (first_frame) {
        first_frame=false;
        display(0);                
      }
      else {
        if (current_frame_index+1 < (int) fdv->size()) {  // still exist frame to display
          if ((elapsed_in_timer.elapsed() >               // time to display             
              (*fdv)[current_frame_index].elapsed)) {             
            // display next frame
            current_frame_index++;                              // current new frame index    
            cumul_elapsed+=(*fdv)[current_frame_index].elapsed; // elapsed time since begining
            //emit infoPlay(whole_time,current_frame_index,     // update information regard- 
            emit infoPlay((int) fdv->size(),current_frame_index,// update information regard- 
                          cumul_elapsed);                       // -ing play action.          
            display(current_frame_index);                       // update display             
            //std::cerr << "Current frame time = " << (*fdv)[current_frame_index]. time << "\n";
          }  
        }
        else { // end of record
          stop(); // raz
          //emit endOfRecord();
        }
      }
      break;
      
    case 2: // pause
      break;
  }          

  return 1;
}
// ============================================================================
// AnimationPlay::display()                                                    
int AnimationPlay::display(int index)
{
  // get Options from currect frame
  GlobalOptions store_options = (*fdv)[index].store_options;
  
  if (isNewTimeFrame(index)) {  
    emit loadNextFrame();
  }
  else {
    // send to GLBox
    emit uploadToGL(&store_options,select_options_gui);  
  }
  // reset timer
  elapsed_in_timer.restart();
  return 1;
}
// ============================================================================
// AnimationPlay::displayFrameIndex(int)                                       
void AnimationPlay::displayFrameIndex(int index)
{
  if (fdv->size() && ((playing==2)||(playing==0))) { // frame && (pause|stop) playing
    display(index); 
    emit infoPlay((int) fdv->size(),index,       // update information regard- 
                  (*fdv)[index].cumul_elapsed);  // -ing play action.          
    current_frame_index = index;
  }
}
// ============================================================================
// AnimationPlay::start()                                                      
// return true if it exist recorded frame                                      
bool AnimationPlay::start()
{
  bool status=true;
  if (playing == 2 ) {  // from pause mode
    playing = 1;       // to play mode   
    emit infoStatus("Playing");
  } 
  else {
    if (fdv->size()) {
      if (playing == 1 ) { // from play mode
         playing = 2;      // to pause mode 
        emit infoStatus("Pause playing");
      }
      else {               // stop mode          
        playing=1;         // switch to play mode
        whole_time = wholeTime();
        status = true;
        my_timer.start(1);
        emit infoStatus("Playing");
      }
    } 
    else { // no recorded frames
      emit infoStatus("None");
      status = false;
    }
  }
  return status; 
}
// ============================================================================
// AnimationPlay::pause()                                                      
void AnimationPlay::pause()
{
  playing=2;
}
// ============================================================================
// AnimationPlay::stop()                                                       
bool AnimationPlay::stop()
{
  playing=0;
  current_frame_index = 0;
  reset=true;
  first_frame=true;
  cumul_elapsed=0;
  my_timer.stop();
  emit infoStatus("None");
  return true;
  
}

// ----------------------------------------------------------------------------
// ---------------------------- R E N D E R -----------------------------------
// ----------------------------------------------------------------------------

// ============================================================================
// Constructor                                                                 
AnimationRender::AnimationRender(FrameDataVector * _fdv):AnimationBase(_fdv)
{
  rendering = 0;
  connect(&my_timer,SIGNAL(timeout()),this,SLOT(runTimeout()));
}
// ============================================================================
// Destructor                                                                  
AnimationRender::~AnimationRender()
{
}
// ============================================================================
// AnimationRender::start()                                                    
bool AnimationRender::start()
{
  bool status=true;
  if (rendering == 2 ) {  // from pause mode
    rendering = 1;       // to play mode   
    emit renderButtonTextSignal("Pause");
  } 
  else {
    if (init()) {
      if (rendering == 1 ) { // from play mode
         rendering = 2;      // to pause mode 
        emit renderButtonTextSignal("Continue");
      }
      else {               // stop mode          
        rendering=1;       // switch to play mode
        status = true;
        my_timer.start(1);
        emit renderButtonTextSignal("Start");
      }
    } 
    else { // no recorded frames
      my_timer.stop();
      emit renderButtonTextSignal("None");
      status = false;
    }
  }
  return status; 
}
// ============================================================================
// AnimationRender::stop()                                                     
void AnimationRender::stop()
{
  my_timer.stop();
}
// ============================================================================
// AnimationRender::reset()                                                    
void AnimationRender::reset()
{ 
  stop();
  rendering = 0; // no rendering
  activated = false; 
  current_frame_render=0;
};  
// ============================================================================
// AnimationRender::initOptionsSlot()                                          
// init rendering animation options                                            
void AnimationRender::initOptionsSlot(const QString & _dirname, const QString & _framename,
      const int _fps, const bool _pngformat)
{
  dirname    = _dirname;
  framename  = _framename;
  fps        = _fps;
  pngformat  = _pngformat;
}      
// ============================================================================
// AnimationRender::init()                                                     
// init rendering animation                                                    
bool AnimationRender::init()
{
  bool ret;
  if (fdv->size() == 0) {  // there is no frame
    QMessageBox::warning(NULL,NULL,
      "Can't start rendering, because there are no animation recorded yet");
      ret=false;
  } 
  else {
    if ( ! isActivated()) {      // if isActivated means a restart rendering
      current_frame_render = 0;
    }
    activated=true;  
    ret=true;
  }
  return ret;
}            
// ============================================================================
// AnimationRender::runTimeout()                                               
// init rendering animation                                                    
void AnimationRender::runTimeout()
{
  if (current_frame_render < fdv->size() &&  (rendering==1)) {
    int nframe = (*fdv)[current_frame_render].elapsed * fps / 1000;
    nframe = (nframe==0?1:nframe);
    std::setfill(0);
    //std::cerr << "frame #"<< std::setw(5) << current_frame_render<<" frame = "<<std::setw(3) << nframe << "\n";
    // ------ render now -----------
    // get Options from currect frame
    GlobalOptions store_options = (*fdv)[current_frame_render].store_options;    
    if (isNewTimeFrame(current_frame_render)) {
      emit loadNextFrame();
    } 
    else {
      emit uploadToGL(&store_options,select_options_gui);           // send to GLBox  
    }
    QImage img;
    emit takeScreensgot(img);                                     // screenshot     
    saveScreenshot(img);
    current_frame_render++;                                       // next frame     
    emit infoRender((int)current_frame_render,(int)fdv->size());  // progess bar    
  }
  if ( current_frame_render == fdv->size() ) {
      my_timer.stop();
      emit renderButtonTextSignal("None");
      
  }
}
// ============================================================================
// AnimationRender::saveScreenshot()                                           
// save screenshot                                                             
bool AnimationRender::saveScreenshot(const QImage img)
{
  QString indx = QString("%1").arg(current_frame_render);
  indx = indx.rightJustify(6,'0');
  QString shot_name=QString("%1/%2.%3.png").arg(dirname).arg(framename).arg(indx);
  
  std::cerr << "Saving frame [" << shot_name << "]\n";
  
  if  ( !img.save( shot_name, "PNG" ) ) { // save the image
    QMessageBox::warning( NULL, 
      "AnimationRender::saveScreenshot failed", 
      "Error saving file" );
    return false;
  }    
  return true;
}
// ============================================================================


