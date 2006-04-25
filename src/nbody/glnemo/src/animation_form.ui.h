// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004-2005                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================
#include <iostream>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qdatetime.h>
// ============================================================================
// AnimationForm::init()
void AnimationForm::init()
{
    //std::cerr << "Here -> AnimationForm::init()\n";
    record=false;
    play=false;
}
// ============================================================================
// AnimationForm::destroy()
void AnimationForm::destroy()
{
}
// ============================================================================
//                                                                     R E C O R D  /  P L A Y      T A B 
// ============================================================================

// AnimationForm::playRecord()
void AnimationForm::playRecord()
{
    emit play_signal();
}

// ============================================================================
// AnimationForm::actionReset()
void AnimationForm::actionReset()
{
    infoRecord(0,0);
        
    // update slider
    timing_slide->setValue(0);
    emit reset_signal();
}

// ============================================================================
// AnimationForm::browseScriptFile()
void AnimationForm::browseScriptFile()
{
    QString fn = QFileDialog::getOpenFileName( edit_filename->text(), QString::null,this);
    if ( !fn.isEmpty() ) {
	edit_filename->setText(fn);
    }
}

// ============================================================================
// AnimationForm::commitScriptFile()
void AnimationForm::commitScriptFile()
{

}

// ============================================================================
// AnimationForm::gotoAnimBegin()
void AnimationForm::gotoAnimBegin()
{

}

// ============================================================================
// AnimationForm::gotoAnimEnd()
void AnimationForm::gotoAnimEnd()
{

}

// ============================================================================
// AnimationForm::playPauseAnim()
void AnimationForm::playPauseAnim()
{
    emit start_play_signal(); 
}	

// ============================================================================
// AnimationForm::stopAnim()
void AnimationForm::stopAnim()
{
   emit stop_play_signal();    
}

// ============================================================================
// AnimationForm::recordAnim()
void AnimationForm::recordAnim()
{
   emit start_record_signal();
}

// ============================================================================
// AnimationForm::sliderAnim()
void AnimationForm::sliderAnim(int frame)
{
   display_frame_index(frame);
}
// ============================================================================
// AnimationForm::infoRecord()
void AnimationForm::infoRecord( int _whole_elapsed, int nframe )
{
    QTime t,n;
    whole_elapsed = _whole_elapsed;
    t=n.addMSecs(_whole_elapsed); 
    frame_duration->setText(t.toString());
    frame_lcd->display(nframe);
}
// ============================================================================
// AnimationForm::infoPlay()
void AnimationForm::infoPlay(int _max_frame, int current_frame, int current_elapsed)
{
    // update slider
    timing_slide->setMinValue(0);
    timing_slide->setMaxValue(_max_frame-1);
    //timing_slide->setValue(current_elapsed);
    timing_slide->setValue(current_frame);
    // update LCD screen
    frame_lcd->display(current_frame);
    // update duration time
    QTime t,n;
    t=n.addMSecs(current_elapsed); 
    frame_duration->setText(t.toString());
}
// ============================================================================
// AnimationForm::statusSlot
// set status
void AnimationForm::statusSlot( QString sstatus)
{
    recp_status->setText(sstatus);
}
// ============================================================================
// AnimationForm::prOptions()
// Playing?rendering Options
void AnimationForm::prOptions()
{
    emit pr_options(gui_options->isChecked());
}
// ============================================================================
//                                                                     R E N D E R I N G     T A B 
// ============================================================================

// ============================================================================
// AnimationForm::selectRenderDir()
bool AnimationForm::selectRenderDir()
{
    bool ret=true;
    QDir render_dir(render_dir_name->text());

    if ( ! render_dir.exists()) {
	int status=QMessageBox::question(this,NULL,
			       "Directory ["+render_dir_name->text()+"] does not exist\n"
			       "Would you like to create it ?","Ok","Cancel");
	if (status == 0) {                 // ok we want to mkdir
	    std::cerr << "PATH = " << render_dir.absPath() << "\n";
	    if ( ! render_dir.mkdir(render_dir.path()))  { // mkdir failed
		QMessageBox::critical(this,NULL,
				      "Unable to create directory ["+render_dir_name->text()+"]","Ok");
		ret = false;
	    } else {
		ret=true;
	    }
	} else {
	    ret=false;
	}
    } 
    return ret;
}
// ============================================================================
// AnimationForm::startRendering()
// !!!!!!!!!!!!
// Start and Pause text MUST be modified by render engine
// !!!!!!!!!!!!!
void AnimationForm::startRendering()
{
    if (selectRenderDir()) {
	emit rendering_options(render_dir_name->text(),
			       render_frame_name->text(),
			       frame_fps->value(),png->isChecked());
	emit start_render_signal();
     }

}
// ============================================================================
// AnimationForm::progressRender()
void AnimationForm::progressRender(int frame_index, int max_frame )
{
    progress_render->setProgress(frame_index,max_frame);
    
}
// ============================================================================
// AnimationForm::resetRendering()
void AnimationForm::resetRendering()
{
    emit reset_render_signal();
    progress_render->setProgress(0);      // reset progress bar
    start_button_render->setText("Start"); // change start rendering button's name
}
// ============================================================================
// AnimationForm::renderButtonTextSlot()
void AnimationForm::renderButtonTextSlot( QString  name)
{
    start_button_render->setText(name);
}
