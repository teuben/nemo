// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique de galaxies                                             
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#ifndef OPTIONSFORM_H
#define OPTIONSFORM_H

#include <qvariant.h>
#include <qdialog.h>
#include "glbox.h"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QSlider;
class QLabel;
class QGroupBox;
class QCheckBox;

class OptionsForm : public QDialog
{
    Q_OBJECT

public:
    OptionsForm( GLBox * glbox, QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~OptionsForm();

    QSlider* slider_part_size;
    QLabel* text_particle_size;
    QGroupBox* blendingBox;
    QCheckBox* checkBox_depthbuffer;
    QCheckBox* checkBox_blending;
    QGroupBox* gridsBox;
    QCheckBox* checkBox_gridX;
    QCheckBox* checkBox_gridZ;
    QCheckBox* checkBox_gridY;
    QGroupBox* hudBox;
    QGroupBox* polyGroup;
    QCheckBox* polyCheck;
    QCheckBox* hud_display;
    
protected:

protected slots:
    virtual void languageChange();

};

#endif // OPTIONSFORM_H
