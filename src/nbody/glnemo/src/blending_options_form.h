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

#ifndef BLENDINGOPTIONSFORM_H
#define BLENDINGOPTIONSFORM_H

#include <qvariant.h>
#include <qdialog.h>

#include "glbox.h"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QSlider;
class QCheckBox;

class BlendingOptionsForm : public QDialog
{
    Q_OBJECT

public:
  BlendingOptionsForm( GLBox * glbox, QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~BlendingOptionsForm();

    QLabel* text_particle_size;
    QSlider* slider_part_size;
    QCheckBox* checkBox_blending;
    QCheckBox* checkBox_depthbuffer;

protected:

protected slots:
    virtual void languageChange();

};

#endif // BLENDINGOPTIONSFORM_H
