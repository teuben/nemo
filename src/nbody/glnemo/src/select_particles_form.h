// ============================================================================
// Copyright Jean-Charles LAMBERT - 2004                                       
// e-mail:   Jean-Charles.Lambert@oamp.fr                                      
// address:  Dynamique des galaxies                                            
//           Laboratoire d'Astrophysique de Marseille                          
//           2, place Le Verrier                                               
//           13248 Marseille Cedex 4, France                                   
//           CNRS U.M.R 6110                                                   
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".        
// ============================================================================

#ifndef SELECTPARTICLESFORM_H
#define SELECTPARTICLESFORM_H

#include <qvariant.h>
#include <qdialog.h>
#include <qlineedit.h>

#include "particles_range.h"

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QLabel;
class QLineEdit;
class QPushButton;

class SelectParticlesForm : public QDialog
{
    Q_OBJECT

public:
    SelectParticlesForm( const char * title, const int nbody,QWidget* parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~SelectParticlesForm();
    const QString getSelectedRange() { 
                return line_edit_select_part->text();};
    
private:    
    QLabel* text_nbody;
    QLabel* text_select_part;
    QLabel* text_nbody_value;
    QLineEdit* line_edit_select_part;
    QPushButton* buttonOk;
    QPushButton* buttonCancel;
    QLabel* text_title;

protected:

protected slots:
    virtual void languageChange();

private slots:
  bool isValidParticlesRange();
  
private:
  int nbody;        
};

#endif // SELECTPARTICLESFORM_H
