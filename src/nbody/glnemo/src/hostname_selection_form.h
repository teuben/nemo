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


#ifndef HOSTNAMESELECTIONFORM_H
#define HOSTNAMESELECTIONFORM_H

#include <qvariant.h>
#include <qdialog.h>

class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QPushButton;
class QLabel;
class QLineEdit;

class HostnameSelectionForm : public QDialog
{
    Q_OBJECT

public:
    HostnameSelectionForm(QWidget* parent = 0, const char* name = 0, 
                          bool modal = FALSE, WFlags fl = TRUE );
    ~HostnameSelectionForm();

    QPushButton* buttonOk;
    QPushButton* buttonCancel;
    QLabel* text_hostname;
    QLineEdit* edit_hostname;

protected:
    QGridLayout* layout4;
    QHBoxLayout* layout3;

protected slots:
    virtual void languageChange();

};

#endif // HOSTNAMESELECTIONFORM_H
