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
//                                                                             
// RunningServerForm class definition                                          
//                                                                             
// ============================================================================
#ifndef RUNNINGSERVERFORM_H
#define RUNNINGSERVERFORM_H

#include <qvariant.h>
#include <qdialog.h>
#include "network_data.h"
#include "globjwin.h"


class QVBoxLayout;
class QHBoxLayout;
class QGridLayout;
class QPushButton;
class QLabel;
class QListBox;
class QListBoxItem;
class GLObjectWindow ;

class RunningServerForm : public QDialog
{
    Q_OBJECT

public:
    RunningServerForm( GLObjectWindow * parent = 0, const char* name = 0, bool modal = FALSE, WFlags fl = 0 );
    ~RunningServerForm();
    
    void fillList();
    QPushButton* cancel_button;
    QLabel* server_title;
    QListBox* server_list;
    QString hostname;
    QPushButton* connect_button;

private slots:
    void hostSelected(QListBoxItem* );
    void acceptHostSelected(QListBoxItem* );
    int findHost(QListBoxItem* );
protected:

protected slots:
    virtual void languageChange();

};

#endif // RUNNINGSERVERFORM_H
