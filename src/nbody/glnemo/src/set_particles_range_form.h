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
//                                                                             
// Definition of the form which set the particles range & color                
//                                                                             
//                                                                             
// ============================================================================
#ifndef SET_PARTICLES_RANGE_FORM_H
#define SET_PARTICLES_RANGE_FORM_H

#include "virtual_particles_select.h"
#include "glbox.h"
#include <qdialog.h>

class QHBoxLayout;
class QPushButton;
class QTable;
class QVBoxLayout;

class SetParticlesRangeForm: public QDialog
{
    Q_OBJECT
      public:
  
  SetParticlesRangeForm(const GLBox * glbox,ParticlesSelectVector * , 
  			const int n_body,const float *pos,
			QWidget *parent = 0, const char *name = "set particles range and color",
			bool modal = FALSE, WFlags f = 0 );
  ~SetParticlesRangeForm();
public:
  void updateData(ParticlesSelectVector * , 
  		  const int n_body,const float *pos);
public slots:
    void setColor();
    void setColor( int row, int col );
    void currentChanged( int row, int col );
    void valueChanged( int row, int col );
 signals:
  void applyData(const int *, const float *, ParticlesSelectVector * );

protected slots:
    void accept();
    void apply();
private:
    QTable *table;
    QPushButton *colorPushButton;
    QPushButton *okPushButton;
    QPushButton *applyPushButton;
    QPushButton *cancelPushButton;

    // nbody
    int nbody;
    // pos;
    const float * pos;

    void fillForm();
protected:
    QVBoxLayout *tableButtonBox;
    QHBoxLayout *buttonBox;

private:
    ParticlesSelectVector * my_psv;
};
#endif
// ============================================================================
