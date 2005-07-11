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
// ============================================================================
#include <string.h>
#define _vectmath_h // put this statement to avoid conflict with C++ vector class
#include <nemo.h>
#include <qmessagebox.h>
#include <iostream>
// ============================================================================
void CSelectNbodyForm::init()
{
}
// ============================================================================
void CSelectNbodyForm::destroy()
{
}
// ============================================================================
void CSelectNbodyForm::setData( QString type, QString host, int _nbody )
{
    // set Network Label
    networkDataLabel->setText(type);
    // set hostname
    hostLabel->setText(host);
    // set nbody
    nbodyValueLabel->setText(QString("%1").arg(_nbody));
    nbody = _nbody;
		
}
// ============================================================================
// CSelectNbodyForm::isValidParticlesRange()
// Check if particles range entered is valid
bool CSelectNbodyForm::isValidParticlesRange()
{
  bool status=true;
  QString select_string=range_edit->text();
  
  if (  select_string!="all") {
    int * int_array = new int[nbody];
    const char * s = select_string;
    int npart = nemoinpi(const_cast<char*> (s) , int_array, nbody);
    if (npart <=0 ) {
      status = FALSE;
      std::cerr << "nemoinpi = [" << select_string << "] npart = "<<npart
	      <<"and nbody["   << nbody          <<"\n";
    }
    delete [] int_array;  // useless anymore
  }
  
  if (status) {
    accept();
  } 
  else {
    QString message="Particles selection string misformated";
    QMessageBox::information( this,"Warning",message,"Ok");
  }
  return true;
}
// ============================================================================
// CSelectNbodyForm::getSelectedRange()
const QString CSelectNbodyForm::getSelectedRange()
{
    return range_edit->text();
}
// ============================================================================
