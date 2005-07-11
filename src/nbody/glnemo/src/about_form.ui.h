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
#include "images/glnemo35.xpm"
// ============================================================================
// AboutForm::init()
// Automatically called by constructor
void AboutForm::init()
{
    // Form's icon
    setIcon( QPixmap( glnemo35_xpm ) );
    // glnemo's pixmap
    pixmapLabel3->setPixmap( QPixmap( glnemo35_xpm ) );
}
// ============================================================================
// AboutForm::destroy()
// Automatically called by destructor
void AboutForm::destroy()
{
}
// ============================================================================
