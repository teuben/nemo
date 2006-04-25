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

void HostnameSelectForm ::init()
{
    // Form's icon
    //setIcon( QPixmap( glnemo35_xpm ) );
    host_edit_list->setEditable(true);
    host_edit_list->insertStrList( hosts_list);
    //first = true;
    clearWState( WState_Polished ); 
}
void HostnameSelectForm::destroy()
{
}
