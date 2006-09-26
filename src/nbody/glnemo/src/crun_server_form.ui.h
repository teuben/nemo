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
#include "network_data.h"
// ============================================================================
void CListRunningServerForm::init()
{
    updateList(); // look for simulation servers
}
// ============================================================================
void CListRunningServerForm::destroy()
{
}

// ============================================================================
// CListRunningServerForm::updateList()
// Scan the network to find simulation server
void CListRunningServerForm::updateList()
{
  int i=0;        
  static const char * p = hosts_list[i];  
  
  p = hosts_list[i];
  while (p) {   
    // Instantiate a new NetworkData object
    NetworkData * new_virtual_data=
	new NetworkData(p);
    // check connexion
    if (!new_virtual_data->isConnected()) { // not connected ?
    } 
    else { // successfull connexion
  
      // Get NBODY
      int nbody   = new_virtual_data->getNbody();
      // set hostname and nbody
      QString qshost    = QString( "%1" ).arg(QString(p));
      QString qsnbody = QString( "%1" ).arg(nbody);
      // create new item
      QListViewItem * item = new QListViewItem(listRunServer);
      item->setText(0,qshost);
      item->setText(1,qsnbody);      
      // insert new item
      listRunServer->insertItem(item);
    }
    // delet virtual object
    delete new_virtual_data;
    i++;
    p = hosts_list[i]; 
  }    
}
// ============================================================================
// CListRunningServerForm::hostSelected
// record selected hostname
bool CListRunningServerForm::hostSelected( QListViewItem * item )
{
    if (item) {
	// get selected hostname
	hostname = item->text(0);
	//std::cerr << "CListRunningServerForm::hostSelected: " << item->text(0) << "\n";
	return true;
    } else {
	hostname = "";
	return false;
    }
}

// ============================================================================
// CListRunningServerForm::acceptHostSelected
// check if a host has been selected: accept() if yes
void CListRunningServerForm::acceptHostSelected( QListViewItem * item )
{
    // select hostname
    if (hostSelected( item) ) {
	accept();	
    }
}
// ============================================================================
// CListRunningServerForm::refreshList()
// refresh server list
void CListRunningServerForm::refreshList()
{
    listRunServer->clear();
    updateList();
}
// ============================================================================
