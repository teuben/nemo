// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2015
//           Yannick Dalbin
// e-mail:   Jean-Charles.Lambert@lam.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           Pôle de l'Etoile, site de Château-Gombert
//           38, rue Frédéric Joliot-Curie
//           13388 Marseille cedex 13 France
//           CNRS U.M.R 7326
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".
// ============================================================================
#include <QtGui>
#include <sstream>
#include <QThread>
#include <QMessageBox>
#include "snapshotnetwork.h"
#include "serveur.h"

#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
Q_PLUGIN_METADATA(IID "fr.glnemo2.networkPlugin")
#endif

namespace glnemo {
  
  // ============================================================================
  // Constructor
  //============================================================================
  SnapshotNetwork::SnapshotNetwork() : SnapshotInterface() {
    
    valid=false;
    part_data = NULL;
    part_data = new ParticlesData();
    interface_type    ="Network";
    connected = false;
    
  }
  
  // ============================================================================
  // Destructor
  //============================================================================
  SnapshotNetwork::~SnapshotNetwork() {
    
    delete cl;
    
  }
  
  // ============================================================================
  // Create a new object of SnapshotNetwork and return it
  //============================================================================
  SnapshotInterface * SnapshotNetwork::newObject(const std::string _adresseIP, const int x) {
    
    port = x;
    filename = _adresseIP;
    obj      = new SnapshotNetwork();
    obj->setFileName(_adresseIP);
    obj->setPort(x);
    return obj;
    
  }
  
  // ============================================================================
  // Action when connection with the server is lost...
  //============================================================================
  void SnapshotNetwork::connectFaild(QString message) {
    
    end_of_data=true;
    valid = false;
    QMessageBox::critical(NULL, tr("Logon Failure"), message);
    
  }
  
  // ============================================================================
  // Action when the server is full
  //============================================================================
  void SnapshotNetwork::serveurFull() {
    
    end_of_data=true;
    valid = false;
    QMessageBox::critical(NULL, tr("Can not Connect"), tr("The server is full. The connection is impossible."));
    
  }
  
  // ============================================================================
  // Verify if the connection is ready and good
  //============================================================================
  bool SnapshotNetwork::isValidData() {
    
    valid=false;
    
    if(!connected) {
      cl = new network::Serveur();
      QObject::connect(cl, SIGNAL(connexionFaild(QString)), this, SLOT(connectFaild(QString)));
      QObject::connect(cl, SIGNAL(serveurFull()), this, SLOT(serveurFull()));
      
      cl->setConnexion(filename, port);
      
      if (cl->isBeenConnected(60000)) {
        nbody_first = cl->getNbBody();
        time_first = cl->getTime();
        valid=true;
      }
      else {
        valid=false;
      }
      cl->deco();
      delete cl;
      cl = NULL;
      
      connected = true;
    }
    else {
      valid=true;
    }
    return valid;
    
  }
  
  // ============================================================================
  // getSnapshotRange()
  //============================================================================
  ComponentRangeVector * SnapshotNetwork::getSnapshotRange() {
    
    crv.clear();
    if (valid) {
      ComponentRange * cr = new ComponentRange();
      cr->setData(0,nbody_first-1);
      cr->setType("all");
      crv.push_back(*cr);
      //ComponentRange::list(&crv);
      delete cr;
      if (first) {
        first       = false;
        crv_first   = crv;
      }
    }
    return &crv;
    
  }
  
  // ============================================================================
  // initLoading()
  //============================================================================
  int SnapshotNetwork::initLoading(GlobalOptions * so) {
    go = so;
    load_vel = so->vel_req;
    select_time = so->select_time;  
    return 1;
  }
  
  // ============================================================================
  // Next frame
  //============================================================================
  int SnapshotNetwork::nextFrame(const int * index_tab, const int nsel) {
    
    //bool retur = 0;
    load_vel = go->vel_req;
    if (valid) {
      //Need new connection cause the new thread of nextFrame
      cl = new network::Serveur();
      QObject::connect(cl, SIGNAL(connexionFaild(QString)), this, SLOT(connectFaild(QString)));
      QObject::connect(cl, SIGNAL(serveurFull()), this, SLOT(serveurFull()));
      cl->setConnexion(filename, port);
      if (cl->isBeenConnected(60000)) {
        cl->getNbBody();
        cl->getTime();
        if (cl->readPart(part_data,index_tab,nsel, load_vel, getSelectPart())) {
          part_data->computeVelNorm();
          if (part_data->rho) {
            part_data->rho->computeMinMax();
          }
          end_of_data=false; // only one frame from an ftm snapshot
          //retur = 1;
        }
        else {
          end_of_data=true;
        }
      }
      else {
        end_of_data=true;
        valid = false;
      }
      cl->deco();
      delete cl;
      cl = NULL;
    }
    return 1;
    
  }
  
  // ============================================================================
  // close
  //============================================================================
  int SnapshotNetwork::close() {
    
    int status=0;
    if (valid) {
      status = 1;
      end_of_data = false;
      valid = false; // added 2009 June 19th
    }
    return status;
    
  }
  
  QString SnapshotNetwork::endOfDataMessage() {
    
    QString message=tr("Snapshot [")+QString(filename.c_str())+tr("] end of snapshot reached!");
    return message;
    
  }
  
}

// You have to export outside of the namespace "glnemo"
// BUT you have to specify the namespace in the export:
// ==> glnemo::SnapshotNetwork

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
Q_EXPORT_PLUGIN2(networkplugin, glnemo::SnapshotNetwork);
#endif
