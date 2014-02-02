// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2014
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
#include <QHostAddress>
#include <QCoreApplication>
#include <iostream>

#include "masterserver.h"
#include "client.h"
#include "sendreseau.h"
#include "global.h"

using namespace network;
using namespace std;

// ============================================================================
// Constructor
// ============================================================================
MasterServer::MasterServer(int _port, int _maxConnexion, std::string _sim_name, const falcON::snapshot * _S) {

    while (!listen(QHostAddress::Any, _port)) {
        _port++;
    }
    port = _port;
    maxConnexion = _maxConnexion;
    nbConnexionFree = maxConnexion;
    sim_name = _sim_name;
    S = _S;
    clients.clear();
    std::cerr<<"Server listening on port "<<_port<<"\n";

}

// ============================================================================
// Destructor
// ============================================================================
MasterServer::~MasterServer() {

    close();
    deleteLater();

}

// ============================================================================
// Wait for new connection
// ============================================================================
void MasterServer::startList() {

    while (isListening()) {

       if (clients.length() > (maxConnexion-nbConnexionFree)) {
           //std::cerr << "<<< client lentgh="<<clients.length()<<"\n";
//            for (int i = 0; i < clients.length(); i++) {
//                if (!clients.at(i)->isRunning()) {
//                    delete clients.at(i);
//                    clients.removeAt(i);
//                }
//            }
            int cpt=0;
            while (cpt<clients.length()) {
                if (!clients.at(cpt)->isRunning()) {
                    delete clients.at(cpt);
                    clients.removeAt(cpt);
                    cpt=0;
                } else {
                    cpt++;
                }
            }
            //std::cerr << " >>> client lentgh="<<clients.length()<<"\n";
        }
       //std::cerr << " >>> client lentgh="<<clients.length()<<"\n";

        if (waitForNewConnection(3000)) {
            QTcpSocket *tmpS = nextPendingConnection();
            Client *cl = new Client(tmpS->socketDescriptor(), S, &mutex, &cond, &nbConnexionFree);
            
            if (!lstConnected.contains(tmpS->peerAddress().toString())) {
                std::cerr<<"New client connected on port "<<tmpS->peerPort()<<"    ["<<tmpS->peerAddress().toString().toStdString()<<"]\n";
                lstConnected.insert(0, tmpS->peerAddress().toString());
            }
            //std::cerr<<"Connection\n";
            
            clients.insert(0, cl);
            nbConnexionFree--;
        }
        //std::cerr<<"No connection\n";

    }

}

// ============================================================================
// Say at all children server that there are new data
// ============================================================================
void MasterServer::updateData() {

    cond.wakeAll();

}
