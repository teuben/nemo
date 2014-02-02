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
#include <QTcpSocket>
#include <QCoreApplication>
#include <iostream>

#include "master_server_thread.h"
#include "masterserver.h"
#include "global.h"

// ============================================================================
// Constructor
// ============================================================================
MasterServerThread::MasterServerThread(std::string _sim_name, int _port, int _max_port, const falcON::snapshot * S) {

    port = _port; //Initialisation du port de connexion
    nbConnexionAutorise = _max_port;
    my_snapshot = S;
    sim_name = _sim_name;
    start();

}

// ============================================================================
// Destructor
// ============================================================================
MasterServerThread::~MasterServerThread() {

    delete server;

}

// ============================================================================
// New data arrived !
// ============================================================================
void MasterServerThread::updateData() {

    server->updateData();

}

// ============================================================================
// Start listening of new connection
// ============================================================================
void MasterServerThread::run() {

    server = new MasterServer(port, nbConnexionAutorise, sim_name, my_snapshot);
    server->startList();

}
