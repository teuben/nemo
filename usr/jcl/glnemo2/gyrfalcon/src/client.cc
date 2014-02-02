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
#include <QTcpSocket>
#include <QHostAddress>
#include <iostream>
#include <cstdlib>
#include <string>
#include "client.h"
#include "global.h"
#include "sendreseau.h"
#include "assert.h"

using namespace std;
using namespace network;

// ============================================================================
// Constructor
// ============================================================================
Client::Client(int _socketDescriptor, const falcON::snapshot *_S, QMutex *_mutex, QWaitCondition *_condition, int *_nbConnectfree) {

    socketDescriptor = _socketDescriptor;
    S = _S;
    mutex = _mutex;
    condition = _condition;

    interval = new QList<int>();

    tablePart = NULL;
    nbBody = S->N_bodies();

    nbConnectFree = _nbConnectfree;
    Freeplace = false;
    initialiseMessage();

    start();

}

// ============================================================================
// Destructor
// ============================================================================
Client::~Client() {

    if (socket->state() == QAbstractSocket::ConnectedState) {
        socket->disconnectFromHost();
    }

    if (!Freeplace) {
        nbConnectFree++;
        Freeplace = true;
    }

    delete [] tablePart;

    interval->clear();
    delete interval;

    delete socket;

}

// ============================================================================
// Initialize for read new data
// ============================================================================
void Client::initialiseMessage() {

    tag = -1;
    tailleMess = -1;
    message = "empty";

}

// ============================================================================
// Wait client for new ask
// ============================================================================
void Client::run() {

    socket = new QTcpSocket();
    socket->setSocketDescriptor(socketDescriptor);

    r.setSocket(socket);
    s.setSocket(socket);

    vitesse = false;

    s.sendData(sizeof(qint16), (char *) &Global::HELLO);

    s.sendData(Global::IDENTSERVEUR.length(), (char *) Global::IDENTSERVEUR.toStdString().c_str());

    while (socket->state() == QAbstractSocket::ConnectedState) {
        readTag();
    }

    if (!Freeplace) {
        *nbConnectFree = *nbConnectFree+1;
        Freeplace = true;
    }

}

// ============================================================================
// Read what type of asking client do
// ============================================================================
void Client::readTag() {

    initialiseMessage();

    r.getData(sizeof(qint16), (char *) &tag);

    if (tag == Global::HELLO) {
        bonjour();
    }
    else if (tag == Global::NBODY) {
        getBody();
    }
    else if (tag == Global::TIME) {
        getTime();
    }
    else if (tag == Global::SELECTBODY) {
        setSelectBody();
    }
    else if (tag == Global::SIZEARRAY) {
        getSizeArray();
    }
    else if (tag == Global::SELECTVEL) {
        vitesse = true;
    }
    else {
        deconnectClient();
    }

}

// ============================================================================
// First connection : Hello !
// ============================================================================
void Client::bonjour() {

    char mess[Global::IDENTCLIENT.size()+1];
    r.getData(Global::IDENTCLIENT.size(), mess);
    mess[Global::IDENTCLIENT.size()]='\0';
    message = mess;

    if (message != Global::IDENTCLIENT) {
        initialiseMessage();
        deconnectClient();
    }

}

// ============================================================================
// Send nb body
// ============================================================================
void Client::getBody() {

    s.sendData(sizeof(int), (char *) &nbBody);

}

// ============================================================================
// Send time of simulation
// ============================================================================
void Client::getTime() {

    float tim = (float) S->time();
    s.sendData(sizeof(float), (char *) &tim);

}

// ============================================================================
// Prepare array to send client
// ============================================================================
void Client::setSelectBody() {

    r.getData(sizeof(int), (char *) &tailleMess);

    char mess[tailleMess+1];
    r.getData(tailleMess, mess);
    mess[tailleMess]='\0';
    select_part = mess;

    NbBodySelect = setInterval(&select_part);

    tablePart = new float[NbBodySelect*3];

}

// ============================================================================
// Add selected particule
// ============================================================================
int Client::addToInterval(int bodyStart, int bodyStop) {

    int i;
    int nb = 0;
    for (i = bodyStart; i <= bodyStop; i++) {
        nb++;
        interval->insert(0, i);
    }

    return nb;

}

// ============================================================================
// Search selected particule
// ============================================================================
int Client::setInterval(std::string *_interval) {

    std::string tmpInterval = *_interval;
    int nbPart = 0;
    std::string inter;

    if (tmpInterval.find("all") != std::string::npos) {
        nbPart += addToInterval(0, S->N_bodies()-1);
    }
    else {
        std::string::size_type stTmp = tmpInterval.find(",");
        while (stTmp != std::string::npos) {
            inter = tmpInterval.substr(0, stTmp);

            std::string::size_type stTmp2 = inter.find(":");
            nbPart += addToInterval(atoi((inter.substr(0, stTmp2)).c_str()), atoi((inter.substr(stTmp2+1)).c_str()));

            tmpInterval = tmpInterval.substr(stTmp+1);
            stTmp = tmpInterval.find(",");
        }

        inter = tmpInterval.substr(0, stTmp);
        std::string::size_type stTmp2 = inter.find(":");
        nbPart += addToInterval(atoi((inter.substr(0, stTmp2)).c_str()), atoi((inter.substr(stTmp2+1)).c_str()));
    }

    return nbPart;

}

// ============================================================================
// Prepare and send the size of array
// ============================================================================
void Client::getSizeArray() {

    mutex->lock();
    condition->wait(mutex);
    mutex->unlock();

    int nbody_out = 0;

    for (int i = 0; i < NbBodySelect; i++) {
        assert(i < interval->length());
        int p_index = interval->at(i);
        assert(p_index < nbBody);
        if (p_index < nbBody) {
            falcON::body b_current = S->bodyNo(p_index);
            falcON::vect x;
            if ( vitesse ) {
                x = b_current.vel();
            }
            else {
                x = b_current.pos();
            }
            assert(nbody_out<NbBodySelect);
            tablePart[nbody_out*3     ] = x[0];
            tablePart[nbody_out*3+1] = x[1];
            tablePart[nbody_out*3+2] = x[2];
            nbody_out++;
        }
    }
    assert(nbody_out==NbBodySelect);

    int tailleTotal = NbBodySelect*3;

    s.sendData(sizeof(int), (char *) &tailleTotal);
    r.getData(sizeof(qint16), (char *) &tag);
    if(tag == Global::POS) {
        getPositions();
    }
    else {
        exit(1);
    }

}

// ============================================================================
// Send the array
// ============================================================================
void Client::getPositions() {

    int tailleTotal = NbBodySelect*3;
    s.sendData(sizeof(float)*tailleTotal, (char *) &tablePart[0]);

    initialiseMessage();
    r.getData(sizeof(qint16), (char *) &tag);
    if (tag == Global::ENDMESSAGES) {
        s.sendData(sizeof(qint16), (char *) &Global::ENDMESSAGES);
        if (vitesse) { vitesse = false; }
    }
    else {
        exit(1);
    }

}

// ============================================================================
// Disconnect client
// ============================================================================
void Client::deconnectClient() {

    if (!Freeplace) {
        *nbConnectFree = *nbConnectFree+1;
        Freeplace = true;
    }

    socket->disconnectFromHost();

}
