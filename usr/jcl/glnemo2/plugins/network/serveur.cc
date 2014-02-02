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
#include <QTextStream>

#if QT_VERSION >= QT_VERSION_CHECK(5, 0, 0)
#include <QtWidgets/QWidget>
#else // QT4
#include <QMessageBox>
#endif
#include <iostream>

#include "serveur.h"
#include "readreseau.h"
#include "sendreseau.h"
#include "global.h"

namespace network {

    // ============================================================================
    // Constructor
    // ============================================================================
    Serveur::Serveur() {

        socket = new QTcpSocket();
        QObject::connect(socket, SIGNAL(connected()), this, SLOT(connexionSuccefull()));
        QObject::connect(socket, SIGNAL(error(QAbstractSocket::SocketError)), this, SLOT(connexionUnsuccefull(QAbstractSocket::SocketError)));

        time = 0.0;
        nbbody = 0;
        
        r.setSocket(socket);
        s.setSocket(socket);

    }

    // ============================================================================
    // Destructor
    // ============================================================================
    Serveur::~Serveur() {

        deco();
        delete socket;

    }

    // ============================================================================
    // Initialize for read new data
    // ============================================================================
    void Serveur::initialiseReception() {

        tag = -1;
        tailleMess = -1;
        message = "empty";

    }

    // ============================================================================
    // Connect to the server
    // ============================================================================
    void Serveur::setConnexion(std::string ip, int port) {

        socket->connectToHost(QString::fromStdString(ip), port);

    }

    // ============================================================================
    // Wait during new connection
    // ============================================================================
    bool Serveur::isBeenConnected(int mill) {

        if(socket->waitForConnected(mill)) {
            return true;
        }
        else {
            return false;
        }

    }

    // ============================================================================
    // Lanch the hello procedure to confirm the correct connection with the server
    // ============================================================================
    void Serveur::connexionSuccefull() {

        r.getData(sizeof(qint16), (char *) &tag);

        if (tag == Global::HELLO) {

            char mess[Global::IDENTSERVEUR.size()+1];
            r.getData(Global::IDENTSERVEUR.size(), mess);
            mess[Global::IDENTSERVEUR.size()]='\0';
            message = mess;

            if (message == Global::IDENTSERVEUR) {
                s.sendData(sizeof(qint16), (char *) &Global::HELLO);
                s.sendData(Global::IDENTCLIENT.length(), (char *) Global::IDENTCLIENT.toStdString().c_str());
                initialiseReception();
            }
            else {
                std::cerr<<"Wrong message : ["<<message.toStdString()<<"].... Good bye !\n";
                initialiseReception();
                deco();
            }
        }
        else if (tag == Global::SERVEURFULL) {
            emit serveurFull();
            initialiseReception();
            deco();
        }
        else {
            std::cerr<<"Wrong tag. Good bye !\n";
            initialiseReception();
            deco();
        }

        emit connectSuccess();

    }

    // ============================================================================
    // Get the number of total particules that the server have
    // ============================================================================
    int Serveur::getNbBody() {

        s.sendData(sizeof(qint16), (char *) &Global::NBODY);

        r.getData(sizeof(int), (char *) &nbbody);
        initialiseReception();

        return nbbody;

    }

    // ============================================================================
    // Get the time of the simulation
    // ============================================================================
    float Serveur::getTime() {

        s.sendData(sizeof(qint16), (char *) &Global::TIME);

        r.getData(sizeof(float), (char *) &time);
        initialiseReception();

        return time;

    }

    // ============================================================================
    // Read the snapshot
    // ============================================================================
    int Serveur::readPart(glnemo::ParticlesData * part_data, const int *index, const int nsel, bool  load_vel, const std::string _selectPart) {
      if (index) {;} // remove compiler warning
        initialiseReception();

        if ( *part_data->nbody < nsel) {
            *part_data->nbody = nsel;
            if (part_data->pos) {
                delete [] part_data->pos;
            }
            part_data->pos = new float[nsel*3];
            if (load_vel) {
                if (part_data->vel) {
                    delete [] part_data->vel;
                }
                part_data->vel = new float[nsel*3];
            }
        }

        s.sendData(sizeof(qint16), (char *) &Global::SELECTBODY);

        std::string tmpSlectPart = _selectPart;
        int tmpTailleMess = tmpSlectPart.size();
        s.sendData(sizeof(int), (char *) &tmpTailleMess);

        s.sendData(tmpTailleMess, (char *) tmpSlectPart.c_str());

        initialiseReception();

        s.sendData(sizeof(qint16), (char *) &Global::TIME);
        r.getData(sizeof(float), (char *) part_data->timu);

        s.sendData(sizeof(qint16), (char *) &Global::SIZEARRAY);
        r.getData(sizeof(int), (char *) &tailleMess);

        s.sendData(sizeof(qint16), (char *) &Global::POS);
        r.getData(sizeof(float)*tailleMess, (char *) &part_data->pos[0]);

        initialiseReception();

        s.sendData(sizeof(qint16), (char *) &Global::ENDMESSAGES);
        r.getData(sizeof(qint16), (char *) &tag);
        if(tag == Global::ENDMESSAGES) {

            initialiseReception();

            if (load_vel) {
                s.sendData(sizeof(qint16), (char *) &Global::SELECTVEL);

                s.sendData(sizeof(qint16), (char *) &Global::SIZEARRAY);
                r.getData(sizeof(int), (char *) &tailleMess);

                s.sendData(sizeof(qint16), (char *) &Global::POS);
                r.getData(sizeof(float)*tailleMess, (char *) &part_data->vel[0]);

                initialiseReception();

                s.sendData(sizeof(qint16), (char *) &Global::ENDMESSAGES);
                r.getData(sizeof(qint16), (char *) &tag);
                if(tag != Global::ENDMESSAGES) {
                    std::cerr<<"Erreur sur le tag : "<<tag<<"\n";
                    initialiseReception();
                    deco();
                }
            }

        }
        else {
            std::cerr<<"Erreur sur le tag : "<<tag<<"\n";
            initialiseReception();
            deco();
        }
        return 1;

    }

    // ============================================================================
    // Return any error of connection
    // ============================================================================
    void Serveur::connexionUnsuccefull(QAbstractSocket::SocketError error) {

        QString mess;

        switch(error) {
                        case QAbstractSocket::ConnectionRefusedError:
            mess = tr("The server is not responding. \n Reconnect you.");
            break;

                        case QAbstractSocket::RemoteHostClosedError:
            mess = tr("Communication with the server was interrupted.");
            break;

                        case QAbstractSocket::HostNotFoundError:
            mess = tr("The server is not found.");
            break;

                        case QAbstractSocket::SocketAccessError:
            mess = tr("Error in attempting to create the socket");
            break;

                        case QAbstractSocket::SocketResourceError:
            mess = tr("Too many current connection. \n Please free resources.");
            break;

                        case QAbstractSocket::SocketTimeoutError:
            mess = "Time out on the socket";
            break;

                        case QAbstractSocket::NetworkError:
            mess = tr("Network problem");
            break;

                        default:
            mess = "Unknown problem : " + QString::number(error);
            break;
        }

        if(mess != "Time out on the socket") {
            emit connexionFaild(mess);
        }

    }

    // ============================================================================
    // Return the state of connection
    // ============================================================================
    bool Serveur::getEtatConnexion() {
        return (socket->state() == QAbstractSocket::ConnectedState);
    }

    // ============================================================================
    // Disconnect the client
    // ============================================================================
    void Serveur::deco() {

        socket->disconnectFromHost();

    }


}
