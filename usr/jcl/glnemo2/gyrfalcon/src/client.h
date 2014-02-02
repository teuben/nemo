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
#ifndef CLIENT_H
#define CLIENT_H

#include <QThread>
#include <QMutex>
#include <QWaitCondition>
#include <body.h>

#include "readreseau.h"
#include "sendreseau.h"

class QTcpSocket;

class Client : public QThread {


    public:
        Client(int socketDesc, const falcON::snapshot *_my_snap, QMutex *_mutex, QWaitCondition *_condition, int *nbConnectfree);
        ~Client();
        void deconnectClient();

    private:
        int socketDescriptor;
        QTcpSocket *socket;
        QMutex *mutex;
        QWaitCondition *condition;
        const falcON::snapshot *S;
        int nbBody;
        int NbBodySelect;
        std::string select_part;
        QList<int> *interval;
        float *tablePart;
        network::readReseau r;
        network::sendReseau s;
        int tailleMess;
        qint16 tag;
        QString message;
        int *nbConnectFree;
        bool Freeplace;
        bool vitesse;

        void run();
        void initialiseMessage();
        void readTag();
        void bonjour();
        void getBody();
        void getTime();
        void setSelectBody();
        void getSizeArray();
        void getPositions();
        int setInterval(std::string *_interval);
        int addToInterval(int bodyStart, int bodyStop);

};

#endif // CLIENT_H
