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
#ifndef SERVEUR_H
#define SERVEUR_H

#include <QThread>
#include <QMutex>
#include<QTcpSocket>
#include "componentrange.h"
#include "particlesdata.h"
#include "readreseau.h"
#include "sendreseau.h"

class QTcpSocket;

namespace network {

    class Serveur  : public QObject {

        Q_OBJECT

        public:
            Serveur();
            ~Serveur();
            void setConnexion(std::string, int);
            bool getEtatConnexion();
            void deco();
            float getTime();
            int getNbBody();
            bool isBeenConnected(int);
            int readPart(glnemo::ParticlesData * pdata, const int *index, const int nsel, const bool load_vel, const std::string _selectPart);

        private:
            QTcpSocket *socket;
            int tailleMess;
            qint16 tag;
            QString message;
            void initialiseReception();

            float time;
            int nbbody;
            readReseau r;
            sendReseau s;

        signals:
            void serveurFull();
            void connexionFaild(QString);
            void connectSuccess();

        private slots:
            void connexionSuccefull();
           void connexionUnsuccefull(QAbstractSocket::SocketError);

    };

}

#endif // SERVEUR_H
