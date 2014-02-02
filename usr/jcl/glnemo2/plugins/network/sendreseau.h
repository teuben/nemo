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
#ifndef SENDRESEAU_H
#define SENDRESEAU_H

#include <QObject>
#include <QTcpSocket>

namespace network {

    class sendReseau : public QObject {

        Q_OBJECT

        public:
            sendReseau();
            void setSocket(QTcpSocket *Serveur);
            void sendData(int taille, char *val);

        private:
            QTcpSocket *client;

    };

}

#endif // SENDRESEAU_H
