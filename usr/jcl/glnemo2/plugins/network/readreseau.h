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
#ifndef READRESEAU_H
#define READRESEAU_H

#include <QObject>
#include <QTcpSocket>

namespace network {

    class readReseau : public QObject {

        Q_OBJECT

        public:
            readReseau();
            readReseau(QTcpSocket *Client);
            void setSocket(QTcpSocket *Client);
            void getData(int size, char * tab);

        private:
            QTcpSocket *client;

    };

}

#endif // READRESEAU_H
