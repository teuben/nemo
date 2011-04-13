// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2011
//           Yannick Dalbin
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           Pôle de l'Etoile, site de Château-Gombert
//           38, rue Frédéric Joliot-Curie
//           13388 Marseille cedex 13 France
//           CNRS U.M.R 6110
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".
// ============================================================================
#include<QTcpSocket>
#include<QByteArray>
#include<iostream>
#include <assert.h>
#include "readreseau.h"

namespace network {

    readReseau::readReseau() {


    }

    // ============================================================================
    // Set socket to read
    // ============================================================================
    void readReseau::setSocket(QTcpSocket* sock) {

        client = sock;

    }

    // ============================================================================
    // Read all of type of data from the socket
    // ============================================================================
    void readReseau::getData(int size, char * tab) {

        int read;
        int totalread=0;
        while (totalread<size && client->state() == QAbstractSocket::ConnectedState) {
            read=client->bytesAvailable();
            if (read>0) {
                if ((totalread+read)>size) {
                    client->read((char *) tab+totalread, size-totalread);
                    totalread += (size-totalread);
                    assert(totalread==size);
                } else {
                    client->read((char *) tab+totalread,read);
                    totalread += read;
                }
            } else {
                client->waitForReadyRead(300);
            }
        }

    }




}
