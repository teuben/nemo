// ============================================================================
// Copyright Jean-Charles LAMBERT - 2007-2012
//           Yannick Dalbin
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           Pôle de l'Etoile, site de Château-Gombert
//           38, rue Frédéric Joliot-Curie
//           13388 Marseille cedex 13 France
//           CNRS U.M.R 7326
// ============================================================================
// See the complete license in LICENSE and/or "http://www.cecill.info".
// ============================================================================
#ifndef MASTER_SERVER_THREAD_H
#define MASTER_SERVER_THREAD_H

#include <QThread>
#include <body.h>

class MasterServer;

class MasterServerThread : public QThread {

    Q_OBJECT

    public:
        MasterServerThread(std::string _sim_name, int _port, int _max_port, const falcON::snapshot * S);
        ~MasterServerThread();
        void updateData();

    private:
        MasterServer *server;
        int port; //Proprit qui permet de savoir sur quel port l'coute doti tre effectu
        int nbConnexionAutorise;
        const falcON::snapshot * my_snapshot;
        std::string sim_name;
        void run();

};

#endif // MASTER_SERVER_THREAD_H
