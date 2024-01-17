/**
 * Simulationarchive Viewer
 *
 * This example allows you load in a Simulationarchive and visualize it.
 * You can use the keyboard to step through the individual snapshots.
 * It work with the web-based visualization as well as with OpenGL.
 *
 */
#include <nemo.h>
#include "rebound.h"

string defv[] = {
    "in=???\n                     input SimulationArchive",
    "server=1234\n                server port",
    "VERSION=0.1\n                17-jan-2024 PJT",
    NULL,
};

string usage="view a SimulationArchive in a web browser app";


struct reb_simulationarchive* sa;
void heartbeat(struct reb_simulation* const r);

int64_t current_snapshot = 0;

int key_callback(struct reb_simulation* r, int key){
    switch (key){
        case 262: // right arrow
            current_snapshot++;
            break;
        case 263: // left arrow
            current_snapshot--;
            break;
        case 268: // home
            current_snapshot = 0;
            break;
        case 269: // end
            current_snapshot = sa->nblobs - 1;
            break;
        case 266: // page up
            current_snapshot -= 10;
            break;
        case 267: // page down
            current_snapshot += 10;
            break;
        default: // unknown key
            return 0; // check default keys
    }

    // Update simulation
    if (current_snapshot < 0){
        current_snapshot = 0;
    }
    if (current_snapshot >= sa->nblobs){
        current_snapshot = sa->nblobs - 1;
    }
    r->status = REB_STATUS_SUCCESS; // will trigger reb_simulation_integrate to exit
    return 1;
}

void nemo_main()
{
    string infile = getparam("in");
    int server = getiparam("server");
    sa = reb_simulationarchive_create_from_file(infile);
    if (!sa) error("Error loading Simulationarchive from file `%s`",infile);

    printf("Simulationarchive loaded from file `%s`.\n",infile);
    printf("Number of snapshots: %ld.\n", sa->nblobs);
    printf("You can step through the Simulationarchive using the following keys in the visualization window:\n");
    printf(" Right arrow: jump to next snapshot\n");
    printf(" Left arrow:  jump to previous snapshot\n");
    printf(" Page down:   jump 10 snapshots foward\n");
    printf(" Page up:     jump 10 snapshots backward\n");
    printf(" Home key:    jump to first snapshot\n");
    printf(" End key:     jump to last snapshot\n\n");

    while(1){
        printf("Loading snapshot %ld.\n", current_snapshot);
        struct reb_simulation* r = reb_simulation_create_from_simulationarchive(sa, current_snapshot);
        if (!r)
            error("Error loading Simulation from Simulationarchive");

        r->key_callback = key_callback;
        r->status = REB_STATUS_PAUSED;

        // This allows you to connect to the simulation using
        // a web browser by pointing it to http://localhost:1234
        reb_simulation_start_server(r, server);

        // Not actually integrating because simulation is paused. 
        reb_simulation_integrate(r, INFINITY);

        if (r->status > 0){ // quit
            reb_simulation_free(r);
            break;
        }

        reb_simulation_free(r);
    }
    reb_simulationarchive_free(sa);
}
