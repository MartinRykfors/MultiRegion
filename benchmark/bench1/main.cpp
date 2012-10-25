/*
#### Benchmark 1 ####
Use std::clock to time how long calculate_maxflow takes.
*/
#include <iostream>
#include <ctime>
#include <cstdlib>
#include "..\..\include\MRGraph_2D_4C.h"
#include "..\..\include\GridGraph_2D_4C.h"

#define NUM_TRIALS 10
#define WIDTH 512
#define HEIGHT 512
#define MAXCAP 600

using namespace std;

typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> MRGrid;
typedef GridGraph_2D_4C<uint16_t, uint16_t, uint16_t> GGrid;

void create_grid(MRGrid& mgrid, GGrid& ggrid) {
    //set random east caps
    for (int x = 0; x < WIDTH-1; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), +1, 0, val);
            ggrid.set_neighbor_cap(ggrid.node_id(x,y),   +1, 0, val);
        }
    }

    //set random north caps
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 1; y < HEIGHT; ++y) {
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), 0, -1, val);
            ggrid.set_neighbor_cap(ggrid.node_id(x,y),   0, -1, val);
        }
    }

    //set random west caps
    for (int x = 1; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), -1, 0, val);
            ggrid.set_neighbor_cap(ggrid.node_id(x,y),   -1, 0, val);
        }
    }

    //set random south caps
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT-1; ++y) {
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), 0, +1, val);
            ggrid.set_neighbor_cap(ggrid.node_id(x,y),   0, +1, val);
        }
    }

    //set random terminal caps
    for (int x = 0; x < WIDTH; ++x) {
        for (int y = 0; y < HEIGHT; ++y) {
            uint16_t val1 = rand()%MAXCAP;
            uint16_t val2 = rand()%MAXCAP;
            mgrid.set_terminal_cap(mgrid.node_id(x,y,0), val1, val2);
            ggrid.set_terminal_cap(ggrid.node_id(x,y),   val1, val2);
        }
    }
}

int main() {

    double g_time[NUM_TRIALS];
    double m_time[NUM_TRIALS];

    for (int i = 0; i < NUM_TRIALS; ++i) {
        MRGrid mgrid(WIDTH,HEIGHT,1);
        GGrid ggrid(WIDTH,HEIGHT);

        create_grid(mgrid, ggrid);
        clock_t start, end;

        start = clock();
        ggrid.compute_maxflow();
        end = clock();
        g_time[i] = (double(end) - double(start))/CLOCKS_PER_SEC;
        cout << "Trial " << (i+1) << " of " << NUM_TRIALS <<
            "\t Gridcut: " << g_time[i] << " s" <<
            "\t flow: " << ggrid.get_flow() << endl;

        start = clock();
        mgrid.compute_maxflow();
        end = clock();
        m_time[i] = (double(end) - double(start))/CLOCKS_PER_SEC;
        cout << "Trial " << (i+1) << " of " << NUM_TRIALS <<
            "\t MRGraph: " << m_time[i] << " s" <<
            "\t flow: " << mgrid.get_flow() << endl;
        cout << endl;
    }
    double m_avg = 0., g_avg = 0.;
    for (int i = 0; i < NUM_TRIALS; ++i) {
        m_avg += m_time[i];
        g_avg += g_time[i];
    }
    m_avg /= NUM_TRIALS;
    g_avg /= NUM_TRIALS;

    cout << "Average times on " << HEIGHT << " x " << WIDTH << " grid:" << endl <<
        "Gridcut: " << g_avg << endl <<
        "MRGraph: " << m_avg << endl <<
        "Speedup: " << (g_avg/m_avg) << endl;
}






