/*
##### Example 3 #####

Evaluate GridCut and MRGraph on a random 2D-grid
The GridCut source code can be obtained at gridcut.com

*/

#include <iostream>
#include <cstdlib>
#include "..\..\include\MRGraph_2D_4C.h"
#include "..\..\include\GridGraph_2D_4C.h"


#define WIDTH 18
#define HEIGHT 18
#define MAXCAP 9


using namespace std;

int main() {
    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> MRGrid;
    typedef GridGraph_2D_4C<uint16_t, uint16_t, uint16_t> GGrid;

    MRGrid mgrid(WIDTH,HEIGHT,1);
    GGrid ggrid(WIDTH,HEIGHT);

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
    for (int y = 0; y < HEIGHT; ++y) {
        mgrid.set_terminal_cap(mgrid.node_id(0,y,0), 400, 0);
        ggrid.set_terminal_cap(ggrid.node_id(0,y),   400, 0);

        mgrid.set_terminal_cap(mgrid.node_id(WIDTH-1,y,0), 0, 400);
        ggrid.set_terminal_cap(ggrid.node_id(WIDTH-1,y),   0, 400);
    }



    ggrid.compute_maxflow();

    for(int y=0; y<HEIGHT; y++) {
        for(int x=0; x<WIDTH; x++) {
            std::cout << ((ggrid.get_segment(ggrid.node_id(x,y)) == 0) ? char(177) : char(178));
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    mgrid.compute_maxflow();
    for(int y=0; y<HEIGHT; y++) {
        for(int x=0; x<WIDTH; x++) {
            std::cout << char(mgrid.get_segment(x,y)+177); //((mgrid.get_segment(x,y) == 0) ? char(176) : char(178));
        }
        std::cout << std::endl;
    }
    std::cout << "Gridcut: " << ggrid.get_flow() << std::endl;
    std::cout << "MRGrid : " << mgrid.get_flow() << std::endl;
}
