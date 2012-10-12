/*
##### Example 4 #####

Demonstration of MGRAPH_VERBOSE on a small 2D grid.
Legend:
    s       Node connected directly to the source
    t       Node connected directly to the sink
    a       Active nodes of the growing tree
    f       Found nodes
    o       Orphaned nodes
    >v<^    Bridge node on the s-tree when a path is found.
    A       Active node being scanned
    O       Orphan being processed
    !       Tainted node (nodes marked active that shouldn't be)
*/

#include <iostream>
#include <cstdlib>
#define MRGRAPH_VERBOSE
#define MRGRAPH_VERBOSE_START 0

#include "..\..\include\MRGraph_2D_4C.h"

#define WIDTH 7
#define HEIGHT 7
#define MAXCAP 10


using namespace std;

int main()
{
    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> MRGrid;

    MRGrid mgrid(WIDTH,HEIGHT,1);
    //set random east caps
    for (int x = 0; x < WIDTH-1; ++x){
        for (int y = 0; y < HEIGHT; ++y){
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), +1, 0, val);
        }
    }

    //set random north caps
    for (int x = 0; x < WIDTH; ++x){
        for (int y = 1; y < HEIGHT; ++y){
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), 0, -1, val);
        }
    }

    //set random west caps
    for (int x = 1; x < WIDTH; ++x){
        for (int y = 0; y < HEIGHT; ++y){
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), -1, 0, val);
        }
    }

    //set random south caps
    for (int x = 0; x < WIDTH; ++x){
        for (int y = 0; y < HEIGHT-1; ++y){
            uint16_t val = rand()%MAXCAP;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0), 0, +1, val);
        }
    }


    for (int y = 0; y < HEIGHT; ++y) {
        mgrid.set_terminal_cap(mgrid.node_id(0,y,0), 40, 0);
        mgrid.set_terminal_cap(mgrid.node_id(WIDTH-1,y,0), 0, 40);
    }

//    //set random terminal caps
//    for (int x = 0; x < WIDTH; ++x) {
//        for (int y = 0; y < HEIGHT; ++y) {
//            uint16_t val1 = rand()%MAXCAP;
//            uint16_t val2 = rand()%MAXCAP;
//            mgrid.set_terminal_cap(mgrid.node_id(x,y,0), val1, val2);
//        }
//    }


    mgrid.compute_maxflow();

    std::cout << std::endl;
    for(int y=0; y<HEIGHT; y++) {
        for(int x=0; x<WIDTH; x++) {
            std::cout << char(mgrid.get_segment(x,y)+176); //((mgrid.get_segment(x,y) == 0) ? char(176) : char(178));
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    mgrid.print_grid("Final trees.");
}
