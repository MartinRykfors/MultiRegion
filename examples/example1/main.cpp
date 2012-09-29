
#include <iostream>
#include <stdint.h>
#include "..\..\include\MRGraph_2D_4C.h"
/*
                    (S) <- source
                     |
                    <8>
                     |
                     o--<3>--o--<1>--o
                     |       |       |
                    <5>     <1>     <4>
                     |       |       |
                     o--<1>--o--<5>--o
                     |       |       |
                    <2>     <2>     <3>
                     |       |       |
                     o--<1>--o--<2>--o
                                     |
                                    <9>
                                     |
                            sink -> (T)

*/
int main() {
    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> Grid;

    Grid* grid = new Grid(3,3,1);

    grid->set_neighbor_cap(grid->node_id(0,0,0),+1, 0,3);
    grid->set_neighbor_cap(grid->node_id(1,0,0),-1, 0,3);
    grid->set_neighbor_cap(grid->node_id(1,0,0),+1, 0,1);
    grid->set_neighbor_cap(grid->node_id(2,0,0),-1, 0,1);

    grid->set_neighbor_cap(grid->node_id(0,0,0), 0,+1,5);
    grid->set_neighbor_cap(grid->node_id(0,1,0), 0,-1,5);
    grid->set_neighbor_cap(grid->node_id(1,0,0), 0,+1,1);
    grid->set_neighbor_cap(grid->node_id(1,1,0), 0,-1,1);
    grid->set_neighbor_cap(grid->node_id(2,0,0), 0,+1,4);
    grid->set_neighbor_cap(grid->node_id(2,1,0), 0,-1,4);

    grid->set_neighbor_cap(grid->node_id(0,1,0),+1, 0,1);
    grid->set_neighbor_cap(grid->node_id(1,1,0),-1, 0,1);
    grid->set_neighbor_cap(grid->node_id(1,1,0),+1, 0,5);
    grid->set_neighbor_cap(grid->node_id(2,1,0),-1, 0,5);

    grid->set_neighbor_cap(grid->node_id(0,1,0), 0,+1,2);
    grid->set_neighbor_cap(grid->node_id(0,2,0), 0,-1,2);
    grid->set_neighbor_cap(grid->node_id(1,1,0), 0,+1,2);
    grid->set_neighbor_cap(grid->node_id(1,2,0), 0,-1,2);
    grid->set_neighbor_cap(grid->node_id(2,1,0), 0,+1,3);
    grid->set_neighbor_cap(grid->node_id(2,2,0), 0,-1,3);

    grid->set_neighbor_cap(grid->node_id(0,2,0),+1, 0,1);
    grid->set_neighbor_cap(grid->node_id(1,2,0),-1, 0,1);
    grid->set_neighbor_cap(grid->node_id(1,2,0),+1, 0,2);
    grid->set_neighbor_cap(grid->node_id(2,2,0),-1, 0,2);

    grid->set_terminal_cap(grid->node_id(0,0,0),8,0);
    grid->set_terminal_cap(grid->node_id(2,2,0),0,9);

    grid->compute_maxflow();

    for(int y=0; y<3; y++) {
        for(int x=0; x<3; x++) {
            std::cout << ((grid->get_segment(x,y) == 0) ? 'S' : 'T');
        }
        std::cout << std::endl;
    }
    std::cout << grid->get_flow() << std::endl;

    delete grid;

}
