//#########  Example 2

#include <iostream>
#include "..\..\include\MRGraph_2D_4C.h"

/*
    s
    |
   <2>
    |
    o-<2>-o-<2>-o-<2>-o-<2>-o-<2>-o-<2>-o-<2>-o-<1>-o-<2>-o
                                                          |
                                                         <2>
                                                          |
                                                          t
*/
int main()
{
    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> Grid;
    Grid grid(10,1,1);

    grid.set_neighbor_cap(grid.node_id(0,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(1,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(2,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(3,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(4,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(5,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(6,0,0), +1, 0, 2);
    grid.set_neighbor_cap(grid.node_id(7,0,0), +1, 0, 1);
    grid.set_neighbor_cap(grid.node_id(8,0,0), +1, 0, 2);

    grid.set_terminal_cap(grid.node_id(0,0,0),2,0);
    grid.set_terminal_cap(grid.node_id(9,0,0),0,2);

    grid.print();
    grid.compute_maxflow();
    std::cout << std::endl << "Segmentation:" << std::endl;
    for(int y=0; y<1; y++) {
        for(int x=0; x<10; x++) {
            std::cout << ((grid.get_segment(x,y) == 0) ? 'S' : 'T');
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "max flow: " << grid.get_flow() << std::endl;
    return 0;
}
