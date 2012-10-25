/*
 ##### Example 6 #####
 Run MRGraph on a 1x1 4-layered grid

            s
            |
            8
            |
            o--
            |\ \
            0 \ \
            |  | \
            o  2  \
           /|  |   |
          / 2 /    1
         |  |/     |
         1  o     /
         |  |    /
          \ 1   /
           \|  /
            o--
            |
            8
            |
            t
*/
#include <iostream>
#include <cstdlib>
#include <exception>
#include "..\..\include\MRGraph_2D_4C.h"
#include "..\..\include\GraphParser.h"

int main()
{


    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> Grid;
    typedef GraphParser<uint16_t, uint16_t, uint16_t> Parser;

    Parser parser("complete4.mrmax");
    Grid grid(parser.get_width(), parser.get_height(), parser.get_layers());
    parser.parse_data(grid);
    /*Grid grid(1,1,4);

    grid.set_interior_cap(grid.node_id(0,0,0), 1, 0);
    grid.set_interior_cap(grid.node_id(0,0,0), 2, 2);
    grid.set_interior_cap(grid.node_id(0,0,0), 3, 1);

    grid.set_interior_cap(grid.node_id(0,0,1),-1, 0);
    grid.set_interior_cap(grid.node_id(0,0,1), 1, 2);
    grid.set_interior_cap(grid.node_id(0,0,1), 2, 1);

    grid.set_interior_cap(grid.node_id(0,0,2),-2, 2);
    grid.set_interior_cap(grid.node_id(0,0,2),-1, 2);
    grid.set_interior_cap(grid.node_id(0,0,2), 1, 1);

    grid.set_interior_cap(grid.node_id(0,0,3),-3, 1);
    grid.set_interior_cap(grid.node_id(0,0,3),-2, 1);
    grid.set_interior_cap(grid.node_id(0,0,3),-1, 1);

    grid.set_terminal_cap(grid.node_id(0,0,0), 8, 0);
    grid.set_terminal_cap(grid.node_id(0,0,3), 0, 8);
*/
    grid.compute_maxflow();

    for (int i = 0; i < 4; ++i) {
        std::cout << (((grid.get_segment(0,0) >> i) & 1) == 0? 'S' : 'T') << std::endl;
    }
    std::cout << std::endl << "max flow: " << grid.get_flow() << std::endl;
    return 0;
}
