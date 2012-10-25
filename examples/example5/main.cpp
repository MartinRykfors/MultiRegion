/*
 ##### Example 5 #####
 Run MRGraph on a 2-layered 2x2 grid

   s--8--o----1----o
         |\        |\
         | 8       | 8
         8  \      1  \
         |   o----1----o
         |   |     |   |
         |   1     |   |
         o---|-8---o   8
          \  |      \  |
           1 |       1 |
            \|        \|
             o----8----o--8--t


*/
#include <iostream>
#include <cstdlib>
#include "..\..\include\MRGraph_2D_4C.h"
#include "..\..\include\GraphParser.h"


int main()
{
    typedef MRGraph_2D_4C<uint16_t, uint16_t, uint16_t> Grid;
    typedef GraphParser<uint16_t, uint16_t, uint16_t> Parser;

    Parser parser("2x2x2grid.mrmax");
    Grid grid(parser.get_width(), parser.get_height(), parser.get_layers());
    parser.parse_data(grid);
/*
    Grid grid(2,2,2);
    // W -- E edges
    grid.set_neighbor_cap(grid.node_id(0,0,0),+1, 0, 1);
    grid.set_neighbor_cap(grid.node_id(1,0,0),-1, 0, 1);

    grid.set_neighbor_cap(grid.node_id(0,1,0),+1, 0, 1);
    grid.set_neighbor_cap(grid.node_id(1,1,0),-1, 0, 1);

    grid.set_neighbor_cap(grid.node_id(0,0,1),+1, 0, 8);
    grid.set_neighbor_cap(grid.node_id(1,0,1),-1, 0, 8);

    grid.set_neighbor_cap(grid.node_id(0,1,1),+1, 0, 8);
    grid.set_neighbor_cap(grid.node_id(1,1,1),-1, 0, 8);

    // N -- S edges
    grid.set_neighbor_cap(grid.node_id(0,0,0), 0, +1, 8);
    grid.set_neighbor_cap(grid.node_id(0,1,0), 0, -1, 8);

    grid.set_neighbor_cap(grid.node_id(1,0,0), 0, +1, 8);
    grid.set_neighbor_cap(grid.node_id(1,1,0), 0, -1, 8);

    grid.set_neighbor_cap(grid.node_id(0,0,1), 0, +1, 1);
    grid.set_neighbor_cap(grid.node_id(0,1,1), 0, -1, 1);

    grid.set_neighbor_cap(grid.node_id(1,0,1), 0, +1, 1);
    grid.set_neighbor_cap(grid.node_id(1,1,1), 0, -1, 1);

    // Down -- Up edges
    grid.set_interior_cap(grid.node_id(0,0,0), +1, 8);
    grid.set_interior_cap(grid.node_id(0,0,1), -1, 8);

    grid.set_interior_cap(grid.node_id(1,0,0), +1, 1);
    grid.set_interior_cap(grid.node_id(1,0,1), -1, 1);

    grid.set_interior_cap(grid.node_id(0,1,0), +1, 1);
    grid.set_interior_cap(grid.node_id(0,1,1), -1, 1);

    grid.set_interior_cap(grid.node_id(1,1,0), +1, 8);
    grid.set_interior_cap(grid.node_id(1,1,1), -1, 8);

    // Set terminal caps
    grid.set_terminal_cap(grid.node_id(0,0,0), 8, 0);
    grid.set_terminal_cap(grid.node_id(1,1,1), 0, 8);
*/
    grid.compute_maxflow();

    char* space = "    ";
    // print upper layer
    std::cout <<        (grid.get_segment(0,0) % 2 == 0 ?'S':'T') << space
              <<        (grid.get_segment(1,0) % 2 == 0 ?'S':'T') << std::endl;
    std::cout << " " << (grid.get_segment(0,1) % 2 == 0 ?'S':'T') << space
              <<        (grid.get_segment(1,1) % 2 == 0 ?'S':'T') << std::endl;
    std::cout << std::endl;
    // print lower layer
    std::cout <<        (grid.get_segment(0,0) / 2 == 0 ?'S':'T') << space
              <<        (grid.get_segment(1,0) / 2 == 0 ?'S':'T') << std::endl;
    std::cout << " " << (grid.get_segment(0,1) / 2 == 0 ?'S':'T') << space
              <<        (grid.get_segment(1,1) / 2 == 0 ?'S':'T') << std::endl;

    std::cout << "Max flow: " <<grid.get_flow() << std::endl;;
}
