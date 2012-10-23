/*
##### Example 7 #####

Importing a grid graph from a .mrmax file;

*/
#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include "../../include/MRGraph_2D_4C.h"
#include "../../include/GraphParser.h"
#include "../../include/GraphWriter.h"

using namespace std;

typedef MRGraph_2D_4C<int, int, int> MRGrid;
typedef GraphParser<int, int, int> Parser;

int main() {
    Parser parser("test1.mrmax");
    MRGrid mgrid(parser.get_width(), parser.get_height(), parser.get_layers());
    parser.parse_data(mgrid);
    mgrid.compute_maxflow();
    write_grid(mgrid, "output.csv");
}
