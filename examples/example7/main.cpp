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

using namespace std;

typedef MRGraph_2D_4C<int, int, int> MRGrid;

int main() {
    MRGrid* mgrid;
    ifstream graphData("test1.mrmax");
    string input;
    while(graphData.good()) {
        graphData >> input;
        if (input == "c") {
            graphData.ignore(10000, '\n');
            continue;
        } else if (input == "p") {
            int width, height, layers;
            graphData >> width;
            graphData >> height;
            graphData >> layers;
            mgrid = new MRGrid(width, height, layers);
            graphData.ignore(10000, '\n');
            continue;
        } else if (input == "n") {
            int x, y, n, w, s, e;
            graphData >> x;
            graphData >> y;
            graphData >> n;
            graphData >> e;
            graphData >> s;
            graphData >> w;
            x=x-1;
            y=y-1;
            mgrid->set_neighbor_cap(mgrid->node_id(x,y,0), +1, 0, e);
            mgrid->set_neighbor_cap(mgrid->node_id(x,y,0), 0, +1, s);
            mgrid->set_neighbor_cap(mgrid->node_id(x,y,0), -1, 0, w);
            mgrid->set_neighbor_cap(mgrid->node_id(x,y,0), 0, -1, n);
            graphData.ignore(10000, '\n');
            continue;
        } else if (input == "t") {
            int x, y, s, t;
            graphData >> x;
            graphData >> y;
            graphData >> s;
            graphData >> t;
            x=x-1;
            y=y-1;
            mgrid->set_terminal_cap(mgrid->node_id(x,y,0), s, t);
            graphData.ignore(10000, '\n');
            continue;
        }
        else if (input == "e") {
            break;
        }
    }

    mgrid->compute_maxflow();
    for(int y=0; y<3; y++) {
        for(int x=0; x<3; x++) {
            cout << ((mgrid->get_segment(x,y) == 0) ? 'S' : 'T');
        }
        cout << endl;
    }
    cout << mgrid->get_flow() << endl;

    delete mgrid;
}
