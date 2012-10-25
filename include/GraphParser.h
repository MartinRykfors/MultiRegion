#ifndef GRAPHPARSER_H
#define GRAPHPARSER_H

#include <string>
#include <fstream>
#include <iostream>
#include "MRGraph_2D_4C.h"

//todo, jätteviktigt. saknar stöd för i-lines. test1.mrmax är fel för den saknar LAYER-koordinater.
//matlab-generatorn behöver också fixa detta.

template <typename tcap_t, typename ncap_t, typename flow_t>
class GraphParser {
    public:
    GraphParser(const char* filename);
    void parse_data(MRGraph_2D_4C<tcap_t, ncap_t, flow_t> & mgraph);
    int get_width() {
        return width;
    };
    int get_height() {
        return height;
    };
    int get_layers() {
        return layers;
    };

private:
    int width, height, layers;
    std::ifstream graph_stream;

};

template <typename tcap_t, typename ncap_t, typename flow_t>
GraphParser<tcap_t, ncap_t, flow_t>::
GraphParser(const char* filename):
    graph_stream(filename) {
    std::string input;
    while(graph_stream.good()) {
        graph_stream >> input;
        if (input == "c") {
            graph_stream.ignore(10000, '\n');
            continue;
        } else if (input == "p") {
            graph_stream >> width;
            graph_stream >> height;
            graph_stream >> layers;
            graph_stream.ignore(10000, '\n');
            break;
        }
    }
}

template <typename tcap_t, typename ncap_t, typename flow_t>
void GraphParser<tcap_t, ncap_t, flow_t>::
parse_data(MRGraph_2D_4C<tcap_t, ncap_t, flow_t>& mgrid) {
    std::string input;
    while(graph_stream.good()) {
        graph_stream >> input;
        if (input == "c") {
            graph_stream.ignore(10000, '\n');
            continue;
        } else if (input == "n") {
            int x, y, d;
            ncap_t n, w, s, e;
            graph_stream >> x;
            graph_stream >> y;
            graph_stream >> d;
            graph_stream >> n;
            graph_stream >> e;
            graph_stream >> s;
            graph_stream >> w;
            x=x-1;
            y=y-1;
            d=d-1;
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,d), +1, 0, e);
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,d), 0, +1, s);
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,d), -1, 0, w);
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,d), 0, -1, n);
            graph_stream.ignore(10000, '\n');
            continue;
        } else if (input == "t") {
            int x, y, d;
            tcap_t s, t;
            graph_stream >> x;
            graph_stream >> y;
            graph_stream >> d;
            graph_stream >> s;
            graph_stream >> t;
            x=x-1;
            y=y-1;
            d=d-1;
            mgrid.set_terminal_cap(mgrid.node_id(x,y,d), s, t);
            graph_stream.ignore(10000, '\n');
            continue;
        } else if (input == "i") {
            int x, y, d;
            graph_stream >> x;
            graph_stream >> y;
            graph_stream >> d;
            x=x-1;
            y=y-1;
            d=d-1;
            for (int i = 0; i < layers; ++i) {
                int offset = i - d;
                if (offset == 0) {
                    continue;
                }
                ncap_t cap;
                graph_stream >> cap;
                mgrid.set_interior_cap(mgrid.node_id(x,y,d),offset,cap);
            }

            graph_stream.ignore(10000, '\n');
            continue;
        }
        else if (input == "e") {
            break;
        }
    }
}

#endif // GRAPHPARSER_H
