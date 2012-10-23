#ifndef GRAPHWRITER_H
#define GRAPHWRITER_H

#include <fstream>
#include "MRGraph_2D_4C.h"

template <typename tcap_t, typename ncap_t, typename flow_t>
void write_grid(MRGraph_2D_4C<tcap_t, ncap_t, flow_t> &mgraph, const char* filename) {
    std::ofstream outfile (filename);
    for(int y=0; y<mgraph.get_height(); y++) {
        for(int x=0; x<mgraph.get_width(); x++) {
            outfile << int(mgraph.get_segment(x,y)) << ",";
        }
        outfile << std::endl;
    }
    outfile.close();
}

#endif // GRAPHWRITER_H
