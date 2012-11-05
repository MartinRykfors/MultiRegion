#include "mex.h"
#include "MRGraph_2D_4C.h"

void print_dim(const mxArray *a) {
    mwSignedIndex m = (mwSignedIndex)mxGetM(a);
    mwSignedIndex n = (mwSignedIndex)mxGetN(a);
    mexPrintf("m x n = %d x %d\n",m,n);
}

void parse_single_layer(MRGraph_2D_4C<int,int,int> &mgrid, const mxArray *prhs[], const int width, const int height) {
    double *S, *T;
    S = mxGetPr(prhs[3]);
    T = mxGetPr(prhs[4]);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            mgrid.set_terminal_cap(mgrid.node_id(x,y,0),
                                   int(S[y+height*x]),
                                   int(T[y+height*x]));
        }
    }
    double* north = mxGetPr(prhs[5]);
    double test = north[1];
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0),0,-1,
                                   int(north[y+height*x]));
        }
    }

    double* east = mxGetPr(prhs[6]);
    for (int x = 0; x < width-1; ++x) {
        for (int y = 0; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0),1,0,
                                   int(east[y+height*x]));
        }
    }

    double* south = mxGetPr(prhs[7]);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height-1; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0),0, 1,
                                   int(south[y+height*x]));
        }
    }

    double* west = mxGetPr(prhs[8]);
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,0),-1,0,
                                   int(west[y+height*x]));
        }
    }
}

void parse_layer(MRGraph_2D_4C<int,int,int> &mgrid,
                 const mxArray *prhs[],
                 const int width,
                 const int height,
                 const int layer_idx) {
    double *S, *T;
    S = mxGetPr(prhs[3]);
    T = mxGetPr(prhs[4]);

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            mgrid.set_terminal_cap(mgrid.node_id(x,y,layer_idx),
                                   int(S[y + height*x + width*height*layer_idx]),
                                   int(T[y + height*x + width*height*layer_idx]));
        }
    }
    double* north = mxGetPr(prhs[5]);
    double test = north[1];
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,layer_idx),0,-1,
                                   int(north[y + height*x + width*height*layer_idx]));
        }
    }

    double* east = mxGetPr(prhs[6]);
    for (int x = 0; x < width-1; ++x) {
        for (int y = 0; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,layer_idx),1,0,
                                   int(east[y + height*x + width*height*layer_idx]));
        }
    }

    double* south = mxGetPr(prhs[7]);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height-1; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,layer_idx),0, 1,
                                   int(south[y + height*x + width*height*layer_idx]));
        }
    }

    double* west = mxGetPr(prhs[8]);
    for (int x = 0; x < width; ++x) {
        for (int y = 1; y < height; ++y) {
            mgrid.set_neighbor_cap(mgrid.node_id(x,y,layer_idx),-1,0,
                                   int(west[y + height*x + width*height*layer_idx]));
        }
    }
}

void parse_interior(MRGraph_2D_4C<int,int,int> &mgrid,
                    const mxArray *prhs[],
                    const int width,
                    const int height,
                    const int layers) {
    double* weights = new double[layers-1];
    double* I = mxGetPr(prhs[9]);
    for (int i = 0; i < layers; ++i) {
        for (int j = 0; j < layers-1; ++j) {
            weights[j] = I[j + i*(layers-1)];
            mexPrintf("wieghts[%d] = %d\n",j,int(weights[j]));
        }
        for (int x = 0; x < width; ++x) {
            for (int y = 0; y < height; ++y) {
                for (int j = 0; j < layers; ++j) {
                    if (j == i) {
                        continue;
                    }
                    int idx = (j < i) ? j : j-1;
                    //mexPrintf("Setting interior at (%d, %d, %d): %d\n",x,y,i,int(weights[idx]));
                    mgrid.set_interior_cap(mgrid.node_id(x,y,i),j-i,int(weights[idx]));
                }
            }
        }
    }
    delete[] weights;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Check for proper number of arguments */
    /*if (nrhs != 10) {
        mexErrMsgTxt("Wrong number of input arguments. (Ten are needed)\n");
    } else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }*/

    const int width  = int(*mxGetPr(prhs[0]));
    const int height = int(*mxGetPr(prhs[1]));
    const int layers = int(*mxGetPr(prhs[2]));
    MRGraph_2D_4C<int,int,int> mgrid(width,height,layers);
    for (int i = 0; i < layers; ++i) {
        parse_layer(mgrid,prhs,width,height,i);
    }
    if (layers > 1) {
        parse_interior(mgrid, prhs, width, height, layers);
    }
    mgrid.compute_maxflow();
    mexPrintf("Maxflow: %d\n",mgrid.get_flow());
    plhs[0] = mxCreateNumericMatrix(height,width, mxDOUBLE_CLASS, mxREAL);
    double* output = mxGetPr(plhs[0]);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            output[y+height*x]=mgrid.get_segment(x,y);
        }
    }
    return;
}
