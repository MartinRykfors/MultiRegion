#ifndef MRGRAPH_2D_4C_H
#define MRGRAPH_2D_4C_H

#ifndef MRGRAPH_VERBOSE_START
#define MRGRAPH_VERBOSE_START 0
#endif

#include <stdint.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <stdexcept>


template
<
typename tcap_t, // Type used to represent capacities of edges between nodes and terminals.
typename ncap_t, // Type used to represent capacities of edges between nodes and their neighbors.
typename flow_t  // Type used to represent the total flow.
>
class MRGraph_2D_4C {
public:
    //for indexing nodes
    typedef uint32_t index_t;

    //for enumerating outgoing edges
    typedef uint8_t dir_t;

    ~MRGraph_2D_4C();

    MRGraph_2D_4C(int width, int height, int depth);

    // Returns the number of layers of the graph
    int get_depth() const {
        return D;
    };

    int get_width() const{
        return ow;
    }

    int get_height() const{
        return oh;
    }

    // Returns the index of the node at point (x,y), label z
    inline index_t node_id(int x,int y,int z) const {
        if (x<0 || x>=ow || y<0 || y>=oh || z < 0 || z>=D) {
            //throw
        }
        return nodeId(x+1,y+1,z);
    };

    // Sets the capacities of the source-->node and node-->sink edges.
    // This function can be called only once per node.
    // Capacities less than 0 are treated as 0
    inline void set_terminal_cap(index_t node_id, tcap_t cap_source, tcap_t cap_sink);

    // Sets the capacity of the edge between node and its neighbor at [offset_x,offset_y].
    // For example, to set the capacity of the edge from node [x,y] to node [x+1,y], call:
    // graph->set_neighbor_cap(graph->node_id(x,y),+1,0,edge_capacity);
    inline void set_neighbor_cap(index_t node_id, int offset_x, int offset_y, ncap_t cap);

    // Sets the capacity of the edge between node and one other node in the pixel clique
    inline void set_interior_cap(index_t node_id, int offset_z, ncap_t cap);

    // Alternative way to set the edge capacities is to use the set_caps function which
    // sets capacities of all edges at once based on values from input arrays.
    // Each array has width*height elements, where each element corresponds to one node.
    // For example, cap_le[x+y*width] is the capacity of the outgoing edge from node [x,y]
    // to node [x-1,y], and cap_sink[x+y*width] is capacity of edge from node [x,y] to sink.
    template<typename type_arg_tcap,typename type_arg_ncap>
    inline void set_caps(const type_arg_tcap* cap_source,
                         const type_arg_tcap* cap_sink,
                         const type_arg_ncap* cap_le,         // [-1, 0]
                         const type_arg_ncap* cap_ge,         // [+1, 0]
                         const type_arg_ncap* cap_el,         // [ 0,-1]
                         const type_arg_ncap* cap_eg);        // [ 0,+1]


    // After the maxflow is computed, this function returns a bitfield
    // corresponding to to which region the pixel (x,y) belongs.
    // This differs form the GridCut interface.
    // Maybe 8 regions is not enough, but then just change the type.
    // Currently implemented as return label(x,y) = {0,1,2}
    inline uint8_t get_segment(int x, int y) const;

    // After the maxflow is computed, this function returns the value of maximum flow.
    inline flow_t get_flow() const;

    // Number crunch
    void compute_maxflow();

    // Print info about the graph
    void print();

    int get_x(index_t v) {
        return ((v/(64*D))%(W/8))*8 + (v%8);
    };
    int get_y(index_t v) {
        return v/(W*D*8)*8 + ((v/8)%8);
    };
    int get_z(index_t v) {
        return (v/64)%D;
    };
    //print a representation of the grid
    void print_grid(std::string msg,        //message to be printed
                    std::string options,    //o = orphans, a = active_t, A = active_s, f = found, q = query node
                    index_t sp_id,          //specific node to be hightlighted
                    char sp_icon);          //icon for specific node
private:

    struct Stack {
        index_t size;
        index_t* buffer;

        void push(index_t v) {
            buffer[size] = v;
            ++size;
        }

        index_t pop() {
            --size;
            return buffer[size];
        }

        bool is_empty() {
            return size==0;
        }

        void clear() {
            size = 0;
        }
    };

    const index_t ow; //image width
    const index_t oh; //image height
    const index_t W;  //grid width
    const index_t H;  //grid height
    const index_t WH; //nodes per z-layer
    const index_t D;  //number of z-layers (depth)
    const index_t N;  //total number of nodes

    const uint8_t DEGREE;   //number of outgoing edges for every node (not counting t-links
    const index_t XOFS;     //x-offset when going between blocks along x
    const index_t YOFS;

    typedef uint8_t label_t;

    static const label_t LABEL_F = 0;
    static const label_t LABEL_S = 1;
    static const label_t LABEL_T = 2;

    const dir_t DIR_NONE;
    const dir_t DIR_TERMINAL;

    const index_t PARENT_ID_TERMINAL;
    const index_t PARENT_ID_NONE;

    dir_t* REVERSE;         //would like this to be a const member
    index_t* neighbors;

    Stack orphans_s;
    Stack orphans_t;
    Stack active_s;
    Stack active_t;
    Stack found_s;
    Stack found_t;

    label_t* label;
    index_t* parent_id;     //nodindex för förälder
    dir_t*   parent_edge;   //kantindex (Parent -> Child för S, Child -> Parent för T
    index_t* dist;
    index_t s_dist;
    index_t t_dist;

    tcap_t* rc_st;          //residualcap för t-länkar
    ncap_t** rc_nbhd;       //residualcap c-f+f_rev
    flow_t maxflow;

    #ifdef MRGRAPH_VERBOSE
    int step;
    #endif

    //these return the index of nodes neighboring v
    inline index_t north(index_t v) const {
        return  (v & 56) == 0 ? v - YOFS : v-8;
    };
    inline index_t east (index_t v) const {
        return (~v & 7)  == 0 ? v + XOFS : v+1;
    };
    inline index_t south(index_t v) const {
        return (~v & 56) == 0 ? v + YOFS : v + 8;
    };
    inline index_t west (index_t v) const {
        return  (v & 7)  == 0 ? v - XOFS : v-1;
    };
    inline index_t down(index_t v, uint8_t n) const {
        return v + ((((v>>6)%D+n)%D)-(v>>6)%D)*64;
    };

    inline index_t nodeId(unsigned int x, unsigned int y, unsigned int z) const {
        return ((W/8)*(y/8)*D+x/8*D+z)*64 + 8*(y%8) + (x%8);
    }

    bool grow(index_t& s, index_t& t, dir_t& st);
    bool grow_s(index_t& s, index_t& t, dir_t& st);
    bool grow_t(index_t& s, index_t& t, dir_t& st);
    void augment(index_t s, index_t t, dir_t st);
    ncap_t min_residual_s(index_t s, ncap_t rcmin);
    ncap_t min_residual_t(index_t t, ncap_t rcmin);
    void augment_s(index_t s, const ncap_t rcmin);
    void augment_t(index_t t, const ncap_t rcmin);
    void rank_relabel_s(const index_t os, const index_t max_dist);
    void rank_relabel_t(const index_t ot, const index_t max_dist);
    dir_t get_origin(index_t v);
    void get_neighbors(index_t v);
    char get_icon(uint8_t conn);

    // hide assigment op, copy ctor, default ctor
    const MRGraph_2D_4C& operator=(const MRGraph_2D_4C &rhs) {};
    MRGraph_2D_4C(MRGraph_2D_4C &other) {};
    MRGraph_2D_4C() {};
};



template <typename tcap_t,typename ncap_t,typename flow_t>
MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
MRGraph_2D_4C(const int w, const int h, const int depth):
    ow(w),
    oh(h),
    W(((w+1)/8+1)*8),
    H(((h+1)/8+1)*8),
    WH(W*H),
    D(depth),
    N(W*D*H),
    DEGREE(4+depth-1),
    XOFS(64*D-7),
    YOFS(W/8*64*D-56),
    DIR_NONE(DEGREE),
    DIR_TERMINAL(DEGREE+1),
    PARENT_ID_TERMINAL(-1),
    PARENT_ID_NONE(-2),
    maxflow(0) {
    label = new label_t[N];
    parent_id = new index_t[N];
    parent_edge = new dir_t[N];
    rc_nbhd = new ncap_t*[N];   // en pekare per nod. i gridcut har de indexen
    // tvärt om -> en pekare per riktning.
    dist = new index_t[N];
    for(unsigned int i = 0; i < N; ++i) {
        rc_nbhd[i] = new ncap_t[DEGREE];
    }
    rc_st = new tcap_t[N];

    for(unsigned int i = 0; i < N; ++i) {
        label[i] = LABEL_F;
        parent_id[i] = PARENT_ID_NONE;
        parent_edge[i] = DIR_NONE;
        rc_st[i] = 0;
        dist[i] = 0;
        for (int j = 0; j < DEGREE; ++j) {
            rc_nbhd[i][j] = 0;
        }
    }

    REVERSE = new dir_t[DEGREE];
    REVERSE[0] = 2;
    REVERSE[1] = 3;
    REVERSE[2] = 0;
    REVERSE[3] = 1;
    for (int i = 4; i < DEGREE; ++i) {
        REVERSE[i] = DEGREE-1-(i-4);
    }
    orphans_s.buffer = new index_t[N];
    orphans_t.buffer = new index_t[N];
    active_s.buffer = new index_t[N];
    active_t.buffer = new index_t[N];
    found_s.buffer = new index_t[N];
    found_t.buffer = new index_t[N];

    neighbors = new index_t[DEGREE];

    orphans_s.clear();
    orphans_t.clear();
    active_s.clear();
    active_t.clear();
    found_s.clear();
    found_t.clear();

    #ifdef MRGRAPH_VERBOSE
    step = 0;
    #endif
}

template <typename tcap_t,typename ncap_t,typename flow_t>
MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
~MRGraph_2D_4C() {
    delete[] label;
    delete[] parent_id;
    delete[] parent_edge;
    delete[] rc_st;
    for(unsigned int i = 0; i < N; ++i) {
        delete[] rc_nbhd[i];
    }
    delete[] rc_nbhd;
    delete[] orphans_s.buffer;
    delete[] orphans_t.buffer;
    delete[] active_s.buffer;
    delete[] active_t.buffer;
    delete[] found_s.buffer;
    delete[] found_t.buffer;

    delete[] neighbors;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
flow_t MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
get_flow() const {
    return maxflow;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
inline uint8_t MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
get_segment(const int x, const int y) const {
    uint8_t result = 0;
    for (unsigned int i = 0; i < D; ++i){
        if(label[node_id(x,y,i)] == LABEL_T) {
            result += (1<<i);
        }
    }
    return result;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
set_terminal_cap(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t v,
                                                              const tcap_t  cap_s,
                                                              const tcap_t  cap_t) {
    if (cap_s > 0 && cap_t > 0) {
        if (cap_s > cap_t) {
            rc_st[v] = cap_s - cap_t;
            label[v] = LABEL_S;
            parent_edge[v] = DIR_TERMINAL;
            parent_id[v] = PARENT_ID_TERMINAL;
            maxflow += cap_t;
            dist[v] = 1;
            active_s.push(v);
            s_dist = 1;
        } else if (cap_s < cap_t) {
            rc_st[v] = cap_t - cap_s;
            label[v] = LABEL_T;
            parent_edge[v] = DIR_TERMINAL;
            parent_id[v] = PARENT_ID_TERMINAL;
            maxflow += cap_s;
            dist[v] = 1;
            active_t.push(v);
            t_dist = 1;
        } else {
            rc_st[v] = 0;
            label[v] = LABEL_F;
            parent_edge[v] = DIR_NONE;
            parent_id[v] = PARENT_ID_NONE;
            maxflow += cap_s;
        }
    } else {
        if (cap_s > 0) {
            rc_st[v] = cap_s;
            label[v] = LABEL_S;
            parent_edge[v] = DIR_TERMINAL;
            parent_id[v] = PARENT_ID_TERMINAL;
            dist[v] = 1;
            active_s.push(v);
            s_dist = 1;
        } else {
            rc_st[v] = cap_t;
            label[v] = LABEL_T;
            parent_edge[v] = DIR_TERMINAL;
            parent_id[v] = PARENT_ID_TERMINAL;
            dist[v] = 1;
            active_t.push(v);
            t_dist = 1;
        }
    }

}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
set_neighbor_cap(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t v,
                                                              const int     offset_x,
                                                              const int     offset_y,
                                                              const ncap_t  cap) {
    dir_t dir;
    if      (offset_x ==  1 && offset_y ==  0)
        dir = 1;
    else if (offset_x == -1 && offset_y ==  0)
        dir = 3;
    else if (offset_x ==  0 && offset_y ==  1)
        dir = 2;
    else if (offset_x ==  0 && offset_y == -1)
        dir = 0;
    else {
        throw std::range_error("Attempted to access unsupported direction");
    }
    rc_nbhd[v][dir] = cap;

}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
set_interior_cap(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t v,
                                                     const int     offset_z,
                                                     const ncap_t  cap) {
    uint8_t node_z = get_z(v);
    int target_z = node_z + offset_z;
    if (target_z < 0 || target_z >= D || offset_z == 0) {
        throw std::range_error("Bad z-offset for interior edge");
    }
    rc_nbhd[v][((offset_z + D) % D) + 3] = cap;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
compute_maxflow() {
    index_t vs, vt;
    dir_t st;
    bool path_found;
    while(1) {
        //Growth loop for s

#ifdef MRGRAPH_VERBOSE
        print_grid("About to start s-growth","A");
#endif

        while (1) {
            path_found = grow_s(vs,vt,st);
            if (!path_found) {
                //we've run out of active s-nodes, so the growth stage is now done
                break;
            }

#ifdef MRGRAPH_VERBOSE
            char arrow = ' ';
            switch (st) {
                case 0: arrow = char(94); break;
                case 1: arrow = char(62); break;
                case 2: arrow = char(118); break;
                case 3: arrow = char(60); break;

            }
            print_grid("Found a path", "", vs,arrow);
#endif

            augment(vs,vt,st);

#ifdef MRGRAPH_VERBOSE
            print_grid("Augmentation done, orphans (o) created.","o");
#endif

            while(!orphans_s.is_empty()) {
                index_t os = orphans_s.pop();
                rank_relabel_s(os, s_dist+1);
            }
            while(!orphans_t.is_empty()) {
                index_t ot = orphans_t.pop();
                rank_relabel_t(ot, t_dist);
            }
        }

        //make all found nodes active for the next s-growth
        std::swap(active_s.buffer,found_s.buffer);
        //swap size info
        std::swap(active_s.size,found_s.size);
        found_s.clear();
        ++s_dist;

        if (active_t.is_empty() && active_s.is_empty()) {
            //there are no nodes that can be made active for the next s-pass, so we are done
            break;
        }

#ifdef MRGRAPH_VERBOSE
        print_grid("About to start t-growth","a");
#endif
        //Growth loop for t
        while (1) {
            path_found = grow_t(vs,vt,st);
            if (!path_found) {
                break;
            }
#ifdef MRGRAPH_VERBOSE
            char arrow = ' ';
            switch (st) {
                case 0: arrow = char(94); break;
                case 1: arrow = char(62); break;
                case 2: arrow = char(118); break;
                case 3: arrow = char(60); break;

            }
            print_grid("Found a path", "", vs,arrow);
#endif

            augment(vs,vt,st);

#ifdef MRGRAPH_VERBOSE
            print_grid("Augmentation done, orphans (o) created.","o");
#endif
            while(!orphans_t.is_empty()) {
                index_t ot = orphans_t.pop();
                rank_relabel_t(ot, t_dist+1);
            }

            while(!orphans_s.is_empty()) {
                index_t os = orphans_s.pop();
                rank_relabel_s(os, s_dist);
            }
        }

        std::swap(active_t.buffer,found_t.buffer);
        std::swap(active_t.size,found_t.size);
        found_t.clear();
        ++t_dist;
        if (active_t.is_empty() && active_s.is_empty()) {
            break;
        }
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
augment(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t s,
        typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t t,
                                                      dir_t   st) {

    ncap_t rcmin = rc_nbhd[s][st];
    rcmin = min_residual_s(s,rcmin);
    rcmin = min_residual_t(t,rcmin);

    //augment with the minimum residual cap (creates orphans)
    augment_s(s,rcmin);
    augment_t(t,rcmin);
    rc_nbhd[s][st] -= rcmin;
    rc_nbhd[t][REVERSE[st]] += rcmin;
    maxflow += rcmin;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
ncap_t MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
min_residual_s(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t s,
                                                             ncap_t  rcmin) {
    while (parent_edge[s] != DIR_TERMINAL) {
        rcmin = (rcmin < rc_nbhd[parent_id[s]][parent_edge[s]] ?
                 rcmin : rc_nbhd[parent_id[s]][parent_edge[s]]);
        s = parent_id[s];
    }
    return rcmin < rc_st[s] ? rcmin : rc_st[s];
}

template <typename tcap_t,typename ncap_t,typename flow_t>
ncap_t MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
min_residual_t(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t t,
                                                             ncap_t  rcmin) {
    while (parent_edge[t] != DIR_TERMINAL) {
        rcmin = (rcmin < rc_nbhd[t][parent_edge[t]] ?
                 rcmin : rc_nbhd[t][parent_edge[t]]);
        t = parent_id[t];
    }
    return rcmin < rc_st[t] ? rcmin : rc_st[t];
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
augment_s(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t s,
                                                   const ncap_t rcmin) {

    //traverse the tree to source and alter the residual caps along the way
    while (parent_edge[s] != DIR_TERMINAL) {
        rc_nbhd[parent_id[s]][parent_edge[s]] -= rcmin;
        rc_nbhd[s][REVERSE[parent_edge[s]]] += rcmin;
        index_t ps = parent_id[s];

        //check for node being orphaned
        if (rc_nbhd[parent_id[s]][parent_edge[s]] == 0) {
            orphans_s.push(s);
            parent_id[s] = PARENT_ID_NONE;
            parent_edge[s] = DIR_NONE;
        }
        s = ps;
    }
    //augment s-link
    rc_st[s] -= rcmin;
    if(rc_st[s] == 0) {
        orphans_s.push(s);
        parent_id[s] = PARENT_ID_NONE;
        parent_edge[s] = DIR_NONE;
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
augment_t(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t t,
                                                   const ncap_t rcmin) {
    while (parent_edge[t] != DIR_TERMINAL) {
        rc_nbhd[t][parent_edge[t]] -= rcmin;
        rc_nbhd[parent_id[t]][REVERSE[parent_edge[t]]] += rcmin;
        index_t pt = parent_id[t];
        if (rc_nbhd[t][parent_edge[t]] == 0) {
            orphans_t.push(t);
            parent_id[t] = PARENT_ID_NONE;
            parent_edge[t] = DIR_NONE;
        }
        t = pt;
    }
    rc_st[t] -= rcmin;

    if(rc_st[t] == 0) {
        orphans_t.push(t);
        parent_id[t] = PARENT_ID_NONE;
        parent_edge[t] = DIR_NONE;
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
rank_relabel_s(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t os,
               const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t max_dist) {
    get_neighbors(os);

    #ifdef MRGRAPH_VERBOSE
    print_grid("About to process orphan (O) ", "o", os,'O');
    #endif


    uint8_t best_rank = uint8_t(-1);
    index_t best_id;
    dir_t best_dir;
    for (dir_t dir = 0; dir < DEGREE; ++dir) {
        if (label[neighbors[dir]] == LABEL_S &&
            rc_nbhd[neighbors[dir]][REVERSE[dir]] &&
            //parent_id[neighbors[dir]] != os &&
            dist[neighbors[dir]] < max_dist) {

            uint8_t cur_rank = ((get_origin(neighbors[dir]) == DIR_TERMINAL) ? 0 : 1) +
                                2*(dist[neighbors[dir]] - dist[os]) + 2;

            if (cur_rank < best_rank) {
                best_id = neighbors[dir];
                best_dir = REVERSE[dir];
                best_rank = cur_rank;
                if (best_rank == 0) {
                    break;
                }
            }
        }
    }
    //adopt
    if (best_rank == 0 || best_rank == 1){
        parent_id[os] = best_id;
        parent_edge[os] = best_dir;
        return;
    }
    //free
    else if (best_rank == uint8_t(-1)) {
        parent_id[os] = PARENT_ID_NONE;
        parent_edge[os] = DIR_NONE;
        label[os] = LABEL_F;
    }
    //relabel
    else {
        parent_id[os] = best_id;
        parent_edge[os] = best_dir;
        dist[os] = dist[best_id] + 1;
        if (dist [os] == s_dist) {
            active_s.push(os);
        }
        else if (dist[os] == s_dist+1) {
            found_s.push(os);
        }
    }

    for (dir_t dir = 0; dir < DEGREE; ++dir) {
        if(parent_id[neighbors[dir]] == os){
            parent_id[neighbors[dir]] = PARENT_ID_NONE;
            parent_edge[neighbors[dir]] = DIR_NONE;
            orphans_s.push(neighbors[dir]);
        }
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
rank_relabel_t(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t ot,
               const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t max_dist) {
    get_neighbors(ot);

    #ifdef MRGRAPH_VERBOSE
    print_grid("About to process orphan (O) ", "o", ot,'O');
    #endif

    uint8_t best_rank = uint8_t(-1);
    index_t best_id;
    dir_t best_dir;
    for (dir_t dir = 0; dir < DEGREE; ++dir) {
        if (label[neighbors[dir]] == LABEL_T &&
            rc_nbhd[ot][dir] &&
            //parent_id[neighbors[dir]] != ot &&
            dist[neighbors[dir]] < max_dist) {

            uint8_t cur_rank = ((get_origin(neighbors[dir]) == DIR_TERMINAL) ? 0 : 1) +
                                2*(dist[neighbors[dir]] - dist[ot]) + 2;
            if (cur_rank < best_rank) {
                best_id = neighbors[dir];
                best_dir = dir;
                best_rank = cur_rank;
                if (best_rank == 0) {
                    break;
                }
            }
        }
    }
    //adopt
    if (best_rank == 0 || best_rank == 1){
        parent_id[ot] = best_id;
        parent_edge[ot] = best_dir;
        return;
    }
    //free
    else if (best_rank == uint8_t(-1)) {
        parent_id[ot] = PARENT_ID_NONE;
        parent_edge[ot] = DIR_NONE;
        label[ot] = LABEL_F;
    }
    //relabel
    else {
        parent_id[ot] = best_id;
        parent_edge[ot] = best_dir;
        dist[ot] = dist[best_id] + 1;
        if (dist[ot] == t_dist) {
            active_t.push(ot);
        }
        else if (dist [ot] == t_dist+1) {
            found_t.push(ot);
        }
    }

    for (dir_t dir = 0; dir < DEGREE; ++dir) {
        if(parent_id[neighbors[dir]] == ot){
            parent_id[neighbors[dir]] = PARENT_ID_NONE;
            parent_edge[neighbors[dir]] = DIR_NONE;
            orphans_t.push(neighbors[dir]);
        }
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::dir_t MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
get_origin(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t v) {
    while (parent_edge[v] != DIR_NONE && parent_edge[v] != DIR_TERMINAL) {
        v = parent_id[v];
    }
    return parent_edge[v];
}

template <typename tcap_t,typename ncap_t,typename flow_t>
bool MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
grow_s(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t& vs,
       typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t& vt,
       typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::dir_t&   st) {
    while(!active_s.is_empty()) {
        index_t a = active_s.pop();
#ifdef MRGRAPH_VERBOSE
        print_grid("Growing s-node", "fA", a,'A');
#endif
        //active_s is what was previously found_s.
        //this means that it may contain "tainted" nodes that shouldnt be active s-nodes.
        //active nodes get tainted by being relabeled to DS+1 or getting freed.
        //found nodes get tainted by being freed.
        if (dist[a] != s_dist || label[a] != LABEL_S) {
#ifdef MRGRAPH_VERBOSE
            if (label[a] != LABEL_S){
                print_grid("Active node is tainted and ignored. (Wrong label)", "f", a,'!');
            }
            else {
                print_grid("Active node is tainted and ignored. (Wrong distance)", "f", a,'!');
            }
#endif
            continue;
        }
        //index_t neighbors[DEGREE];
        get_neighbors(a);
        for (dir_t dir = 0; dir < DEGREE; ++dir) {
            //scan all neighbors with residual connection from node a
            if(rc_nbhd[a][dir]) {
                if(label[neighbors[dir]] == LABEL_F) {
                    label[neighbors[dir]] = LABEL_S;
                    parent_id[neighbors[dir]] = a;
                    parent_edge[neighbors[dir]] = dir;
                    dist[neighbors[dir]] = dist[a] + 1;
                    found_s.push(neighbors[dir]);
                } else if (label[neighbors[dir]] == LABEL_T) {
                    vs = a;
                    vt = neighbors[dir];
                    st = dir;
                    active_s.push(a);
                    return true;
                }
            }
        }
    }
    return false;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
bool MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
grow_t(typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t& vs,
       typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t& vt,
       typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::dir_t&   st) {
    while(!active_t.is_empty()) {
        index_t a = active_t.pop();
#ifdef MRGRAPH_VERBOSE
        print_grid("Growing t-node", "fa", a,'A');
#endif
        if (dist[a] != t_dist || label[a] != LABEL_T) {
#ifdef MRGRAPH_VERBOSE
            if (label[a] != LABEL_T){
                print_grid("Active node is tainted and ignored. (Wrong label)", "f", a,'!');
            }
            else {
                print_grid("Active node is tainted and ignored. (Wrong distance)", "f", a,'!');
            }
#endif
            continue;
        }
        //index_t neighbors[DEGREE];
        get_neighbors(a);
        for (dir_t dir = 0; dir < DEGREE; ++dir) {
            if(rc_nbhd[neighbors[dir]][REVERSE[dir]]) {
                if(label[neighbors[dir]] == LABEL_F) {
                    label[neighbors[dir]] = LABEL_T;
                    parent_id[neighbors[dir]] = a;
                    parent_edge[neighbors[dir]] = REVERSE[dir];
                    dist[neighbors[dir]] = dist[a] + 1;
                    found_t.push(neighbors[dir]);
                } else if (label[neighbors[dir]] == LABEL_S) {
                    vt = a;
                    vs = neighbors[dir];
                    st = REVERSE[dir];
                    active_t.push(a);
                    return true;
                }
            }
        }
    }
    return false;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
get_neighbors(const typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t a) {
    neighbors[0] = north(a);
    neighbors[1] = east(a);
    neighbors[2] = south(a);
    neighbors[3] = west(a);
    for(int i = 4; i < DEGREE; ++i) {
        neighbors[i] = down(a,i-3);
    }
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
print() {
    std::cout
        <<"ow:   " << ow << std::endl
        <<"oh:   " << oh << std::endl
        <<"W:    " << W << std::endl
        <<"H:    " << H << std::endl
        <<"WH:   " << WH << std::endl
        <<"D:    " << D << std::endl
        <<"N:    " << N << std::endl
        <<"YOFS: " << YOFS << std::endl
        <<"XOFS: " << XOFS << std::endl
        <<"DEG:  " << int(DEGREE) << std::endl;
}

template <typename tcap_t,typename ncap_t,typename flow_t>
void MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
print_grid(std::string msg, std::string options = "",
           typename MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::index_t special_id = 0, char special_icon = char(32)) {
    #ifdef MRGRAPH_VERBOSE
    ++step;
    if (step < MRGRAPH_VERBOSE_START) {return;}
    std::cout << step << ": " << msg;
    if (special_id != 0) {
        std::cout << "\t (" << special_id << ")";
    }
    std::cout << std::endl;
    #endif
    bool draw_active_s = (options.find_first_of('A') != options.npos);
    bool draw_active_t = (options.find_first_of('a') != options.npos);
    bool draw_found = (options.find_first_of('f') != options.npos);
    bool draw_orphans = (options.find_first_of('o') != options.npos);
    bool draw_dist = (options.find_first_of('d') != options.npos);
    bool query_conn = (options.find_first_of('q') != options.npos);

    if (query_conn) {
        std::cout << "Which node? ";
        index_t queried_node;
        std::cin >> queried_node;
        std::cout << std::endl << "Residual outgoing edges of " << queried_node << ": ";
        (rc_nbhd[queried_node][0]) ? (std::cout << char(94)) : (std::cout << char(32));
        (rc_nbhd[queried_node][1]) ? (std::cout << char(62)) : (std::cout << char(32));
        (rc_nbhd[queried_node][2]) ? (std::cout << char(118)) : (std::cout << char(32));
        (rc_nbhd[queried_node][3]) ? (std::cout << char(60)) : (std::cout << char(32));
        std::cout << std::endl << "Dist: " << dist[queried_node] << std::endl;
        #ifdef MRGRAPH_VERBOSE
        std::string input;
        std::getline(std::cin,input);
        if (input.size() != 0){
            std::cout << "Input: ";
            --step;
            print_grid(input, input);
        }
        #endif
        return;
    }

    if (draw_dist) {
        for( int y = 0; y < oh; ++y) {
            for (int x = 0; x < ow; ++x) {
                index_t v = node_id(x,y,0);
                if (label[v] != LABEL_F) {
                    std::cout << dist[v];
                } else {
                    std::cout << char(176);
                }
            }
            std::cout << std::endl;
        }
        std::cout << "s_dist: " << s_dist << "\tt_dist: " << t_dist << std::endl;
        #ifdef MRGRAPH_VERBOSE
        std::string input;
        std::getline(std::cin,input);
        if (input.size() != 0){
            std::cout << "Input: ";
            --step;
            print_grid(input, input);
        }
        #endif
        return;
    }

    for( int y = 0; y < oh; ++y) {
        for (int x = 0; x < ow; ++x) {
            index_t v = node_id(x,y,0);
            if (v == special_id) {
                std::cout << special_icon;
                continue;
            }
            bool found_special = false;
            if (draw_active_s) {

                for (int i = 0; i < active_s.size; ++i) {
                    if (active_s.buffer[i] == v) {
                        std::cout << 'a';
                        found_special = true;
                    }
                }
                if (found_special) {
                    continue;
                }
            }

            if (draw_active_t) {
                for (int i = 0; i < active_t.size; ++i) {
                    if (active_t.buffer[i] == v) {
                        std::cout << 'a';
                        found_special = true;
                    }
                }
                if (found_special) {
                    continue;
                }
            }


            if (draw_found) {
                if (label[v] == LABEL_S) {
                    for (int i = 0; i < found_s.size; ++i) {
                        if (found_s.buffer[i] == v) {
                            std::cout << 'f';
                            found_special = true;
                        }
                    }
                    if (found_special) {
                        continue;
                    }
                } else {
                    for (int i = 0; i < found_t.size; ++i) {
                        if (found_t.buffer[i] == v) {
                            std::cout << 'f';
                            found_special = true;
                        }
                    }
                    if (found_special) {
                        continue;
                    }
                }
            }

            if (draw_orphans) {
                if (label[v] == LABEL_S) {
                    for (int i = 0; i < orphans_s.size; ++i) {
                        if (orphans_s.buffer[i] == v) {
                            std::cout << 'o';
                            found_special = true;
                        }
                    }
                    if (found_special) {
                        continue;
                    }
                } else {
                    for (int i = 0; i < orphans_t.size; ++i) {
                        if (orphans_t.buffer[i] == v) {
                            std::cout << 'o';
                            found_special = true;
                        }
                    }
                    if (found_special) {
                        continue;
                    }
                }
            }

            uint8_t conn = 0;
            index_t nbhrs[4];
            nbhrs[0] = north(v);
            nbhrs[1] = east(v);
            nbhrs[2] = south(v);
            nbhrs[3] = west(v);
            //get_neighbors(v,neighbors);
            if (label[v] == LABEL_S) {
                if (parent_edge[v] == DIR_TERMINAL) {
                    std::cout << 's';
                    continue;
                }
                dir_t par = REVERSE[parent_edge[v]];
                conn += (1<<par);
                for (dir_t dir = 0; dir < 4; ++dir) {
                    if (parent_id[nbhrs[dir]] == v) {
                        conn += (1<<dir);
                    }
                }
            } else if (label[v] == LABEL_T) {
                if (parent_edge[v] == DIR_TERMINAL) {
                    std::cout << 't';
                    continue;
                }
                dir_t par = parent_edge[v];
                conn += (1<<par);
                for (dir_t dir = 0; dir < 4; ++dir) {
                    if (parent_id[nbhrs[dir]] == v) {
                        conn += (1<<dir);
                    }
                }
            } else {
                std::cout << char(176);
                continue;
            }

            std::cout << get_icon(conn);
        }
        std::cout << std::endl;
    }
    #ifdef MRGRAPH_VERBOSE
    std::string input;
    std::getline(std::cin,input);
    if (input.size() != 0){
        std::cout << "Input: ";
        --step;
        print_grid(input, input);
    }
    #endif
}


template <typename tcap_t,typename ncap_t,typename flow_t>
char MRGraph_2D_4C<tcap_t,ncap_t,flow_t>::
get_icon(uint8_t conn) {
    switch(conn) {
    case  0:
        return '?';
    case  1:
        return char(245);//179
    case  2:
        return char(169);//196
    case  3:
        return char(192);
    case  4:
        return char(244);//179
    case  5:
        return char(179);
    case  6:
        return char(218);
    case  7:
        return char(195);
    case  8:
        return char(170);//196
    case  9:
        return char(217);
    case 10:
        return char(196);
    case 11:
        return char(193);
    case 12:
        return char(191);
    case 13:
        return char(180);
    case 14:
        return char(194);
    case 15:
        return char(197);
    }

    return '!';
}


#endif // MRGRAPH_2D_4C_H
