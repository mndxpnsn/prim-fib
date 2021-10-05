/*
 * main.cpp
 *
 *  Created on: Aug 9, 2020
 *      Author: d-w-h
 *
 *      Implementation of Prim's algorithm
 *      using a fibonacci heap.
 */

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "../inc/user_types.hpp"
#include "../inc/functions.hpp"

int main(int argc, char* argv[])
{
    //Declarations
    int s = 2; //Start vertex. The minimum index for vertices is 1
    int n = 2499; //Number of vertices
    int num_edges = 3125; //Number of edges

    //Create edges
    srand(time(NULL));
    std::vector< edge > edges;
    for(int i = 0; i < num_edges; ++i) {
        int start_vert = rand() % n + 1;
        int end_vert = rand() % n + 1;
        float weight =  (( (float) rand() ) / RAND_MAX + 1.0) * 3.14159;

        edge edge_elem;
        edge_elem.start_vertex = start_vert;
        edge_elem.end_vertex = end_vert;
        edge_elem.weight = weight;
        edges.push_back(edge_elem);
    }

    //Compute minimum spanning tree
    mst_props min_span_props = mst(n, edges, s);

    //Print results
    float tot_num_ops_est = 8 * n + 3 * num_edges + 5 * n * log(n) / log(2);
    print_mst(n, min_span_props.node_arr);
    std::cout << "size of minimum spanning tree: " << min_span_props.mst_weight << std::endl;
    std::cout << "total number of operations measured: " << tot_num_ops << std::endl;
    std::cout << "estimated number of operations 8V + 3E + 5VlgV: " << tot_num_ops_est << std::endl;
    std::cout << "ratio complexities: " << (tot_num_ops / tot_num_ops_est) << std::endl;
    std::cout << "done" << std::endl;

    return 0;
}
