/*
 * functions.hpp
 *
 *  Created on: Aug 13, 2020
 *      Author: d-w-h
 */

#ifndef FUNCTIONS_HPP_
#define FUNCTIONS_HPP_

#include "user_types.hpp"

void free_node_ref(node** v_ref, int size);
void print_mst(int size_heap, node** node_arr);
mst_props mst(int n, std::vector<edge>& edges, int s);

#endif /* FUNCTIONS_HPP_ */
