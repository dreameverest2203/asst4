#include "bfs.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <numeric>
#include <execution>
#include <algorithm>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1
#define CHUNK 500

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
#pragma omp parallel 
{
    const int num_threads = omp_get_num_threads();
    const int thread_id = omp_get_thread_num();
    vertex_set local_frontier;
    vertex_set_init(&local_frontier, g->num_nodes);
    #pragma omp for
    for (int i=0; i<frontier->count; i++) {
        
        int node = frontier->vertices[i];
        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                            ? g->num_edges
                            : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];
            
            if (distances[outgoing] == NOT_VISITED_MARKER && __sync_bool_compare_and_swap(distances + outgoing, NOT_VISITED_MARKER, distances[node] + 1)) {
                local_frontier.vertices[local_frontier.count++] = outgoing;
            }
        }
    }
    #pragma omp critical 
        {
            memcpy(new_frontier->vertices + new_frontier->count, local_frontier.vertices, sizeof(int)*local_frontier.count);
            new_frontier->count += local_frontier.count;
        }
    vertex_set_clear(&local_frontier); 
}
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    #pragma omp for
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
    
}

// bool check(Graph g, vertex_set* frontier, int vertex) {
//     for(int i = 0; i < frontier->count; i++) {
//         int node = frontier->vertices[i];
//         int start_edge = g->incoming_starts[node];
//         int end_edge = (node == g->num_nodes - 1)
//                             ? g->num_edges
//                             : g->incoming_starts[node + 1];
//         for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
//             if(g->incoming_edges[neighbor]==vertex) return true;
//         }
//     }
//     return false;
// }

// void bottom_up_step(
//     Graph g,
//     vertex_set* frontier,
//     vertex_set* new_frontier,
//     int* distances)
// {
//     for(int i = 0; i <g->num_nodes; i++) {
//         if(distances[i]==NOT_VISITED_MARKER && check(g, frontier, i)) {
//             new_frontier->vertices[new_frontier->count++] = i;
//         }
//     }
// }



void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances,
    bool* old_frontier_bool)
{
#pragma omp parallel
{
    const int n_threads = omp_get_num_threads();
    int new_dist = distances[frontier->vertices[0]] + 1;
    int* per_thread_storage = (int*)malloc(sizeof(int)* g->num_nodes);
    int nCount = 0;
    #pragma omp for
    for (int node=0; node<g->num_nodes; node++) {
        if(old_frontier_bool[node]){
            int start_edge = g->incoming_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->incoming_starts[node + 1];
     
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                if(!old_frontier_bool[g->incoming_edges[neighbor]]){
                    distances[node] = new_dist;
                    per_thread_storage[nCount] = node;
                    nCount++;
                    break;
                }
            }
        }
    }
    #pragma omp critical
    {
        memcpy(new_frontier->vertices + new_frontier->count, per_thread_storage, sizeof(int)*nCount);
        new_frontier->count += nCount;
    }
    free(per_thread_storage);
}
}



void bfs_bottom_up(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    bool* old_frontier_bool = (bool*)malloc(sizeof(bool) * graph->num_nodes);

    
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
        old_frontier_bool[i] = true;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    old_frontier_bool[ROOT_NODE_ID] = false;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step(graph, frontier, new_frontier, sol->distances, old_frontier_bool);

        #pragma omp parallel for
        for (int i=0; i<new_frontier->count; i++){
            old_frontier_bool[new_frontier->vertices[i]] = false;
        }

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }

    free(old_frontier_bool);
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    bool* node_unvisited = (bool*)malloc(sizeof(bool) * graph->num_nodes);

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
        node_unvisited[i] = true;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;
    node_unvisited[ROOT_NODE_ID] = false;
    int remaining_nodes = graph->num_nodes - 1;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        if(remaining_nodes*0.1 < frontier->count) bottom_up_step(graph, frontier, new_frontier, sol->distances, node_unvisited);
        else top_down_step(graph, frontier, new_frontier, sol->distances);    

        #pragma omp parallel for
        for (int i=0; i<new_frontier->count; i++){
            node_unvisited[new_frontier->vertices[i]] = false;
        }

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        remaining_nodes -= new_frontier->count;
    }

    free(node_unvisited);

}
