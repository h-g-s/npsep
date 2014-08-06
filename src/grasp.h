#ifndef GRASP_H_INCLUDED
#define GRASP_H_INCLUDED

#include "cgraph.h"
#include "clique.h"

typedef struct _Grasp Grasp;


/**
 * creates with default parameters
 * considering conflict graph _cgraph
 * and weights w
 **/
Grasp *grasp_create( const CGraph *_cgraph, const int _minW );

/**
 * continuous parameter [0 ... 1] shuch that
 *    0 - greedy
 *    1 - random
 **/
void grasp_set_alpha( Grasp *grasp, const double _alpha );

/* how to chose alpha, alpha coice: */
#define GRASP_ALPHA_FIXED    0
#define GRASP_ALPHA_RANDOM   1
#define GRASP_ALPHA_REACTIVE 2
void grasp_set_alpha_choice( Grasp *grasp, const int choice );

void grasp_set_max_no_improvement( Grasp *grasp, const int _max_ni_it );

/**
 * run
 **/
void grasp_run( Grasp *grasp );

/**
 * gets all stored solutions
 **/
CliqueSet *grasp_solution_set( Grasp *grasp );

/**
 * return best/worst weight
 **/
int grasp_get_best_weight( Grasp *grasp );

/**
 * returns workst weight
 **/
int grasp_get_worst_weight( Grasp *grasp );

/**
 * frees memory
 **/
void grasp_free( Grasp **grasp );

/* internal use functions */

void grasp_fill_nodes_left( Grasp *grasp );

void grasp_build_candidate_list( Grasp *grasp );

int grasp_select_from_candidate_list( const Grasp *grasp );

/* after adding to the clique node with "neighbors", update nodes left */
void grasp_update_nodes_left( Grasp *grasp, const int n, const int neighbors[], const int nodeToEnter );

void grasp_iteration( Grasp *grasp );

/* "private" functions */
void grasp_ini_reactive_info( Grasp *grasp );
/* if random or reactive, updates alpha */
void grasp_select_alpha( Grasp *grasp );
/* after some interations, update reactive probabilities */
void grasp_recompute_reactive_probabilities( Grasp *grasp );
/* prepare roullete considering probabilibies */
void grasp_prepare_roullete( Grasp *grasp );
/* prepare initial probabilities */
void grasp_fill_initial_probabilities( Grasp *grasp );

#endif

