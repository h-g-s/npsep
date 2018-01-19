/** Conflict Graph
    prepared to store conflict graphs
    and also large cliques

    Haroldo Gambini Santos - haroldo@iceb.ufop.br
    2010
**/

#ifndef CGRAPH_H_INCLUDED
#define CGRAPH_H_INCLUDED

typedef struct _CGraph CGraph;

typedef struct _NeighIterator NeighIterator;

/**
 * starts a new conflict graph prepared to work
 * initially with "columns" columns
 **/
CGraph *cgraph_create( int initialColumns );

CGraph* cgraph_clone(const CGraph *cg);

/**
 * after filling conflicts, call this function
 * to update min and max degrees
 **/
void cgraph_update_min_max_degree( CGraph *cgraph );

/**
 * in graphs with large cliques
 * degree can be overestimated */
void cgraph_recompute_degree( CGraph *cgraph );

/***
 * creates a preprocessed conflict graph,
 * which will have a subset of nodes
 * nindexes indicates for a given node i
 * the new index in the preprocessed
 * conflict graph
 **/
CGraph *cgraph_preprocess( const CGraph *cgraph, int nindexes[] );

/**
 * creates as subgraph induced by vertices
 * in vector nindexes some indices i will have nindexes[i] == -1
 * these nodes will not be included in the new subgraph
 **/
CGraph *cgraph_create_induced_subgraph( const CGraph *cgraph, const int nindexes[] );

/**
 * adds all nodes in conflicts
 * as conflicting nodes for node
 * if necessary, increases the number of columns
 **/
void cgraph_add_node_conflicts( CGraph *cgraph, const int node, const int conflicts[], const int size );

int cgraph_degree( const CGraph *cgraph, const int node );

int cgraph_min_degree( const CGraph *cgraph );

int cgraph_max_degree( const CGraph *cgraph );

int cgraph_get_original_node_index( const CGraph *cgraph, const int node );

const int *cgraph_get_original_node_indexes( const CGraph *cgraph );

/**
 * adds all nodes in conflicts
 * as conflicting nodes for node
 * if necessary, increases the number of columns
 * this version is faster because of there is a conflict (i,j) it does not cares for inserting (j,i)
 **/
void cgraph_add_node_conflicts_no_sim( CGraph *cgraph, const int node, const int conflicts[], const int size );

/**
 * adds a clique of conflicts
 **/
void cgraph_add_clique( CGraph *cgraph, int clique[], const int size );

/**
 * answers if two nodes are
 * conflicting
 **/
int cgraph_conflicting_nodes( const CGraph *cgraph, const int i, const int j );

/**
 * fills all conflicting nodes
 * returns the number of conflicting nodes or
 * aborts if more than maxSize neighs are found
 **/
int cgraph_get_all_conflicting( const CGraph *cgraph, int node, int neighs[], int maxSize );

/**
 * fills first n conflicting nodes
 *
 * v is a temporary vector for processing and
 * vcap is its capacity
 * if vcap is not enough aborts with error
 */
int cgraph_get_n_conflicting( const CGraph *cgraph, int node, int neighs[], int n,
                              int v[], const int vcap );

/**
 * returns the size of the conflict graph
 **/
int cgraph_size( const CGraph *cgraph );

/**
 * creates a new cgraph considering content from a file
 * if w != null reads weights
 **/
CGraph *cgraph_load( const char *fileName );

/**
 * read dimensions from file
 **/
void cgraph_load_dimensions( const char *fileName, int *nodes, int *edges );

/* node names and weights are optional information */

int cgraph_get_node_weight( const CGraph *cgraph, int node );

const int *cgraph_get_node_weights( const CGraph *cgraph );

void cgraph_set_node_weight( CGraph *cgraph, int node, int weight );

const char *cgraph_get_node_name( const CGraph *cgraph, int node );

void cgraph_set_node_name( CGraph *cgraph, int node, const char *name );

/**
 * converts a double weight to an integer weight
 * */
int cgraph_weight( const double w );

/**
 * save the graph in the file
 **/
void cgraph_save( CGraph *cgraph, const char *fileName );

/**
 * for debugging purposes
 **/
void cgraph_print( CGraph *cgraph, const int w[] );

void cgraph_print_summary( CGraph *cgraph, const char *name );

void cgraph_set_low_degree( CGraph *cgraph, const int lowDegree );

/**
 * releases all conflict graph memory
 **/
void cgraph_free( CGraph **cgraph );

/* creates an object which will be used to
 * iterate through neighbors of nodes
 * returning first those with the smallest
 * degree
 */
NeighIterator *nit_create( );

/* starts iterating in a node */
void nit_start( NeighIterator *nit, const CGraph *cgraph, int node, const int costs[] );

/* returns next neighbor or
 * INT_MAX if this is the end of
 * the list */
int nit_next( NeighIterator *nit );

void nit_free( NeighIterator **nit );

/**
 * DEBUG
 **/
#ifdef DEBUG
void cgraph_check_node_cliques( const CGraph *cgraph );
void cgraph_check_neighs( const CGraph *cgraph );
void cgraph_check_preproc( const CGraph *ppgraph, const CGraph *cgraph );
#endif

#endif

