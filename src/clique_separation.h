#ifndef CLIQUE_SEPARATION_H_INCLUDED
#define CLIQUE_SEPARATION_H_INCLUDED

typedef struct _CliqueSeparation CliqueSeparation;

#include "cgraph.h"
#include "clique.h"
#include "clique_extender.h"

/* origGraph is the original conflict graph,
   containing all nodes. separation will always
   process a subset of these nodes */
CliqueSeparation *clq_sep_create( const CGraph *origGraph );

/* separates clique inequalities for fractional solution x */
void clq_sep_separate( CliqueSeparation *sep, const double x[] );

/* if clique will be extended using reduced costs, an array for
   these values should be informed before each clq_sep_separate call */
void clq_sep_set_rc( CliqueSeparation *sep, const double rc[] );

void clq_sep_set_verbose( CliqueSeparation *sep, const char verbose );

int clq_sep_get_verbose( CliqueSeparation *sep );

void clq_sep_set_min_viol( CliqueSeparation *sep, const double viol );

double clq_sep_get_min_viol( CliqueSeparation *sep );

void clq_sep_set_extend_method( CliqueSeparation *sep, const int extendC );

int clq_sep_get_extend_method( CliqueSeparation *sep );

/* returns separated clique inequalities */
const CliqueSet *clq_sep_get_cliques( CliqueSeparation *sep );

void clq_sep_free( CliqueSeparation **clqSep );

void clq_sep_params_print( const CliqueSeparation *clqsp );

void clq_sep_set_params_parse_cmd_line( CliqueSeparation *clqsp, const int argc, const char **argv );

void clq_sep_params_help_cmd_line();

int clq_sep_get_max_passes( const CliqueSeparation *clqSep );

double clq_sep_get_max_time_bk( const CliqueSeparation *clqSep );

void clq_sep_set_max_time_bk( CliqueSeparation *clqSep, const double maxTime );

/* "private" functions, do not use unless you know what you are doing */

/* enumerate all cliques in cgraph for containing nodes with low degree (param enumUsage),
   these nodes will receive value 1 in iv */
void clq_set_enum_all_cliques_low_degree_nodes( CliqueSeparation *sep, int iv[], const CGraph *origGraph, const CGraph *ppGraph );

#endif

