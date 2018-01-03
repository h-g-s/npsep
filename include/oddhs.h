#ifndef ODD_HOLE_SEP_H
#define ODD_HOLE_SEP_H

#include "cgraph.h"
#include "clique.h"

#define ODDH_SEP_DEF_MIN_VIOL 0.02

typedef struct _OddHoleSep OddHoleSep;
typedef OddHoleSep * OddHoleSepPtr;

/**
 * Initialize Odd-Holes
 * search object.
 **/
OddHoleSep *oddhs_create( );

/**
 * finds Odd Holes (OHs), which preferably
 * (but not necessarily), correspond to
 * violated cuts (VCs). those which are
 * not violated are also kept since they
 * can be expanded (or transformed in larger
 * cliques in the case of OHs of size 3).
 * returns the number of odd holes found
 * aggressiveness indicates how much processing
 * will be dedicate to the search:
 * - 0 will perform a basic search
 * - 1 more aggressive search
 **/
int oddhs_search_odd_holes( OddHoleSep *oddhs, const int cols, const double x[],
        const double rc[], const CGraph *conflicts );

/**
 * computes the left hand side of
 * the constraint which will be created by the odd hole
 ***/
double oddhs_lhs( const int size, const int idx[], const double x[] );

/**
 * returns the rhs of a oddhs constraint of size s
 * (s/2)
 **/
double oddhs_rhs( const int size );

/**
 * computes the violation of a odd hole
 **/
double oddhs_viol( const int size, const int idx[], const double x[] );

/**
 * Returns the start of the idx-th discovered odd hole
 * in the last search.
 * All odd-holes are in sequence, so that the end of the idx-the
 * odd hole is the start of the (idx+1)-th.
 * Indexes are related to the original indexes of variables.
 **/
int *oddhs_get_odd_hole( OddHoleSep *oddhs, const int idx );


/**
 * Returns the total number of discovered odd holes in
 * the last search.
 **/
int oddhs_get_odd_hole_count( OddHoleSep *oddhs );

/**
 * all odd holes of size 3 are considered to be cliques
 **/
const CliqueSet *oddhs_get_cliques( OddHoleSep *oddhs );

/**
 * the inequality for a discovered odd hole
 * may be extended with the addition of
 * wheel centers - this function returns the
 * number of computed wheel centers for a
 * discovered odd hole
 */
int oddhs_get_nwc_doh( OddHoleSep *oddhs, const int doh );


/**
 * the inequality for a discovered odd hole
 * may be extended with the addition of
 * wheel centers - this function returns the
 * computed wheel centres for a discovered odd hole
 */
const int *oddhs_get_wc_doh( OddHoleSep *oddhs, const int doh );

/**
 * Frees memory related to Odd-Hole
 * search object.
 **/
void oddhs_free( OddHoleSepPtr *oddhs );

#endif
