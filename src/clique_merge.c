#include <ctype.h>
#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "clique_merge.h"
#include "macros.h"
#include "clique_extender.h"
#include "vint_set.h"
#include "memory.h"

enum CliqueType
{
    NotAClique = 0,
    NotDominated = 1,
    Dominated = 2,
    MoreThanAClique = 3 // do not try to dominate this constraint
};

#define MAX_SIZE_CLIQUE_TO_BE_EXTENDED 256

// global configurations and stats
int clqMergeVerbose = 1;
double clqMergeSecsCheckClique = 0.0;
double clqMergeSecsExtendAndDominate = 0.0;
double clqMergeSecsAddAndRemove = 0.0;
double clqMergeSecsExtend = 0.0;
int clqMergeNExtended = 0;
int clqMergeNDominatedFull = 0;
int clqMergeNDominatedEqPart = 0;

static char dominates( int nClq1, const int clq1[], int nClq2, const int clq2[], char *iv )
{
    if ( nClq1 < nClq2 )
        return False;

    // filling incidence vector
    for ( int i=0 ; (i<nClq1) ; ++i )
    {
#ifdef DEBUG
        assert( iv[clq1[i]]==False );
#endif
        iv[clq1[i]] = True;
    }

    char res = True;
    // checking if clq2 is contained in clq1
    for ( int i=0 ; (i<nClq2) ; ++i )
    {
        if (iv[clq2[i]]==False)
        {
            res = False;
            break;
        }
    }


    // clearing iv
    for ( int i=0 ; (i<nClq1) ; ++i )
        iv[clq1[i]] = False;
    return res;
}

static void add_clique(
        LinearProgram *mip,
        CliqueSet *newCliques,      // new cliques
        int size, const int el[],  // clique being added
        enum CliqueType cliqueState[], const int clqRow[],  // which are the clique rows and their state
        const int **elRows, const int sizeRows[],  // cliques in original problem
        char *iv, // temporary incidence vector
        char *ivrt, // temporary incidence vector
        int *rtc, // list of rows checked
        const int **colClqs, // column cliques
        const int nColClqs[] // number of cliques of a column
        )
{
#ifdef DEBUG
    for ( int i=0 ; (i<lp_cols(mip)) ; ++i )
    {
        assert( iv[i] == False );
    }
    for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
    {
        assert( ivrt[i] == False );
    }

#endif
    int add = 0;
#pragma omp critical
    {
        int add = 
            clq_set_add( newCliques, size, el, size );
    }
    if (add==0)
        return;


    int nrtc = 0;
    for ( int i=0 ; (i<size) ; ++i )
    {
        int col = el[i];
        if (col >= lp_cols(mip)) // complimetary variable
            col -= lp_cols(mip);
#ifdef DEBUG            
        assert( col >= 0 && col<lp_cols(mip) );
#endif        
        // checking non-dominated cliques that column appear
        for ( int j=0 ; (j<nColClqs[col]) ; ++j )
        {
            int ir = colClqs[col][j];
#ifdef DEBUG
            assert( cliqueState[ir]!=NotAClique );
#endif
            
            // skipping already dominated, already included or larger rows (cannot be dominated by this one )
            if ( cliqueState[ir]==Dominated || ivrt[ir] || sizeRows[ir]>size )
                continue;

            ivrt[ir] = True;
            rtc[nrtc++] = ir;;
        }
    }

    for ( int i=0 ; (i<nrtc) ; ++i )
    {
        int rowClique = rtc[i];
        ivrt[rowClique] = False;

        if (dominates(size, el, sizeRows[rowClique], elRows[rowClique], iv ))
        {
#pragma omp critical
{
            cliqueState[rowClique] = Dominated;
            if (clqMergeVerbose>=2)
            {
                char domname[256]="";
                lp_row_name( mip, rowClique, domname );
                printf("\t\tdominates %s\n", domname );
            }
            
} // critical
        } // dominates 
    } // all rows to check
 
}

static int check_cliques( LinearProgram *mip, enum CliqueType cliqueState[], int cliques[], int **elRow, int *sizeRow, int *allEl, int *nAllEl, const char *_sense )
{
    clock_t startcq = clock();
    if (clqMergeVerbose>=1)
        printf( "checking candidates for clique extension.\n" );
    int rows = lp_rows(mip);
    int nCliques = 0;
    *nAllEl = 0;

    int *idx;
    ALLOCATE_VECTOR( idx, int, lp_cols(mip) );
    double *coef;
    ALLOCATE_VECTOR( coef, double, lp_cols(mip) );
    
    int minClqSize = INT_MAX;
    int maxClqSize = 0;
    double avClgSize = 0.0;

    for ( int i=0 ; (i<rows) ; ++i )
    {
        int nz = lp_row( mip, i, idx, coef );
        if ( nz<=1 || nz>MAX_SIZE_CLIQUE_TO_BE_EXTENDED )
            continue;
        char sense = _sense[i];
        double rhs = lp_rhs( mip, i );
        // all constraints as <=
        double mult = (sense == 'G') ? -1.0 : 1.0;

        /* all variables should be positive, integers */
        char varsOk = True;
        double minCoef = DBL_MAX;

        // to check if original and complimentary variables are involved
        int nOnes = 0, nMinusOne = 0;

        for ( int j=0 ; (j<nz) ; ++j )
        {
            if ( !lp_is_integer(mip, idx[j]) || lp_col_lb(mip,idx[j]) <= -1e-5 || ( coef[j]<=-1e8 && lp_col_ub(mip,idx[j])>=1.0+1e-8 ) )
            {
                varsOk = False;
                break;
            }

            // binaries, checking for original and complimentary variables
            if ( lp_col_lb( mip, idx[j] ) >= -1e-8 && lp_col_ub( mip, idx[j] ) <= 1.0+1e-8 )
            {
                if ( fabs(mult*coef[j]-1.0)<=1e-8 ) 
                    ++nOnes;
                else
                    if ( fabs(mult*coef[j]+1.0)<=1e-8 )
                        ++nMinusOne;
            }

            minCoef = MIN( mult * coef[j], minCoef );
        }

        if (!varsOk)
            continue;

        // first case
        cliqueState[i] = ( 2*minCoef >= mult*rhs+1e-5 && rhs >= 1e-5 ) ? NotDominated : NotAClique;

        if (cliqueState[i]==NotDominated)
        {
            cliques[nCliques++] = i;

            elRow[i] = allEl;
            sizeRow[i] = nz;
            memcpy( elRow[i], idx, sizeof(int)*nz );
            allEl += nz;
            
            minClqSize = MIN( minClqSize, nz );
            maxClqSize = MAX( maxClqSize, nz );
            avClgSize += nz;
        } // found a clique candidate
        else
        {
            // checking for clique involving normal variables and complementary variables
            if ( (nOnes+nMinusOne==nz) && (fabs(rhs-(1.0-nMinusOne))<=1e-8) )
            {
                cliqueState[i]=NotDominated;
                
                cliques[nCliques++] = i;

                elRow[i] = allEl;
                sizeRow[i] = nz;
                memcpy( elRow[i], idx, sizeof(int)*nz );
                for ( int j=0 ; (j<nz) ; ++j )
                {
                    if (coef[j]<-1e-8)
                        idx[j] += lp_cols(mip);
                }
                allEl += nz;
                minClqSize = MIN( minClqSize, nz );
                maxClqSize = MAX( maxClqSize, nz );
                avClgSize += nz;
            }
        }

        if ( clqMergeVerbose>=2 && cliqueState[i]==NotDominated )
        {
            char rname[256]="";
            lp_row_name( mip, i, rname );
            printf( "\trow %s: ", rname );
            for ( int j=0 ; (j<nz) ; ++j )
            {
                char cname[256]="";
                lp_col_name( mip, idx[j], cname );
                printf(" %+g %s ", coef[j], cname );
            }
            char strSense[3] = "";
            switch (sense)
            {
                case 'E' :
                {
                    strcpy( strSense, "=") ;
                    break;
                }
                case 'L' :
                {
                    strcpy( strSense, "<=") ;
                    break;
                }
                case 'G' :
                {
                    strcpy( strSense, ">=") ;
                    break;
                }
            } // sense

            printf("%s %g\n", strSense, rhs );
        } // if verbose

    } // all rows

    free( idx );
    free( coef );

    clqMergeSecsCheckClique = ((double)clock()-startcq) / ((double)CLOCKS_PER_SEC);

    if (clqMergeVerbose>=1)
        printf("model checked in %.4f seconds. %d candidate cliques for extension/merging. clique sizes range:[%d...%d], av %.2f.\n", clqMergeSecsCheckClique, nCliques, minClqSize, maxClqSize, avClgSize/((double)nCliques) );
        

    return nCliques;
}

/* to sort cliques per size */
struct CliqueSize
{
    int clqIdx;
    int size;
};

static int cmp_clq_size( const void *v1, const void *v2 )
{
    const struct CliqueSize *cs1 = (const struct CliqueSize *) v1;
    const struct CliqueSize *cs2 = (const struct CliqueSize *) v2;

    return cs2->size-cs1->size;
}

static void addName( char ***rNames, int *lines, int *chars, int *linesCap, int *charsCap, const char *name  )
{
    int len = strlen(name);

    if ( (*chars) + len+2 >= (*charsCap) )
    {
        // line sizes must be stored to fix the pointers
        int *lsizes;
        ALLOCATE_VECTOR( lsizes, int, (*lines) );
        for ( int i=0 ; i<(*lines) ; ++i )
            lsizes[i] = strlen( (*rNames)[i] );

        (*charsCap) = MAX((*charsCap)*2, (*chars) + len+2);

        (*rNames)[0] = xrealloc( (*rNames)[0], sizeof(char)*(*charsCap) );
        for ( int i=1 ; i<(*lines) ; ++i )
            (*rNames)[i] = (*rNames)[i-1] + lsizes[i-1] + 1;

        free( lsizes );
    }

    if ( (*lines) + 1 >= (*linesCap) )
    {
        *linesCap = MAX( 2*(*linesCap), (*linesCap)+10  );
        *rNames = xrealloc( *rNames, sizeof(char*)*(*linesCap) );
    }

    if ((*lines)>0)
        (*rNames)[*lines] = (*rNames)[(*lines)-1] + strlen((*rNames)[(*lines)-1]) + 1;

    strcpy( (*rNames)[*lines], name );

    ++(*lines);
    (*chars) += len+1;
}

void addRow( 
        int *nrRows, int *nrCapRows, int *nrNz, int *nrCapNz, int **nrStart, int **nrIdx, double **nrCoef, char **nrSense, double **nrRhs,
        int nz, const int idx[], const double coef[], char sense, double rhs,
        const char *rowName, char ***rowNames, int *nrLines, int *nrCap, int *nrChars, int *nrCharsCap )
{
    if ( (*nrRows)+2 >= (*nrCapRows) )
    {
        (*nrCapRows) = MAX( 2*(*nrCapRows),  (*nrCapRows)+32 );

        (*nrStart) = xrealloc( (*nrStart), sizeof(int)*((*nrCapRows)+1) );
        (*nrSense) = xrealloc( (*nrSense), sizeof(char)*((*nrCapRows)) );
        (*nrRhs) = xrealloc( (*nrRhs), sizeof(double)*(*nrCapRows) );
    }

    if ( ((*nrNz)+nz) >= (*nrCapNz))
    {
        (*nrCapNz) = MAX( (*nrCapNz)*2, ((*nrNz)+nz) );
        (*nrIdx) = xrealloc( (*nrIdx), sizeof(int)*(*nrCapNz) );
        (*nrCoef) = xrealloc( (*nrCoef), sizeof(double)*(*nrCapNz) );
    }

    /* adding */
    int *pidx = (*nrIdx) + (*nrNz);
    double *pcoef = (*nrCoef) + (*nrNz);
    memcpy( pidx, idx, sizeof(int)*nz );
    memcpy( pcoef, coef, sizeof(double)*nz );
    (*nrSense)[(*nrRows)] = sense;
    (*nrRhs)[(*nrRows)] = rhs;

    (*nrStart)[(*nrRows)+1] = (*nrStart)[(*nrRows)] + nz;
    (*nrNz) += nz;
    
    addName( rowNames, nrLines, nrChars, nrCap, nrCharsCap, rowName );
    
    ++(*nrRows);
}

/* tries to extend every clique in mip using
 * conflict graph cgraph, dominated cliques are removed */
void merge_cliques( LinearProgram *mip, CGraph *cgraph, int maxExtensions, int maxItBk )
{
    enum CliqueType *cliqueState; // if it is a clique or not and if it is dominated
    int *cliques;       // list of rows which are cliques
    int *rc = NULL;               // reduced costs (dummy values by now)

    int **idx = NULL;              // temporary area
    double **coef = NULL;          // to get indexes from rows

    ALLOCATE_VECTOR( cliques, int, lp_rows(mip) );
    ALLOCATE_VECTOR( cliqueState, enum CliqueType, lp_rows(mip) );
    FILL( cliqueState, 0, lp_rows(mip), NotAClique );

    CliqueSet *newCliques = clq_set_create(); // where new cliques will be stored

    int **elRow = NULL;  // clique elements per row
    int *sizeRow = NULL; // size of clique row
    int *allEl = NULL;   // to store all non-zeros
    int nAllEl = 0;      // total number of non zeros in cliques

    int *nColClqs = NULL;    // number of cliques that a column is involved
    int *allCollClqs = NULL; // contiguous vector to store column cliques
    int **colClqs = NULL;    // cliques that a column is involved

    ALLOCATE_VECTOR_INI( elRow, int *, lp_rows(mip) );
    ALLOCATE_VECTOR_INI( sizeRow, int, lp_rows(mip) );
    ALLOCATE_VECTOR( allEl, int, lp_nz(mip) );

    /* working variables of each thread */
    char **ivt = NULL;      // incidence vector per thread
    char **ivrt = NULL;     // incidence vector per thread
    int **rtc = NULL;       // rows to check per thread
    IntSet *currClique = NULL;  // to store the current clique
    struct CliqueSize **clqsSize = NULL; // to sort per clique size
    int *clqSizeCap = NULL; // current capacity in number of cliques for each thread

    int clqSizeCapIni = 4096;

    char *sense = NULL; // row senses
    ALLOCATE_VECTOR( sense, char, lp_rows(mip) );
    for ( int i=0 ; i<lp_rows(mip) ; ++i )
        sense[i] = toupper(lp_sense( mip, i ));

    char **clqNames;
    int clqnLines = 0;
    int clqnChars = 0;
    int clqnLinesCap = 8191;
    int clqnCharsCap = 64*clqnLinesCap;

    ALLOCATE_VECTOR( clqNames, char *, clqnLinesCap );
    ALLOCATE_VECTOR_INI( clqNames[0], char , clqnCharsCap );

    /* new rows */
    int *nrStart = NULL;
    int *nrIdx = NULL;
    double *nrCoef = NULL;
    double *nrRhs = NULL;
    char *nrSense = NULL;
    int nrRows = 0;
    int nrNz = 0;
    int nrCapRows = 0;
    int nrCapNz = 0;
    
    /* row names */
    char **nrNames;
    int nrnLines = 0;
    int nrnChars = 0;
    int nrnLinesCap = 8192;
    int nrnCharsCap = 64*nrnLinesCap;
    ALLOCATE_VECTOR( nrNames, char *, nrnLinesCap );
    ALLOCATE_VECTOR_INI( nrNames[0], char , nrnCharsCap );


    const int nCliques = check_cliques( mip, cliqueState, cliques, elRow, sizeRow, allEl, &nAllEl, sense );
    if ( nCliques == 0 )
        goto TERMINATE;

    /* filling cliques per col */
    {
        int totnz = 0;
        ALLOCATE_VECTOR_INI( nColClqs, int, lp_cols(mip) );
        for ( int i=0 ; (i<nCliques) ; ++i )
        {
            int ir = cliques[i];
            assert( sizeRow[ir] >= 2 );
            for ( int j=0 ; (j<sizeRow[ir]) ; ++j )
            {
                int col = elRow[ir][j];
                ++nColClqs[col];
                ++totnz;
            }
        }
        ALLOCATE_VECTOR( allCollClqs, int, totnz );
        ALLOCATE_VECTOR( colClqs, int *, lp_cols(mip) );
        // start of each column
        colClqs[0] = allCollClqs;
        for ( int i=1 ; i<lp_cols(mip) ; ++i )
            colClqs[i] = colClqs[i-1] + nColClqs[i-1];

        FILL( nColClqs, 0, lp_cols(mip), 0 );
        // filling cliques of each col
        for ( int i=0 ; (i<nCliques) ; ++i )
        {
            int ir = cliques[i];
            assert( sizeRow[ir] >= 2 );
            for ( int j=0 ; (j<sizeRow[ir]) ; ++j )
            {
                int col = elRow[ir][j];
                colClqs[col][nColClqs[col]++] = ir;
            }
        }
    }


    nrCapRows = nCliques*MIN(maxExtensions,5) + 1024;
    nrCapNz = lp_nz(mip) * 2;
    ALLOCATE_VECTOR( nrStart, int, nrCapRows+1 );
    ALLOCATE_VECTOR( nrIdx, int, nrCapNz );
    ALLOCATE_VECTOR( nrCoef, double, nrCapNz );
    ALLOCATE_VECTOR( nrRhs, double, nrCapRows );
    ALLOCATE_VECTOR( nrSense, char, nrCapRows );
    nrStart[0] = 0;

    // just to fill parameter
    ALLOCATE_VECTOR( rc, int, lp_cols(mip)*2 );
    FILL( rc, 0, lp_cols(mip)*2, 1.0 );


    /* allocating structures for threads */

    int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
{
    if (omp_get_thread_num()==0)
        nthreads = omp_get_num_threads();
}
#endif
    assert( nthreads>= 1);
    
    ALLOCATE_VECTOR( idx, int *, nthreads );
    ALLOCATE_VECTOR( coef, double *, nthreads );
    
    ALLOCATE_VECTOR( currClique, IntSet, nthreads );
    ALLOCATE_VECTOR_INI( clqsSize, struct CliqueSize *, nthreads );
    ALLOCATE_VECTOR_INI( clqSizeCap, int, nthreads );
    FILL( clqSizeCap, 0, nthreads, clqSizeCapIni );

    ALLOCATE_VECTOR( ivt, char * , nthreads );
    /* incidence vector for rows to check */
    ALLOCATE_VECTOR( ivrt, char * , nthreads );
    /* vector of rows to check */
    ALLOCATE_VECTOR( rtc, int * , nthreads );

    for ( int i=0 ; (i<nthreads) ; ++i )
    {
        vint_set_init( &currClique[i] );
        ALLOCATE_VECTOR( clqsSize[i], struct CliqueSize, clqSizeCap[i] );
        ALLOCATE_VECTOR_INI( ivt[i], char, lp_cols(mip)*2 );
        ALLOCATE_VECTOR_INI( ivrt[i], char, lp_rows(mip) );
        ALLOCATE_VECTOR_INI( rtc[i], int, lp_rows(mip) );
        ALLOCATE_VECTOR( idx[i], int, lp_cols(mip) );
        ALLOCATE_VECTOR( coef[i], double, lp_cols(mip) );
    }

    clock_t startExtend = clock();
#pragma omp parallel for
    for ( int iclq=0 ; iclq<nCliques ; ++iclq )
    {
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
        assert( thread < nthreads );
#endif
        int row = cliques[iclq];

        int *cidx = idx[thread];
        double *ccoef = coef[thread];
        if ( cliqueState[row] != Dominated )
        {
            int nz = lp_row( mip, row, cidx, ccoef );
            IntSet *clq = &currClique[thread];
            for ( int ii=0 ; (ii<nz) ; ++ii )
                if (ccoef[ii]<=-1e-5)
                    cidx[ii] += lp_cols(mip);
            vint_set_clear( clq );
            vint_set_add( clq, cidx, nz );

            /* extending */
            {
                CliqueExtender *clqe = clqe_create();
//                printf("aaa cgraph: %p\n", cgraph ); fflush(stdout); fflush(stderr);
//                printf("rc %p cgraph size %d\n", (void*) rc, cgraph_size(cgraph) ); fflush(stdout); fflush(stderr);
                clqe_set_max_it_bk( clqe, maxItBk );
                clqe_set_costs( clqe, rc, cgraph_size(cgraph) );

                clock_t startext = clock();
                int status = clqe_extend( clqe, cgraph, clq, lp_cols(mip), CLQEM_EXACT );
                clqMergeSecsExtend += ((double)clock()-startext)/(double)CLOCKS_PER_SEC;
                if (status>0)
                {
                    if (clqMergeVerbose>=2)
                    {
                        char rname[256] = "";
                        printf("-> constraint %s extended: \n", lp_row_name(mip, row, rname) );
                    }
                    ++clqMergeNExtended;

                    // to sort cliques found per size
                    struct CliqueSize *clqbs = clqsSize[thread];

                    const CliqueSet *clqs = clqe_get_cliques( clqe );
                    assert(clq_set_number_of_cliques(clqs)>=1);

                    int nCliquesToInsert = MIN( maxExtensions, clq_set_number_of_cliques(clqs) );
                    if (nCliquesToInsert)
                    {
                        cliqueState[row] = Dominated;
                        if (nCliquesToInsert==1)
                        {
                            clqbs[0].clqIdx = 0;
                            clqbs[0].size = clq_set_clique_size( clqs, 0 );
                        }
                        else
                        {
                            /* resize */
                            if (clq_set_number_of_cliques(clqs)>clqSizeCap[thread])
                            {
                                clqSizeCap[thread] = MAX( clqSizeCap[thread]*2, clq_set_number_of_cliques(clqs) );
                                clqbs = clqsSize[thread] = xrealloc( clqsSize[thread], sizeof(struct CliqueSize)*clqSizeCap[thread] );
                            }
                            /* sorting per size */
                            for ( int ic=0 ; (ic<clq_set_number_of_cliques(clqs)) ; ++ic )
                            {
                                clqbs[ic].clqIdx = ic;
                                clqbs[ic].size = clq_set_clique_size( clqs, ic );
                            }
                            /* maybe needed to sort cliques */

                            qsort( clqbs, clq_set_number_of_cliques(clqs), sizeof(struct CliqueSize), cmp_clq_size );
                            //printf( "from %d to %d\n", clqbs[0].size, clqbs[clq_set_number_of_cliques(clqs)-1].size );
                        }
                    }

                    for ( int ic=0 ; (ic<nCliquesToInsert) ; ++ic )
                    {
                        int idxclique = clqbs[ic].clqIdx;
                        int size = clq_set_clique_size( clqs, idxclique );
                        //printf("   adding %d\n", size );
                        const int *el = clq_set_clique_elements( clqs, idxclique );

                        char *iv = ivt[thread];
                        char *ivr = ivrt[thread];
                        int *rtct = rtc[thread];
                        add_clique( mip, newCliques, size, el, cliqueState, cliques,(const int **) elRow, sizeRow, iv, ivr, rtct, (const int **) colClqs, nColClqs );
#pragma omp critical
                        {
                            char origrname[256] = "";
                            lp_row_name( mip, row, origrname );
                            char extn[64] = "";
                            if (nCliquesToInsert>1)
                            {
                                sprintf(extn, "xt(%d)", ic );
                                strcat( origrname, extn );
                            }
                            if (clqMergeVerbose>=2)
                            {
                                printf("\t%s : ", origrname );
                            }
                            addName( &clqNames, &clqnLines, &clqnChars, &clqnLinesCap, &clqnCharsCap, origrname );
                        }
                        
                        if (clqMergeVerbose>=2)
                        {
                            for ( int k=0 ; (k<size) ; ++k )
                            {
                                char strneg[3] = "";
                                int icol = el[k];
                                if (icol>=lp_cols(mip))
                                {
                                    strcpy( strneg, "!");
                                    icol -= lp_cols(mip);
                                }
                                char coln[256];
                                lp_col_name( mip, icol, coln );
                                printf("%s%s ", strneg, coln );
                            }
                            printf("\n");
                        }
                    }
                }
                clqe_free( &clqe );
            }
        } // non dominated clique constraints
    } // all clique constraints
    
    clqMergeSecsExtendAndDominate = ((double)clock()-(double)startExtend)/((double)CLOCKS_PER_SEC);

    clock_t startRemAdd = clock();
    /* removing dominated cliques */
    {
        int nToRemove = 0;
        int *toRemove;
        ALLOCATE_VECTOR( toRemove, int, lp_rows(mip) );

        for ( int i=0 ; (i<lp_rows(mip)) ; ++i )
        {
            if ( cliqueState[i] == Dominated )
            {
                if (sense[i]=='E')
                {
                    ++clqMergeNDominatedEqPart;
                    // adding >= part, <= part will be added separately
                    int nz = lp_row( mip, i, idx[0], coef[0] );
                    double rhs = lp_rhs( mip, i );
                    char rname[256];
                    lp_row_name( mip, i, rname );
                    char nrname[512];
                    sprintf( nrname, "%sEp", rname );
                    addRow(  &nrRows, &nrCapRows, &nrNz, &nrCapNz, &nrStart, &nrIdx, &nrCoef, &nrSense, &nrRhs,
                        nz, idx[0], coef[0], 'G', rhs, nrname, &nrNames, &nrnLines, &nrnLinesCap, &nrnChars, &nrnCharsCap );
                }
                else
                {
                    ++clqMergeNDominatedFull;
                }
                toRemove[nToRemove++] = i;
            }
        }

        if (nToRemove)
            lp_remove_rows( mip, nToRemove, toRemove );

        free( toRemove );
    }

    /* adding stronger cliques */
    if (clq_set_number_of_cliques(newCliques)>0)
    {
        int nNewCliques = clq_set_number_of_cliques(newCliques);

        for ( int ic=0 ; (ic<nNewCliques) ; ++ic )
        {
            int size = clq_set_clique_size( newCliques, ic );
            const int *el = clq_set_clique_elements( newCliques, ic );
            double rhs = 1.0;
            for ( int i=0 ; (i<size) ; ++i )
            {
                if (el[i]>=lp_cols(mip))
                {
                    coef[0][i] = -1.0;
                    rhs -= 1.0;
                    idx[0][i] = el[i]-lp_cols(mip);
                }
                else
                {
                    idx[0][i] = el[i];
                    coef[0][i] = 1.0;
                }
            }
            // adding
            addRow(  &nrRows, &nrCapRows, &nrNz, &nrCapNz, &nrStart, &nrIdx, &nrCoef, &nrSense, &nrRhs,
                size, idx[0], coef[0], 'L', rhs,
                clqNames[ic], &nrNames, &nrnLines, &nrnLinesCap, &nrnChars, &nrnCharsCap );
        }
    }
    
    clqMergeSecsAddAndRemove = ((double)clock()-(double)startRemAdd)/((double)CLOCKS_PER_SEC);
    
    /* flushing rows */
    lp_add_rows( mip, nrRows, nrStart, nrIdx, nrCoef, nrSense, nrRhs, (const char**) nrNames );
    
    
    if (clqMergeVerbose)
    {
        printf("%d extended, %d dom full, %d dom eq.\n", clqMergeNExtended, clqMergeNDominatedFull, clqMergeNDominatedEqPart );
        fflush( stdout );
    }

TERMINATE:
    if (nrStart)
        free( nrStart );
    if (nrIdx)
        free( nrIdx );
    if (nrCoef)
        free( nrCoef );
    if (nrRhs)
        free( nrRhs );
    if (nrSense)
        free( nrSense );
    if (nrNames)
    {
        free( nrNames[0] );
        free( nrNames );
    }
 

    if (clqNames)
    {
        free( clqNames[0]);
        free( clqNames);
    }

    if (sense)
        free( sense );
    free( cliques );
    free( cliqueState );
    if (rc)
        free( rc );
    if (idx)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free( idx[i] );
        free( idx );
    }
    if (coef)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free( coef[i] );
        free( coef );
    }
    if (currClique)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            vint_set_clean( &currClique[i] );
        free( currClique );
    }
    if (clqsSize)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free( clqsSize[i] );
        free( clqsSize );
    }

    if (elRow)
        free( elRow );
    if (sizeRow)
        free( sizeRow );
    if (allEl)
        free(allEl);

    clq_set_free( &newCliques );

    if (ivt)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free(ivt[i]);
        free( ivt );
    }
    if (ivrt)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free(ivrt[i]);
        free( ivrt );
    }
    if (rtc)
    {
        for ( int i=0 ; (i<nthreads) ; ++i )
            free(rtc[i]);
        free( rtc );
    }

    if (clqSizeCap)
        free(clqSizeCap);

    if (nColClqs)
        free(nColClqs);
    if (allCollClqs )
        free(allCollClqs);
    if (colClqs)
        free( colClqs);
        
        
}

