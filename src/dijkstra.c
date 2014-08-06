#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <spaths.h>
#include <limits.h>

void savePaths( ShortestPathsFinder *spf, int source );

int main( int argc, char **argv )
{
   if ( argc < 4 )
   {
      printf("uso: \n");
      printf("\tdijkstra graphFile executions debug\n");
      exit( EXIT_FAILURE );
   }

   int nrExecucoes = atoi( argv[2] );
   int depurar = atoi( argv[3] );

   srand(0);

   ShortestPathsFinder *spf = spf_load_gr( argv[1] );


   printf("dijkstra implementation by Haroldo, 2010.\n" );
   printf("loaded graph with %d nodes and %d arcs.\n", spf_nodes(spf), spf_arcs(spf) );

   clock_t start = clock();

   if (depurar)
   {
      printf("debug mode.\n");
      int source = 0;
      printf("source: %d\n", source );
      spf_find( spf, source );

      savePaths( spf, source );
      printf("file spaths.txt saved.\n");
      goto TERMINATE;
   }
   else
   {
      printf("running %d executions ...", nrExecucoes ); fflush( stdout );
      int exec;
      for ( exec=0 ; (exec<nrExecucoes) ; ++exec )
      {
         int source = rand()%spf_nodes( spf );
         spf_find( spf, source );
      }
      clock_t end = clock();
      printf(" done in %.5f seconds.\n", ((double)end-start)/((double)CLOCKS_PER_SEC) );
   }

TERMINATE:

   spf_free( &spf );
}

void savePaths( ShortestPathsFinder *spf, int source )
{
   int path[ spf_nodes( spf )+1 ]; 
   FILE *f = fopen("spaths.txt", "w");
   if (!f)
   {
      fprintf( stderr, "Could not open file spaths.txt\n" );
      exit( EXIT_FAILURE );
   }

   if ( spf_fw_ran( spf ) )
   {
      int i;
      for ( i=0 ; (i<spf_nodes(spf)) ; ++i )
      {
         if (i==source)
            continue;
         int size = spf_get_path_fw( spf, source, i, path );
         int dist = spf_fw_get_dist( spf, source, i );
         if ( dist == SP_INFTY_DIST )
         {
            fprintf(f, "[%d,%d](infty)\n", source+1, i+1 );
         }
         else
         {
            fprintf(f, "[%d,%d](%d)", source+1, i+1, dist );
            {
               int j;
               for ( j=0 ; (j<size) ; ++j )
                  fprintf(f, " %d", path[j]+1 );
            }
            fprintf(f, "\n");
         }
      }
   }
   else
   {
      int i;
      for ( i=0 ; (i<spf_nodes(spf)) ; ++i )
      {
         if (i==source)
            continue;
         int size = spf_get_path( spf, i, path );
         int dist = spf_get_dist( spf, i );
         if ( dist == SP_INFTY_DIST )
         {
            fprintf(f, "[%d,%d](infty)\n", source+1, i+1 );
         }
         else
         {
            fprintf(f, "[%d,%d](%d)", source+1, i+1, dist );
            int i;
            for ( i=0 ; (i<size) ; ++i )
               fprintf(f, " %d", path[i]+1 );
            fprintf(f, "\n");
         }
      }
   }

   fclose(f);
}
