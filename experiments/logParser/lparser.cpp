#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

int main( int argc, const char **argv )
{
   FILE *f = fopen( argv[1], "r" );
   if (!f)
   {
      fprintf( stderr, "Could not open file %s ", &(argv[1]) );
      exit( EXIT_FAILURE );
   }

   char line[256];

   double lastLB = 0;

   while ( fgets( line, 256, f ) )
   {
      if ( strstr( line, "round" ) )
      {
         char *s = strstr( line, "ated " );
         assert(s);
         s += 5;

         char scuts[256];
         strcpy( scuts, s );
         (*strstr( scuts, " ")) = '\0';

         int nCuts = atoi( scuts );

         char stime[256];
         strcpy( stime, strstr(s,"time: " )+6 );
         (*strstr( stime, "\n")) = '\0';
         double nSecs = atof( stime );

         char slb[256];
         strcpy( slb, strstr(s,"now: " )+5 );
         (*strstr( slb, " ")) = '\0';
         double lb = atof( slb );

         double impr = 0;

         static double sumImpr = 0.0;

         if ( lastLB > 1e-5 )
         {
            impr = ( ( lb - lastLB ) / lastLB ) * 100.0;
            sumImpr += impr;
         }

         static int iteration = 0;

         printf( "%d %d %g %g %g\n", ++iteration, nSecs, nCuts, lb, sumImpr );

         lastLB = lb;
      }
   }

   fclose(f);

   exit( EXIT_SUCCESS );
}
