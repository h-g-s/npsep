#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int splitString( char **columns, const char *str, const char delimiter,
      const int maxColumns, const int columnSize, const char multDel )
{
   int sizeColumn;
   int ncolumn = 0;
   const char *send = str + strlen(str);
   const char *s = str;
   if (str[0] == '\0')
      return 0;
   const char *ns = s;
PROCESS_COLUMN:
   if ( ncolumn+1 == maxColumns )
      return ncolumn;

   /* finds the next delimiter */
FIND_DELIMITER:
   if ( ns == send )
      goto FOUND_COLUMN;
   if ( *ns != delimiter )
   {
      ns++;
      goto FIND_DELIMITER;
   }
FOUND_COLUMN:
   sizeColumn = ns - s;
   if ((!multDel)||(sizeColumn>0))
   {
      if (sizeColumn)
         memcpy( columns[ncolumn], s, sizeColumn );
      columns[ncolumn][sizeColumn] = '\0';
      ncolumn++;
   }
   if ( ns == send )
      return ncolumn;
   ++ns;
   s = ns;
   if ( ns != send )
      goto PROCESS_COLUMN;

   return ncolumn;
}

char* applyInversion(char* str) {
   int head, tail;
   char t;
   for ( head=0,tail=(strlen(str)-1) ; head<tail ; head++,tail-- ) {
      t = str[head];
      str[head] = str[tail];
      str[tail] = t;
   }
   return str;
}

char *getFileName(char *destiny, const char *fileWithPath) {
   int i,pos;
   /* Returning till found a slash */
   for ( i=(strlen(fileWithPath)-1),pos=0 ; (i>=0) ; i--,pos++ ) {
      if (fileWithPath[i]=='/') break;
      destiny[pos] = fileWithPath[i];
   }
   destiny[pos]='\0';
   applyInversion(destiny);
   /* Now, removing the . */
   const int endIdx = strlen(destiny)-1;
   for ( i=endIdx ; (i>=0) ; i-- )
   {
      if (destiny[i]=='.') 
         break;
   }
   if (destiny[i]=='.') 
      destiny[i]='\0';

   return destiny;
}

char *removeLeadingSpaces( char *dest, const char *str )
{
   char *startDest = dest;
   const char *send = str + strlen( str );
   const char *s = str;
   while ( (*s==' ') && (s<send) )
      s++;
   while ( (s<send) )
   {
      *dest = *s;
      ++dest;
      ++s;
   }
   *dest = '\0';

   return startDest;
}

int digitsInLine( const char *str, int STR_SIZE )
{
   int result = 0;

   {
      int i;
      for ( i=0 ; (i<STR_SIZE) ; ++i )
      {
         switch (str[i])
         {
            case '\0' :
            case '\n' :
               return result;
               break;
            default:
               if ( isdigit( str[i] ) )
                  result++;
               break;
         };
      }
   }

   return result;
}

char* getParamName(char* target, const char* str) {
   unsigned int i;
   unsigned int size;
   unsigned int pos;
   pos = 0;
   size = strlen(str);
   for ( i=0 ; i<size ; i++ ) {
      if (str[i] != '-') {
         if (str[i] != '=') {
            target[pos++] = str[i];
         } else {
            break;
         }
      }
   }
   target[pos] = '\0';
   return target;
}

char* getParamValue(char* target, const char* str) {
   unsigned int i;
   unsigned int size;
   size = strlen(str);
   unsigned int destSize = 0;
   for ( i=0 ; i<size ; i++) {
      if (str[i]=='=') break;
   }
   i++;
   for ( ; i<size ; i++) {
      target[destSize] = str[i];
      destSize++;
   }
   target[destSize] = 0;
   return target;
}

void strRemoveEndOfLine( char *str )
{
   char *p;

   // removing spaces at the end */
   p = str + strlen(str)-1;
   while ( p>=str )
   {
      if ((*p==' ') || (*p=='\r') || (*p=='\n'))
         *p = '\0';
      else
         break;

      --p;
   }
}

