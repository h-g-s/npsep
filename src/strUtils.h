#include <ctype.h>

/**
 * removes spaces at the beginning
 * of the string
 **/
char *removeLeadingSpaces( char *dest, const char *str );

/**
 * number of numeric characters in
 * line
 **/
int digitsInLine( const char *str, int STR_SIZE );

/**
 * receives a fileName with path
 * and returns only the fileName
 **/
char *getFileName(char *destiny, const char *fileWithPath);

/**
 * splits a string usign a delimiter char
 * returns how many columns were found
 * multDel = 1 indicates that multiple delimiters
 * grouped together represent only one
 **/
int splitString( char **columns, const char *str, const char delimiter,
                 const int maxColumns, const int columnSize, const char multDel );

char* getParamName(char* target, const char* str);

char* getParamValue(char* target, const char* str);

/* removes spaces and other meaningless characthers from end of line */
void strRemoveEndOfLine( char *str );

