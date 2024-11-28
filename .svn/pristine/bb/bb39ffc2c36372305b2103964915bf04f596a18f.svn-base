/*

    File:           text_input.c

    Created:        Thu Dec 3 11:01:45 EST 1992

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Useful routines to deal with text files.

                    Read the next non-blank non-comment line
		    from an ascii file. Uses fgets() to read
		    the line, so the terminating newline
		    character is read.
		    Comment lines are those lines which have
		    '#' as their first character.

    Revisions:      16/12/1994 SJW
		    Combined nextline.c, flushline.c,
		    added parseline().

    $Id: text_input.c 3494 2012-11-08 09:41:17Z gor171 $

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

int iscommentline(char *line);
int isblankline(char *line);

/** Read the next non-blank non-comment line
  * from an ascii file. Uses fgets() to read
  * the line, so the terminating newline
  * character is read.
  * Comment lines are those lines which have
  * '#' as their first character.
  *
  * @param line pointer to storage for line read
  * @param n max length of line to read
  * @param fp stream from which to read
  *
  * @return non-zero if line was read successfilly.
  */
int nextline(char *line, long int n, FILE *fp)
{
    char *s;

    do {
	s = fgets(line, (int)n, fp);
    } while( s && (iscommentline(line)||isblankline(line)) );

    if( s == NULL )
	return(0);

    return(1);
}

int
iscommentline(char *line)
{
    return( line[0] == '#' );
}

int
isblankline(char *line)
{
    while( *line && isspace(*line) ) line++;
    if( *line )
        return(0);
    return(1);
}

/** Parse a line and split it into fields as delimited
  * by whitespace.
  *
  * @param line pointer to storage for line read
  * @param str returned array of string fields.
  * @param max mamimum number of fields.
  *
  * @return number of fields read.
  */
int parseline(char *line, char **str, int max)
{
    int quote;
    int n = -1;

    /* Skip whitespace */
    while((*line && isspace(*line)) || (*line && *line == '+'))
    	line++;

    /* Loop to find string separators */
    for(n=0; *line && n < max; n++) {
        /* Check for start of quote */
    	if( *line == '\'' || *line == '"' )
    	    quote = *line++;
    	else
    	    quote = 0;
    	/* Store string start */
    	str[n] = line;
    	/* Loop to end of string */
    	while( *line && ((quote&&(*line!=quote)) || (!quote && !isspace(*line) )) )
    		line++;
    	if (*line) *line++ = 0;
    	while((*line && isspace(*line)) || (*line && *line == '+'))
    		line++;
    }
    return( *line ? n+1 : n );
}
/** Read stream to end of line.
  *
  * @param fp stream from which to read
  */
void
flushline(FILE *fp)
{
	int c;

	while( (c=fgetc(fp)) != EOF && c != '\n' )
		/* loop */ ;
}

