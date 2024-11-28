/*

    File:           quit.c

    Created:        Mon May 24 14:44:35 EST 1993

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Print an error message and exit.

    Arguments:      Variable arguments, in the style of printf()

    Revisions:      none

    $Id: quit.c 3435 2012-10-10 00:06:55Z gor171 $

*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include "sjwlib.h"


/** Print an error message and exit the process.
  *
  * @param format printf-like format string.
  * @param ellipsis list of variable arguments appropriate to
  *        the format specified.
  */
void quit(char *format, ...)
{
    va_list args;

    va_start(args,format);
    fprintf(stderr,"\n\n");
    vfprintf(stderr,format,args);
    va_end(args);
    fprintf(stderr,"\n\n");
	exit(0);
}
