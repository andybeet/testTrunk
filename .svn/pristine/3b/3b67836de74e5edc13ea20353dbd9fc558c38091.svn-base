/*

    File:           win32.c

    Created:        1997

    Author:         Adam Davidson
                    CSIRO Division of Oceanography

    Purpose:        This file contains routines standard in Unix
					C, but not in Win32 versions.

    Revisions:      None

*/

#include <ctype.h>
#include <string.h>


#define min(x,y) ((x)<(y) ? (x) : (y))



int strcasecmp(const char *s1, const char *s2) {

	int i,n1,n2;

	n1 = (int)strlen(s1);
	n2 = (int)strlen(s2);

	if (n1 != n2) return(-1);

	for (i=0;(s1[i] != '\0');i++) {
		if (tolower(s1[i]) == tolower(s2[i])) {
		}
		else {
			return(-1);
		}
	}

	return(0);
}

int strncasecmp(const char *s1, const char *s2, size_t n) {

	unsigned int i;


	if (n > (int) min(strlen(s1),strlen(s2)) ) return(-1);

	for (i=0;i<=(n-1);i++) {
		if (tolower(s1[i]) == tolower(s2[i])) {
		}
		else {
			return(-1);
		}
	}
	return(0);
}

