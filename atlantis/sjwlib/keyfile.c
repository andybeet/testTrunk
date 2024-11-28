/*

    File:           keyfile.c

    Created:        Tue Nov 8 13:28:29 EST 1994

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Routines to handle ascii keyword
		    files. These files are commonly
		    used by me to read key - parameter
		    pairs in various programs

    Arguments:      See each routine below.

    Revisions:      30/04/98 JRW
		    Permit parameters to be independent of case.
    $Id: keyfile.c 3417 2012-09-24 04:34:46Z gor171 $

*/


#include <ctype.h>
#include <string.h>
#ifndef WIN32
#include <unistd.h>
#include <strings.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sjwlib.h>


int keyprm_verbose = 0;
//void (*keyprm_errfn)(char *format, ...) = quit;
int keyprm_case_sensitive = 1;

#define STRCMP(a,b) ((keyprm_case_sensitive) ? strcmp(a,b) : strcasecmp(a,b))
#define STRNCMP(a,b,l) ((keyprm_case_sensitive) ? strncmp(a,b,l) : strncasecmp(a,b,l))

//keyprm_errfn
/** Set the behavior of the key finding routines
  * if a key is not found.
  *
  * @param fn pointer to error handling function.
  */
void set_keyprm_errfn(void (*fn) (char *, ...))
{
    keyprm_errfn = fn;
}

/** Set the case sensitivity for key finding routines.
  *
  * @param c 1 if case sensitive, else 0.
  */
void set_keyprm_case(int c)
{
    keyprm_case_sensitive = c;
}


/** Skip forward from the current file position to
  * the start of the next line beginning with key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @return non-zero if successful.
  */
int skipToKeyStart(FILE *fp, char *key)
{
    char buf[MAXLINELEN];
    size_t len = strlen(key);
    long fpos;
    char *s;
    int rewound = 0;

    do {
	fpos = ftell(fp);
	s = fgets(buf,MAXLINELEN,fp);
	if( s==NULL && !rewound ) {
	    fpos = 0L;
	    fseek(fp,fpos,0);
	    s = fgets(buf,MAXLINELEN,fp);
	    rewound = 1;
	}
    }
    while( s && STRNCMP(key,buf,len) != 0 );

    if( s == NULL ) {
	(*keyprm_errfn)("skipToKeyStart: key %s not found\n",key);
	return(0);
    }

    /* seek to start of line */
    if( fseek(fp,fpos,0) < 0 )
	return(0);

    return(1);
}

/** Skip forward from the current file position to
  * the next line beginning with key, positioned at the character
  * immediately after the key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @return non-zero if successful.
  */
int skipToKeyEnd(FILE *fp, char *key)
{
    int i;
    size_t len = strlen(key);
    long fpos;
    char *s;
    int rewound = 0;
    char buf[MAXLINELEN];

    do {
	fpos = ftell(fp);
	s = fgets(buf,MAXLINELEN,fp);

	if( s==NULL && !rewound ) {
	    fpos = 0L;
	    fseek(fp,fpos,0);
	    s = fgets(buf,MAXLINELEN,fp);
	    rewound = 1;
	}

/* Truncate the string at the first space after the key length.
 */
	for (i=(int)len; s && i<(int)strlen(buf); ++i)
	   if (isspace(buf[i])) {
	      buf[i]='\000';
	      break;
	   }
    }
    while( s && STRCMP(key,buf) != 0 );

    if( s == NULL ) {
	(*keyprm_errfn)("skipToKeyEnd: key %s not found\n",key);
	return(0);
    }

    /* seek to character after key */
    if( fseek(fp,fpos+(int)len,0) < 0 )
	return(0);

    return(1);
}


/** Read an integer value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned integer value.
  * @return non-zero if successful.
  */
int readkeyprm_i(FILE *fp, char *key, int *p)
{
    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fscanf(fp,"%d",p) != 1 )
	quit("readkeyprm_i: Can't read %s\n",key);

    /* Echo to stdout if verbose enough */
    if( keyprm_verbose > 1 )
	printf("%s  %d\n",key,*p);

    return(1);
}

/** Read an double value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned double value.
  * @return non-zero if successful.
  */
int
readkeyprm_d(FILE *fp, char *key, double *p)
{
    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
		return(0);

    if( fscanf(fp,"%lf",p) != 1 )
	quit("readkeyprm_d: Can't read %s\n",key);

    /* Echo to stdout if verbose enough */
    if( keyprm_verbose > 1 )
	printf("%s  %.10g\n",key,*p);

    return(1);
}


/** Read an array of double value following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned array of double value.
  * @param size pointer to returned array size.
  * @return non-zero if successful.
  */
int
readkeyprm_darray(FILE *fp, char *key, double **p, int *size)
{
int i=0;

    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fscanf(fp,"%d",size) != 1  || *size < 0) {
	keyprm_errfn("readkeyprm_darray: Parameter %s array size format bad\n",
			key);
	return -1;

    }

    if (*size > 0) {
       double val = 0.0;
       double lastval = 0.0;

/* allocate memory for the array if required */

       if( *p==NULL && (*p=(double *)malloc((size_t)(*size)*sizeof(double))) == NULL ) {
           keyprm_errfn("readkeyprm_darray: %s: Can't read array value\n",key);
           return(-1);
       }

/**** Now populate the array, if we run out of values, then pad
 **** to end with the last value
 ****/
       for (i=0; i<(*size); ++i) {
	   if(fscanf(fp,"%lf",&val) != 1 ) {
	      if ( i == 0 ) {
	         keyprm_errfn("readkeyprm_darray: %s: Can't read any array values\n", key);
		 return -1;
	      }
	      else {
/* not enough values - fill remainder with last value */
                 warn("readkeyprm_darray: %s: Filling to end of array with last value, %f\n",key, lastval);
                 for(; i<(*size); ++i)
                    (*p)[i] = lastval;
		 break;
	      }
	   }
	   else {
/* A single value was located, so populate the array */
		   if((*p) == NULL)
			   quit("readkeyprm_darray array is null. Should have been a double array of length %d\n", *size);

	      (*p)[i] = val;
	      lastval = val;
	   }
       }
    }


       /* Echo to stdout if verbose enough */
    if( keyprm_verbose > 1 ) {
        printf("%s  %d\n",key, *size);
	for (i=0; i<(*size); ++i)
	   printf("%f\n", (*p)[i]);
    }

    return 1;
}



/** Read an array of int values following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned array of int values.
  * @param size pointer to returned array size.
  * @return non-zero if successful.
  */
int
readkeyprm_iarray(FILE *fp, char *key, int **p, int *size)
{
int i=0;

    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fscanf(fp,"%d",size) != 1  || *size < 0) {
	keyprm_errfn("readkeyprm_iarray: Parameter %s array size format bad\n",
			key);
	return -1;
    }

    if (*size > 0) {
       int val = 0;
       int lastval = 0;

/* allocate memory for the array if required */

       if( *p==NULL && (*p=(int *)malloc((size_t)(*size)*sizeof(int))) == NULL ) {
           keyprm_errfn("readkeyprm_iarray: %s: Can't read array value\n",key);
           return(-1);
       }

/**** Now populate the array, if we run out of values, then pad
 **** to end with the last value
 ****/
       for (i=0; i<(*size); ++i) {
	   if(fscanf(fp,"%d",&val) != 1 ) {
	      if ( i == 0 ) {
	         keyprm_errfn("readkeyprm_iarray: %s: Can't read any array values\n", key);
		 return -1;
	      }
	      else {
/* not enough values - fill remainder with last value */
                 warn("readkeyprm_iarray: %s: Filling to end of array with last value, %d\n",key, lastval);
                 for(; i<(*size); ++i)
                    (*p)[i] = lastval;
		 break;
	      }
	   }
	   else {
/* A single value was located, so populate the array */
	      (*p)[i] = val;
	      lastval = val;
	   }
       }
    }


       /* Echo to stdout if verbose enough */
    if( keyprm_verbose > 1 ) {
        printf("%s  %d\n",key, *size);
	for (i=0; i<(*size); ++i)
	   printf("%d\n", (*p)[i]);
    }

    return 1;
}



/** Read an array of string values following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned array of int values.
  * @param size pointer to returned array size.
  * @return non-zero if successful.
  */
int
readkeyprm_sarray(FILE *fp, char *key, char ***p, int *size)
{
	int i=0;

    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fscanf(fp,"%d",size) != 1  || *size < 0) {
    	keyprm_errfn("readkeyprm_sarray: Parameter %s array size format bad\n", key);
    	return -1;
    }

    if (*size > 0) {
      char val[1024];

/* allocate memory for the array if required */

       if((*p=(char **)malloc((size_t)(*size)*sizeof(char *))) == NULL ) {
           keyprm_errfn("readkeyprm_sarray: %s: Can't read array value\n",key);
           return(-1);
       }

       /**** Now populate the array, if we run out of values, then pad
        **** to end with the last value
        ****/
       for (i=0; i<(*size); ++i) {
			if(fscanf(fp,"%s",val) != 1 ) {
				if ( i == 0 ) {
					keyprm_errfn("readkeyprm_sarray: %s: Can't read any array values\n", key);
					return -1;
				}
				else {
					/* not enough values - fill remainder with last value */
					quit("readkeyprm_sarray: %s: Not enough values in array value %s\n", key);
				}
			}
			else {
				/* A single value was located, so populate the array */
				(*p)[i] = (char *)malloc(sizeof(char) * strlen(val) + 1);
				strcpy((*p)[i], val);
			}
		}
	}

    return 1;
}


/** Read a character string following a key.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param p pointer to returned character string.
  * @return non-zero if successful.
  */
int readkeyprm_s(FILE *fp, char *key, char *p)
{
    char line[MAXLINELEN];
    char *s;
    char *r;

    /* Read parameter */
    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fgets(line,MAXLINELEN,fp) != line )
	quit("readkeyprm_s: Can't read %s\n",key);

    /* Strip leading space */
    for(s=line; *s && isspace(*s); s++) /* loop */;

    /* Copy out result */
    for( r=p; *s && (*s != '\n'); *r++ = *s++) /* loop */ ;
    *r = 0;

    /* Echo to stdout if verbose enough */
    if( keyprm_verbose > 1 )
	printf("%s  %s\n",key,p);

    return(1);
}

/** Read a time value and convert to seconds.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param val pointer to returned time value.
  * @return non-zero if successful.
  */
int time_sec_param(FILE *fp, char *key, double *val)
{
    char buf[MAXLINELEN];

    if( !skipToKeyEnd(fp,key) )
	return(0);

    if( fgets(buf,MAXLINELEN,fp) == NULL || !strtosecs(buf,val) )
	quit("timeparam: Can't read %s\n",key);

    return(1);
}


#if 0
/* scanf-like routine which matches a key at the start of
  * the line * before doing the read.
  *
  * @param fp pointer to stdio FILE structure.
  * @param key keyname to locate in file.
  * @param format scanf format character string.
  * @param ellipsis remaining arguments as per scanf.
  * @return non-zero if successful.
  *
  * @see scanf.
  */
#include <stdarg.h>

int
keyfscanf(FILE *fp, char *key, char *format, ...)
{
    va_list args;
    int n;

    va_start(args,format);
    skipToKeyEnd(fp,key);
    n = vfscanf(fp,format,args);
    va_end(args);
    return(n);
}
#endif
