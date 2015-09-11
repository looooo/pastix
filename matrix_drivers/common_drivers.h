/*
 *  File: common_drivers.h
 *
 *  Definition common to all drivers.
 */
#ifndef COMMON_DRIVER_H
#define COMMON_DRIVER_H

#include <assert.h>

#define STR_SIZE 256
#define FGETS(line, BUFSIZ, infile)                                     \
    {                                                                   \
        if (NULL == fgets(line, BUFSIZ, infile))                        \
        {                                                               \
            fprintf(stderr, "ERROR: %s:%d fgets\n", __FILE__, __LINE__); \
            exit(1);                                                    \
        }                                                               \
    }

#define MALLOC_ERROR(_str_)                                             \
    {                                                                   \
        fprintf(stderr, __FILE__ ":%d: " _str_ "\n", __LINE__ );        \
        exit(-1);                                                       \
    }

#define memFree_null( _ptr_ ) { free(_ptr_); _ptr_ = NULL; }

/*
  Function: myupcase

  Rewrites *s* to upper case.

  Parameters:
    s - string to rexwrite in upcase.
*/
void myupcase(char *S);

/*
  Function:  mysubstr

  Copy len element, from *S[pos]* in *s*.

  Parameters:
    s   - destination
    S   - Source
    pos - sarting position
    len - Size of the string to copy.
*/
void mysubstr(char *s, const char *S, const pastix_int_t pos, const pastix_int_t len);

/*
  Function:  mysubstr2

  Copy the number placed between a and b in fmt.

  Parameters:
    fmt - String in which there is a and b
    a   - first element
    b   - last element
    val - the integer between a and b

*/
void mysubstr2(const char *fmt, const char a, const char b, pastix_int_t *val);

#endif /* not COMMON_DRIVER_H */
