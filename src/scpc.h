#ifndef __SCPC_H__
#define __SCPC_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

#ifdef WIN32 // Windows
  #define PATH_SPLIT '\'
#else // Unix
  #define PATH_SPLIT '/'
#endif

/* from R: print(log(4), d = 22) */
#define LN4 1.386294361119890572454

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* returns 1 if the base is adenine surely, 0 otherwise */
#define IsAdenine(a) a == 136

/* returns 1 if the base is guanine surely, 0 otherwise */
#define IsGuanine(a) a == 72

/* returns 1 if the base is cytosine surely, 0 otherwise */
#define IsCytosine(a) a == 40

/* returns 1 if the base is thymine surely, 0 otherwise */
#define IsThymine(a) a == 24

/* returns 1 if the base is a purine surely, 0 otherwise */
#define IsPurine(a) a > 63

/* returns 1 if the base is a pyrimidine surely, 0 otherwise */
#define IsPyrimidine(a) a < 64

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) (a & b) < 16

/* returns 1 if both bases are the same surely, 0 otherwise */
#define SameBase(a, b) KnownBase(a) && a == b

#endif // endif __COMMON_H