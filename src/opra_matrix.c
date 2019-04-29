#include "scpc.h"
#include "opra_matrix.h"

// with individuals as cols and sites as rows
// s stands for the number of sites
// n stands for the number of individuals
void GlobalDeletionDNA(unsigned char *x, int *n, int *s)
{
    int i, j, k, cs = *s;

    for (j = 0; j < *s; j++)
    {
        k = i = *n * j;
	    while (i < *n * (j + 1))
        {
            if (KnownBase(x[i]))
                i++;
            else{
				memcpy(x+k, x+k+(*n), (*n)*(cs-(j+1)));
				(*s)--;j--;
                break;
	        }
	    }
    }
}

// \param[in] x input matrix
// \param[in] n number of rows of m
// \param[in] s number of columns of m
//
// Performs in-place transposition of matrix m.
void transpose(unsigned char *x, const int n, const int s)
{
  for (unsigned start = 0; start <= s * n - 1; ++start)
  {
    unsigned next = start;
    unsigned i = 0;
    do
    {
      ++i;
      next = (next % n) * s + next / n;
    } while (next > start);

    if (next >= start && i != 1)
    {
      const unsigned char tmp = x[start];
      next = start;
      do
      {
        i = (next % n) * s + next / n;
        x[next] = (i == start) ? tmp : x[i];
        next = i;
      } while (next > start);
    }
  }
}