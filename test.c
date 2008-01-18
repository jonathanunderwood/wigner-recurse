#include <stdlib.h>
#include <stdio.h>
#include "recurse.h"


int
main ()
{
  double *a;
  int i, two_jmin, two_jmax, imax;

  wigner3j_family_j (6000, 2000, 200, -350, &a, &two_jmin, &two_jmax);

  imax = (two_jmax-two_jmin)/2;
  printf ("2jmin:%d 2jmax:%d imax:%d\n", two_jmin, two_jmax, imax);

  for(i=0; i<=imax; i++)
    printf ("%g\n", a[i]);

  free (a);

  return 0;
}
