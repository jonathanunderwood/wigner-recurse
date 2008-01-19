#include <stdlib.h>
#include <stdio.h>
#include "recurse.h"


void
hammer_3j_j (const int two_jmax)
{
  double *a;
  int two_j2, two_j3, two_m2, two_m3, tjmax, tjmin;
  
  for (two_j2 = 0; two_j2 <= two_jmax; two_j2++)
    for (two_m2 = -two_j2; two_m2 <=two_j2; two_m2 += 2)
      for (two_j3 = 0; two_j3 <= two_jmax; two_j3++)
	for (two_m3 = -two_j3; two_m3 <=two_j3; two_m3+=2)
	  {
	    printf ("%d %d %d %d\n", two_j2, two_j3, two_m2, two_m3);
	    wigner3j_family_j (two_j2, two_j3, two_m2, two_m3, &a, &tjmax, &tjmin);
	  }
}

int
main ()
{
  double *a;
  int i, two_jmin, two_jmax, imax;

  //hammer_3j_j (10);

  wigner3j_family_j (2, 2, 0, 0, &a, &two_jmin, &two_jmax);

  imax = (two_jmax-two_jmin)/2;
  printf ("2jmin:%d 2jmax:%d imax:%d\n", two_jmin, two_jmax, imax);

  for(i=0; i<=imax; i++)
    printf ("%g\n", a[i]);

  free (a);

  return 0;
}

