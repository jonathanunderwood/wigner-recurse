#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_sf_coupling.h>

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
	    free (a);
	  }
}

void
check_3j_family_j (const int two_j1, const int two_j2, 
		   const int two_m1, const int two_m2)
{
  double *a;
  int two_jmin, two_jmax, imax, i;

  wigner3j_family_j (two_j1, two_j2, two_m1, two_m2, &a, &two_jmin, &two_jmax);

  imax = (two_jmax-two_jmin)/2;

  printf ("two_j1: %d two_j2: %d two_m1: %d two_m2: %d\n", 
	  two_j1, two_j2, two_m1, two_m2);

  for(i=0; i<=imax; i++)
    {
      int j=two_jmin+i;
      printf ("\ttwo_j3: %d\t %g\n", j, a[i]);
    }

  free (a);
}

void
check_3j_family_j_gsl (const int two_j1, const int two_j2, 
		       const int two_m1, const int two_m2)
{
  double *a;
  int two_jmin, two_jmax, imax, i;
  int two_m3 = -(two_m1 + two_m2);

  wigner3j_family_j (two_j1, two_j2, two_m1, two_m2, &a, &two_jmin, &two_jmax);

  imax = (two_jmax-two_jmin)/2;

  printf ("two_j1: %d two_j2: %d two_m1: %d two_m2: %d\n", 
	  two_j1, two_j2, two_m1, two_m2);
  printf ("--->%d\n", imax);
  for(i=0; i<=imax; i++)
    {
      int two_j3=two_jmin+i*2;
      
      double gsl=gsl_sf_coupling_3j(two_j1, two_j2, two_j3,
				    two_m1, two_m2, two_m3);
      printf ("\ttwo_j3: %d\t\twigner: %g\tgsl: %g\tdiff: %g\n", 
	      two_j3, a[i], gsl, a[i]-gsl);
    }

  free (a);
}

int
main ()
{
/*   double *a; */
/*   int i, two_jmin, two_jmax, imax; */

  //hammer_3j_j (2);

  check_3j_family_j_gsl (1, 1, -1, 1);
  check_3j_family_j_gsl (9, 1, -1, 1);
  check_3j_family_j_gsl (9, 8, -1, 0);
  check_3j_family_j_gsl (90, 80, -30, 10);


/*   free (a); */
/*   wigner3j_family_j (2, 2, -2, 2, &a, &two_jmin, &two_jmax); */
/*   free (a); */
/*   wigner3j_family_j (1, 2, -1, 0, &a, &two_jmin, &two_jmax); */


  return 0;
}

