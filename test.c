#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf.h>

#include "recurse.h"


/* void */
/* gough_test (void) */
/* { */
/*   int tj; */
/*   int N = 1e9; */
/*   gsl_set_error_handler_off(); */

/*   for (tj = 0; tj < N; tj = (tj +1)*100-1) */
/*     { */
/*       gsl_sf_result rc, rw; */
/*       int ec = gsl_sf_coupling_3j_e (0, tj, tj, 0, tj, -tj, &rc); */
/*       int ew = gsl_sf_wigner_3j_e (0, tj, tj, 0, tj, -tj, &rw); */
/*       double exact = 1.0 / sqrt(tj + 1.0); */

/*       printf("j=%d/2 c3j(0 j j| 0 j -j)= % .18e %.5e (status=%d)\n", tj, rc.val, rc.err, ec); */
/*       printf("j=%d/2 w3j(0 j j| 0 j -j)= % .18e %.5e (status=%d)\n", tj, rw.val, rw.err, ew); */
/*       printf("j=%d/2 exact             = % .18e\n",tj, exact); */
/*       printf("j=%d/2 wigner-exact      = % .18e\n",tj, rw.val - exact); */
/*     } */
/* } */

void
check_3j_family_j_exact (const int two_jmax)
{
  int two_j;

  for (two_j = 0; two_j <=two_jmax; two_j = (two_j + 1)*100-1)
    {
      double exact = 1.0 / sqrt (two_j + 1.0);
      double gsl = gsl_sf_coupling_3j(two_j, two_j, 0, two_j, -two_j, 0);
      double *a;
      int tjmax, tjmin;

      wigner3j_family_j (two_j, two_j, two_j, -two_j, &a, &tjmax, &tjmin);

      printf ("two_j: %d\trecursive: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, a[0], exact, a[0]-exact);

      printf ("two_j: %d\t calc: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, gsl, exact, gsl-exact);

      
      free (a);
    }

  for (two_j = 0; two_j <=two_jmax; two_j = (two_j + 1)*100-1)
    {
      double exact = 1.0 / sqrt (two_j + 1.0);
      double gsl = gsl_sf_coupling_3j(two_j, two_j, 0, 0, 0, 0);
      double *a;
      int tjmax, tjmin;

      wigner3j_family_j (two_j, two_j, 0, 0, &a, &tjmax, &tjmin);

      printf ("two_j: %d\trecursive: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, a[0], exact, a[0]-exact);

      printf ("two_j: %d\t calc: % .18e\t exact: % .18e\tdiff: % .18e\n",
	      two_j, gsl, exact, gsl-exact);

      
      free (a);
    }

}

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


  check_3j_family_j_exact(1e9);

/*   free (a); */
/*   wigner3j_family_j (2, 2, -2, 2, &a, &two_jmin, &two_jmax); */
/*   free (a); */
/*   wigner3j_family_j (1, 2, -1, 0, &a, &two_jmin, &two_jmax); */


  return 0;
}

