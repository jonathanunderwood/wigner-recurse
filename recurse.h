/* This file is part of wigner-recurse, code for calculating Wigner 3j
 * and 6j symbols using recursion.
 *
 * Copyright (C) 2008-2016 Jonathan * G. Underwood.
 *
 * wigner-recurse is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * wigner-recurse is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with wigner-recurse.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef __WIGNER_RECURSE_H__
#define __WIGNER_RECURSE_H__

#ifndef INLINE
# if __GNUC__ && !__GNUC_STDC_INLINE__
#  define INLINE extern inline
# else
#  define INLINE inline
# endif
#endif

int wigner3j_family_j (const int two_j2, const int two_j3,
                       const int two_m2, const int two_m3,
                       double **family, int *two_jmin, int *two_jmax);

int wigner3j_family_m (const int two_j1, const int two_j2, const int two_j3,
                       const int two_m1, double **family, int *two_mmin,
                       int *two_mmax);

#define __MAX(a,b) ((a) > (b) ? (a) : (b))
#define __MIN(a,b) ((a) < (b) ? (a) : (b))

INLINE size_t
wigner3j_family_m_dim (const int two_j1, const int two_j2, const int two_j3,
                       const int two_m1)
/* For a given set of input parameters that would be used for a 3j
   family calculation using M recursion return the dimension of the
   array needed to hold the family of values. */
{
  int two_mmin = __MAX(-two_j2, -two_j3 - two_m1);
  int two_mmax = __MIN(two_j2, two_j3 - two_m1);

  return 1 + (two_mmax - two_mmin) / 2;
}

INLINE size_t
wigner3j_family_j_dim (const int two_j2, const int two_j3,
                       const int two_m2, const int two_m3)
/* For a given set of input parameters that would be used for a 3j
   family calculation using J recursion return the dimension of the
   array needed to hold the family of values. */
{
  int a = abs (two_j2 - two_j3);
  int b = abs (two_m2 + two_m3);

  int two_jmin = __MAX(a, b);
  int two_jmax = two_j2 + two_j3;

  return 1 + (two_jmax - two_jmin) / 2;
}

#undef __MAX
#undef __MIN

#endif /* __WIGNER_RECURSE_H__ */

