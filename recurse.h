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

#ifndef __RECURE_H__
#define __RECURE_H__

void wigner3j_family_j (const int two_j2, const int two_j3, 
			const int two_m2, const int two_m3,
			double **family, int *two_jmin, int *two_jmax);

#endif /* __RECURE_H__ */

