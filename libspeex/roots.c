/* Copyright (C) 1981-1999 Ken Turkowski. <turk@computer.org>
 *
 * All rights reserved.
 *
 * Warranty Information
 *  Even though I have reviewed this software, I make no warranty
 *  or representation, either express or implied, with respect to this
 *  software, its quality, accuracy, merchantability, or fitness for a
 *  particular purpose.  As a result, this software is provided "as is,"
 *  and you, its user, are assuming the entire risk as to its quality
 *  and accuracy.
 *
 * This code may be used and freely distributed as long as it includes
 * this copyright notice and the above warranty information.


   Code slightly modified by Jean-Marc Valin (2002)

   Speex License:

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
   
   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 */

#include <math.h>

/*******************************************************************************
 * FindPolynomialRoots
 *
 * The Bairstow and Newton correction formulae are used for a simultaneous
 * linear and quadratic iterated synthetic division.  The coefficients of
 * a polynomial of degree n are given as a[i] (i=0,i,..., n) where a[0] is
 * the constant term.  The coefficients are scaled by dividing them by
 * their geometric mean.  The Bairstow or Newton iteration method will
 * nearly always converge to the number of figures carried, fig, either to
 * root values or to their reciprocals.  If the simultaneous Newton and
 * Bairstow iteration fails to converge on root values or their
 * reciprocals in maxiter iterations, the convergence requirement will be
 * successively reduced by one decimal figure.  This program anticipates
 * and protects against loss of significance in the quadratic synthetic
 * division.  (Refer to "On Programming the Numerical Solution of
 * Polynomial Equations," by K. W. Ellenberger, Commun. ACM 3 (Dec. 1960),
 * 644-647.)  The real and imaginary part of each root is stated as u[i]
 * and v[i], respectively, together with the corresponding constant,
 * conv[i], used in the convergence test.  This program has been used
 * successfully for over a year on the Bendix G15-D (Intercard System) and
 * has recently been coded for the IBM 709 (Fortran system).
 *
 * ACM algorithm #30 - Numerical Solution of the Polynomial Equation
 * K. W. Ellenberger
 * Missle Division, North American Aviation, Downey, California
 * Converted to C, modified, optimized, and structured by
 * Ken Turkowski
 * CADLINC, Inc., Palo Alto, California
 *******************************************************************************/

#define MAXN 20

void
poly_roots(
	const float		*a,			/* Coefficients */
	float			*u,			/* Real component of each root */
	float			*v,			/* Imaginary component of each root */
	float			*conv,		/* Convergence constant associated with each root */
	register long	n,			/* Degree of polynomial (order-1) */
	long			maxiter,	/* Maximum number of iterations */
	long			fig			/* The number of decimal figures to be computed */
)
{
	int i;
	register int j;
	float h[MAXN + 3], b[MAXN + 3], c[MAXN + 3], d[MAXN + 3], e[MAXN + 3];
	/* [-2 : n] */
	float K, ps, qs, pt, qt, s, rev, r=0;
	int t;
	float p=0, q=0;

	/* Zero elements with negative indices */
	b[2 + -1] = b[2 + -2] =
	c[2 + -1] = c[2 + -2] =
	d[2 + -1] = d[2 + -2] =
	e[2 + -1] = e[2 + -2] =
	h[2 + -1] = h[2 + -2] = 0.0;

	/* Copy polynomial coefficients to working storage */
	for (j = 0; j <= n; j++)
		h[2 + j] = *a++;						/* Note reversal of coefficients */

	t = 1;
	K = pow(10.0, (double)(fig));				/* Relative accuracy */

	for (; h[2 + n] == 0.0; n--) {				/* Look for zero high-order coeff. */
		*u++ = 0.0;
		*v++ = 0.0;
		*conv++ = K;
	}

INIT:
	if (n == 0)
		return;

	ps = qs = pt = qt = s = 0.0;
	rev = 1.0;
	K = pow(10.0, (double)(fig));

	if (n == 1) {
		r = -h[2 + 1] / h[2 + 0];
		goto LINEAR;
	}

	for (j = n; j >= 0; j--)					/* Find geometric mean of coeff's */
		if (h[2 + j] != 0.0)
			s += log(fabs(h[2 + j]));
	s = exp(s / (n + 1));

	for (j = n; j >= 0; j--)					/* Normalize coeff's by mean */
		h[2 + j] /= s;

	if (fabs(h[2 + 1] / h[2 + 0]) < fabs(h[2 + n - 1] / h[2 + n])) {
REVERSE:
		t = -t;
		for (j = (n - 1) / 2; j >= 0; j--) {
			s = h[2 + j];
			h[2 + j] = h[2 + n - j];
			h[2 + n - j] = s;
		}
	}
	if (qs != 0.0) {
		p = ps;
		q = qs;
	} else {
		if (h[2 + n - 2] == 0.0) {
			q = 1.0;
			p = -2.0;
		} else {
			q = h[2 + n] / h[2 + n - 2];
			p = (h[2 + n - 1] - q * h[2 + n - 3]) / h[2 + n - 2];
		}
		if (n == 2)
			goto QADRTIC;
		r = 0.0;
	}
ITERATE:
	for (i = maxiter; i > 0; i--) {

		for (j = 0; j <= n; j++) {				/* BAIRSTOW */
			b[2 + j] = h[2 + j] - p * b[2 + j - 1] - q * b[2 + j - 2];
			c[2 + j] = b[2 + j] - p * c[2 + j - 1] - q * c[2 + j - 2];
		}
		if ((h[2 + n - 1] != 0.0) && (b[2 + n - 1] != 0.0)) {
			if (fabs(h[2 + n - 1] / b[2 + n - 1]) >= K) {
				b[2 + n] = h[2 + n] - q * b[2 + n - 2];
			}
			if (b[2 + n] == 0.0)
				goto QADRTIC;
			if (K < fabs(h[2 + n] / b[2 + n]))
				goto QADRTIC;
		}

		for (j = 0; j <= n; j++) {				/* NEWTON */
			d[2 + j] = h[2 + j] + r * d[2 + j - 1];/* Calculate polynomial at r */
			e[2 + j] = d[2 + j] + r * e[2 + j - 1];/* Calculate derivative at r */
		}
		if (d[2 + n] == 0.0)
			goto LINEAR;
		if (K < fabs(h[2 + n] / d[2 + n]))
			goto LINEAR;

		c[2 + n - 1] = -p * c[2 + n - 2] - q * c[2 + n - 3];
		s = c[2 + n - 2] * c[2 + n - 2] - c[2 + n - 1] * c[2 + n - 3];
		if (s == 0.0) {
			p -= 2.0;
			q *= (q + 1.0);
		} else {
			p += (b[2 + n - 1] * c[2 + n - 2] - b[2 + n] * c[2 + n - 3]) / s;
			q += (-b[2 + n - 1] * c[2 + n - 1] + b[2 + n] * c[2 + n - 2]) / s;
		}
		if (e[2 + n - 1] == 0.0)
			r -= 1.0;							/* Minimum step */
		else
			r -= d[2 + n] / e[2 + n - 1];		/* Newton's iteration */
	}
	ps = pt;
	qs = qt;
	pt = p;
	qt = q;
	if (rev < 0.0)
		K /= 10.0;
	rev = -rev;
	goto REVERSE;

LINEAR:
	if (t < 0)
		r = 1.0 / r;
	n--;
	*u++ = r;
	*v++ = 0.0;
	*conv++ = K;

	for (j = n; j >= 0; j--) {					/* Polynomial deflation by lin-nomial */
		if ((d[2 + j] != 0.0) && (fabs(h[2 + j] / d[2 + j]) < K))
			h[2 + j] = d[2 + j];
		else
			h[2 + j] = 0.0;
	}

	if (n == 0)
		return;
	goto ITERATE;

QADRTIC:
	if (t < 0) {
		p /= q;
		q = 1.0 / q;
	}
	n -= 2;

	if (0.0 < (q - (p * p / 4.0))) {			/* Two complex roots */
		*(u + 1) = *u = -p / 2.0;
		u += 2;
		s = sqrt(q - (p * p / 4.0));
		*v++ = s;
		*v++ = -s;
	} else {									/* Two real roots */
		s = sqrt(((p * p / 4.0)) - q);
		if (p < 0.0)
			*u++ = -p / 2.0 + s;
		else
			*u++ = -p / 2.0 - s;
		*u++ = q / u[-1];
		*v++ = 0.0;
		*v++ = 0.0;
	}
	*conv++ = K;
	*conv++ = K;

	for (j = n; j >= 0; j--) {					/* Polynomial deflation by quadratic */
		if ((b[2 + j] != 0.0) && (fabs(h[2 + j] / b[2 + j]) < K))
			h[2 + j] = b[2 + j];
		else
			h[2 + j] = 0.0;
	}
	goto INIT;
}


#undef MAXN
