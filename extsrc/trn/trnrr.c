/* tranrr.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

typedef double doublereal ;
typedef int    integer ;
///#include "f2c.h"

/* *DECK  TRNRR  PTRNRR */

/* Subroutine */ int trnrr_(integer *nin, integer *ntl, integer *np, integer *
	mp, doublereal *a, doublereal *b, doublereal *cp, doublereal *cpq)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m, ia, na, nb, nad;
    static doublereal sum;

/* MY */
/*      GO TO 10 */

/*      ENTRY PTRNRR (NP,MP,A,B,CP,CPQ) */
/*      RETURN */
/* MY */
    /* Parameter adjustments */
    --cpq;
    --cp;
    --b;
    --a;

    /* Function Body */
/* L10: */
    nb = 0;
    i__1 = *mp;
    for (m = 1; m <= i__1; ++m) {
	nad = *nin - 1;
	i__2 = *np;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    na = nad;
	    sum = 0.;
	    i__3 = i__;
	    for (j = 1; j <= i__3; ++j) {
/* L16: */
		sum += a[na + j] * cp[nb + j];
	    }
	    if (i__ == *np) {
		goto L30;
	    }
	    na = na + i__ + i__;
	    ia = i__ + 1;
	    i__3 = *np;
	    for (j = ia; j <= i__3; ++j) {
		sum += a[na] * cp[nb + j];
/* L20: */
		na += j;
	    }
L30:
	    cpq[nb + i__] = sum;
/* L32: */
	    nad += i__;
	}
/* L40: */
	nb += *np;
    }

/* L50: */
    na = 0;
    nad = 1;
    i__1 = *mp;
    for (m = 1; m <= i__1; ++m) {
	nb = 0;
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = 0.;
	    i__3 = *np;
	    for (j = 1; j <= i__3; ++j) {
/* L60: */
		sum += cp[na + j] * cpq[nb + j];
	    }
	    b[nad] = sum;
	    ++nad;
	    if (nad > *ntl) {
		goto L90;
	    }
/* L70: */
	    nb += *np;
	}
/* L80: */
	na += *np;
    }
L90:
    return 0;

} /* trnrr_ */

