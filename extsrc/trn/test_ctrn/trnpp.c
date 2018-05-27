/* trnpp.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdio.h>
typedef double doublereal ;
typedef int    integer ;
///#include "f2c.h"

/* *DECK  TRNPP  PTRNPP */
/* MY */
/* Subroutine */ int trnpp_(integer *nin, integer *nout, integer *np, integer 
	*mp, doublereal *a, doublereal *b, integer *ndd, doublereal *cp, 
	doublereal *cpq)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m, ia, na, nb, nad;
    static doublereal sum;

///    {   int ii ;
///        fprintf ( stdout, "nin, nout, np, mp: %5d%5d%5d%5d\n", *nin, *nout, *np, *mp ) ;
///        for ( ii = 0 ; ii < ((*np)*(*np+1))/2 ; ii++ ) {
///            fprintf ( stdout, "  a: %5d%16.2e\n", ii, a[ ii ] ) ;
///        }
///        for ( ii = 0 ; ii < (*np)*(*np) ; ii++ ) {
///            fprintf ( stdout, " cp: %5d%16.2e\n", ii, cp[ ii ] ) ;
///        }
///    }


/* MY */
/*     DIMENSION A(2),B(2),CP(2),CPQ(2) */
/* MY */
/*      GO TO 10 */

/*      ENTRY PTRNPP (NP,MP,A,B,NDD,CP,CPQ) */
/*      RETURN */
/* MY */
    /* Parameter adjustments */
    --cpq;
    --cp;
    --b;
    --a;

///    {   int ii ;
///        fprintf ( stdout, "nin, nout, np, mp: %5d%5d%5d%5d\n", *nin, *nout, *np, *mp ) ;
///        for ( ii = 0 ; ii < ((*np)*(*np+1))/2+1 ; ii++ ) {
///            fprintf ( stdout, "  a: %5d%16.2e\n", ii, a[ ii ] ) ;
///        }
///        for ( ii = 0 ; ii < (*np)*(*np)+1 ; ii++ ) {
///            fprintf ( stdout, " cp: %5d%16.2e\n", ii, cp[ ii ] ) ;
///        }
///    }

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
		sum += a[na + j ] * cp[nb + j ];
///                fprintf ( stdout, "sum-add1: %3d%3d, %3d%3d: %16.2e,%16.2e,%16.2e\n",
///		          na, na + j -1, nb, nb + j -1, a[na + j ] * cp[nb + j ], sum ) ;
	    }
	    if (i__ == *np) {
		goto L30;
	    }
	    na = na + i__ + i__;
	    ia = i__ + 1;
	    i__3 = *np;
	    for (j = ia; j <= i__3; ++j) {
		sum += a[na ] * cp[nb + j ];
///              fprintf ( stdout, "sum-add2: %3d%3d, %3d%3d: %16.2e,%16.2e,%16.2e\n",
///		          na, na -1, nb, nb + j -1, a[na ] * cp[nb + j ], sum ) ;
/* L20: */
		na += j;
	    }
L30:
///            fprintf ( stdout, "cpq-add: %3d%3d: %16.2e\n", nb, nb+i__-1, sum ) ;
	    cpq[nb + i__ ] = sum;
/* L32: */
	    nad += i__;
	}
/* L40: */
	nb += *np;
    }

///    {   int ii ;
///        for ( ii = 1 ; ii <= (*np)*(*np) ; ii++ ) {
///            fprintf ( stdout, "cpq: %5d%16.2e\n", ii, cpq[ ii ] ) ;
///        }
///    }

/* L50: */
    na = 0;
    nad = *nout;
    i__1 = *mp;
    for (m = 1; m <= i__1; ++m) {
	nb = 0;
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum = 0.;
	    i__3 = *np;
	    for (j = 1; j <= i__3; ++j) {
/* L60: */
		sum += cp[na + j ] * cpq[nb + j ];
	    }
	    b[nad ] = sum;
	    nad += *ndd;
/* L70: */
	    nb += *np;
	}
/* L80: */
	na += *np;
    }
/* L90: */

///    {   int ii ;
///        for ( ii = 1 ; ii <= ((*np)*(*np+1))/2 ; ii++ ) {
///            fprintf ( stdout, "  b: %5d%16.2e\n", ii, b[ ii ] ) ;
///        }
///    }
    return 0;

} /* trnpp_ */

