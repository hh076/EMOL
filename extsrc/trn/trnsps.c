/* trnsps.f -- translated by f2c (version 20100827).
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

/* *DECK  TRNSPS */

/* Subroutine */ int trnsps_(integer *nr, integer *nc, doublereal *b)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer m, n;
    static doublereal t;
    static integer na, nb, mm;

/*     DIMENSION B(2) */

    /* Parameter adjustments */
    --b;

    /* Function Body */
    if (*nc <= 1) {
	goto L25;
    }
    nb = *nr - 1;
    i__1 = *nc;
    for (m = 2; m <= i__1; ++m) {
	na = m;
	i__2 = m;
	for (n = 2; n <= i__2; ++n) {
	    t = b[na];
	    b[na] = b[nb + n];
	    b[nb + n] = t;
/* L10: */
	    na += *nr;
	}
/* L20: */
	nb += *nr;
    }

L25:
    if (*nr == *nc) {
	goto L50;
    }
    nb = *nc * *nr;
    mm = *nc + 1;
    i__1 = *nr;
    for (m = mm; m <= i__1; ++m) {
	na = m;
	i__2 = *nc;
	for (n = 1; n <= i__2; ++n) {
	    b[na] = b[nb + n];
/* L30: */
	    na += *nr;
	}
/* L40: */
	nb += *nr;
    }
L50:
    return 0;
} /* trnsps_ */

