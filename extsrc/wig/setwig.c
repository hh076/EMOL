#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"

/// /* Common Block Declarations */
/// 
/// struct {
///     doublereal rt[57], fasq[154], den[57], fac[57], gamma[114], dfac[114];
///     integer factmx, fasqmx, gammax;
/// } comwig_;
/// 
/// #define comwig_1 comwig_

comwig_ comwig_1 ;

integer setwig_(void)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    double sqrt(doublereal), atan(doublereal);

    /* Local variables */
    static integer i__, k;

/*     GENERATED BY FPL90 V(AIX OS   Dec23 08)        NEWS OS  Dec23 08 */
/* input */
/* output */
    ret_val = 0;
    for (k = 1; k <= 154; ++k) {
	comwig_1.fasq[k - 1] = 0.f;
/* L1001: */
    }
    comwig_1.fac[0] = 1.;
    for (k = 1; k <= 57; ++k) {
	if (! (k != 1)) {
	    goto L1004;
	}
	comwig_1.fac[k - 1] = comwig_1.fac[k - 2] * (k - 1);
L1004:
	comwig_1.den[k - 1] = sqrt(comwig_1.fac[k - 1]);
	comwig_1.fasq[(k << 1) + 38] = comwig_1.den[k - 1];
	comwig_1.rt[k - 1] = sqrt((doublereal) (k - 1));
/* L1002: */
    }
/* L1003: */
    comwig_1.gamma[0] = 1.;
    comwig_1.gamma[1] = 1.;
    for (i__ = 1; i__ <= 112; ++i__) {
	comwig_1.gamma[i__ + 1] = comwig_1.gamma[i__ - 1] * i__ / 2;
/* L1005: */
    }
    comwig_1.gamma[0] = sqrt(atan(comwig_1.gamma[0])) * 2;
    for (i__ = 1; i__ <= 112; i__ += 2) {
	comwig_1.gamma[i__ + 1] *= comwig_1.gamma[0];
/* L1006: */
    }
/*     write ( *, * ) 'fac,fasq,rt,den' */
/*     do i=1,57 */
/*         write ( *, '(i5,5d16.8)') */
/*    &    i, fac(i), fasq(i), sqrt(fac(i)),rt(i),den(i) */
/*     enddo */
/*     write ( *, * ) 'gamma' */
/*     do i=1,57 */
/*     enddo */
/* L1000: */
    return ret_val;
} /* setwig_ */

