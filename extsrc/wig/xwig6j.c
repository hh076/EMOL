#include <stdio.h>
#include <stdlib.h>
#include "common.h"

/// /* Common Block Declarations */
/// 
/// struct {
///     doublereal rt[57], fasq[154], den[57], fac[57], gamma[114], dfac[114];
///     integer factmx, fasqmx, gammax;
/// } comwig_;
/// 
/// #define comwig_1 comwig_

extern comwig_ comwig_1 ;
extern doublereal delta_(integer *, integer *, integer *);
extern doublereal xwig6j_(integer *, integer *, integer *, integer *, 
	                                            integer *, integer *);

int max ( int i, int j )
{
    int n = i ;
    if ( j > n ) { n = j ; }
    return n ;
}
int min ( int i, int j )
{
    int n = i ;
    if ( j < n ) { n = j ; }
    return n ;
}

/* Table of constant values */

static integer c_n1 = -1;

/* *DECK XWIG6J */
doublereal xwig6j_(integer *i1, integer *i2, integer *i3, integer *i4, 
	integer *i5, integer *i6)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static doublereal a, g;
    static integer ia, ib, ic, id, ie, ik, il, im, in, kk, kl, km, iq, ip, ir,
	     jk, iv, jl, jm, jn, lk, ll, lm;

/*     FUNCTION QWIG6J ( RT, FASQ, DEN, FAC ) */
/*     DIMENSION RT(57), FASQ(154), DEN(57), FAC(57) */
/*     SAVE */

    ia = *i1;
    ib = *i2;
    ic = *i3;
    id = *i4;
    ie = *i5;
    iv = *i6;
    if (ia == 1) {
	goto L1;
    }
    if (ib == 1) {
	goto L2;
    }
    if (ic == 1) {
	goto L3;
    }
    if (id == 1) {
	goto L4;
    }
    if (ie == 1) {
	goto L5;
    }
    if (iv == 1) {
	goto L6;
    }
    a = delta_(&ia, &ib, &ic) * delta_(&id, &ie, &ic) * delta_(&id, &ib, &iv) 
	    * delta_(&ia, &ie, &iv);
    if (a <= 0.) {
	goto L20;
    } else {
	goto L110;
    }
L900:
    fprintf ( stderr, "/// ***** WIG6J RANGE EXCEEDS 56  IP=%5d  PAR HAM = %3d%3d%3d%3d%3d%3d ***** ///\n",
              ip,ia,ib,ic,id,ie,iv ) ;
/* 900 WRITE(6,999) IP,IA,IB,IC,ID,IE,IV */
/* 900 WRITE(0,999) IP,IA,IB,IC,ID,IE,IV */
/* 999 FORMAT(///5X,33H***** WIG6J RANGE EXCEEDS 56  IP=,I5,6H   PAR, */
/*    1   3HAM=,6I3,6H ****+///) */
L20:
    ret_val = 0.f;
    goto L1000;
L1:
    if (ib != ic || ie != iv) {
	goto L20;
    }
    if (delta_(&id, &ie, &ic) == 0.f) {
	goto L20;
    }
    i__1 = id + ie + ic;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ib] * 
	    comwig_1.rt[ie]);
    goto L1000;
L2:
    if (ia != ic || id != iv) {
	goto L20;
    }
    if (delta_(&id, &ie, &ic) == 0.f) {
	goto L20;
    }
    i__1 = id + ie + ic;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ia] * 
	    comwig_1.rt[id]);
    goto L1000;
L3:
    if (ia != ib || id != ie) {
	goto L20;
    }
    if (delta_(&ia, &ie, &iv) == 0.f) {
	goto L20;
    }
    i__1 = ia + ie + iv;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ia] * 
	    comwig_1.rt[id]);
    goto L1000;
L4:
    if (iv != ib || ic != ie) {
	goto L20;
    }
    if (delta_(&ia, &ib, &ic) == 0.f) {
	goto L20;
    }
    i__1 = ia + ib + ic;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ib] * 
	    comwig_1.rt[ic]);
    goto L1000;
L5:
    if (iv != ia || ic != id) {
	goto L20;
    }
    if (delta_(&ia, &ib, &ic) == 0.f) {
	goto L20;
    }
    i__1 = ia + ib + ic;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ia] * 
	    comwig_1.rt[ic]);
    goto L1000;
L6:
    if (ie != ia || ib != id) {
	goto L20;
    }
    if (delta_(&ia, &ib, &ic) == 0.f) {
	goto L20;
    }
    i__1 = ia + ib + ic;
    i__2 = i__1 / 2 + 1;
    ret_val = (doublereal) pow_ii(&c_n1, &i__2) / (comwig_1.rt[ia] * 
	    comwig_1.rt[ib]);
    goto L1000;
L110:
    ik = (ia + ib + ic) / 2;
    il = (id + ie + ic) / 2;
    im = (id + ib + iv) / 2;
    in = (ia + ie + iv) / 2;
    kk = (ia + ib + id + ie - 2) / 2;
    kl = (ia + ic + id + iv - 2) / 2;
    km = (ib + ic + ie + iv - 2) / 2;
/* Computing MAX */
    i__1 = max(ik,il), i__1 = max(i__1,im);
    iq = max(i__1,in);
/* Computing MIN */
    i__1 = min(kk,kl);
    ip = min(i__1,km);
    g = 0.f;
    if (ip > 56) {
	goto L900;
    }
    i__1 = ip;
    for (ir = iq; ir <= i__1; ++ir) {
	jk = ir - ik + 1;
	jl = ir - il + 1;
	jm = ir - im + 1;
	jn = ir - in + 1;
	lk = kk - ir + 1;
	ll = kl - ir + 1;
	lm = km - ir + 1;
/* L15: */
	g -= (doublereal) pow_ii(&c_n1, &ir) * comwig_1.fac[ir] / (
		comwig_1.fac[jk - 1] * comwig_1.fac[jl - 1] * comwig_1.fac[jm 
		- 1] * comwig_1.fac[jn - 1] * comwig_1.fac[lk - 1] * 
		comwig_1.fac[ll - 1] * comwig_1.fac[lm - 1]);
    }
    ret_val = g * a / (comwig_1.den[ik] * comwig_1.den[il] * comwig_1.den[im] 
	    * comwig_1.den[in]);
L1000:
/*     write(6,*)'wig6j:',xwig6j,ia,ib,ic,id,ie,iv */
    return ret_val;
/*     ENTRY QWIG6J (RT,FASQ,DEN,FAC) */
/*     RETURN */
} /* xwig6j_ */

/* *DECK XWIG9X */

doublereal xwig9x_(integer *ia, integer *ib, integer *ic, integer *id, 
	integer *ie, integer *iv, integer *ig, integer *ih, integer *ii)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal ret_val;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer kl, ir, ks;

    ret_val = 0.f;
/* Computing MIN */
    i__1 = *ia + *ii, i__2 = *id + *ih, i__1 = min(i__1,i__2), i__2 = *ib + *
	    iv;
    kl = min(i__1,i__2) - 1;
/* Computing MAX */
    i__4 = (i__1 = *ia - *ii, abs(i__1)), i__5 = (i__2 = *id - *ih, abs(i__2))
	    , i__4 = max(i__4,i__5), i__5 = (i__3 = *ib - *iv, abs(i__3));
    ks = max(i__4,i__5) + 1;
    if (ks > kl || (kl - ks) % 2 != 0) {
	return ret_val;
    }
    i__1 = kl;
    for (ir = ks; ir <= i__1; ir += 2) {
/* L1: */
	ret_val += (doublereal) ir * xwig6j_(ia, ib, ic, iv, ii, &ir) * 
		xwig6j_(id, ie, iv, ib, &ir, ih) * xwig6j_(ig, ih, ii, &ir, 
		ia, id);
    }
    ret_val = -ret_val * (doublereal) pow_ii(&c_n1, &ks);
    return ret_val;
} /* xwig9x_ */

