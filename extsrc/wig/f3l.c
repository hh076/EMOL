#include "common.h"

/* Common Block Declarations */

/// struct {
///     doublereal rt[57], fasq[154], den[57], fac[57], gamma[114], dfac[114];
///     integer factmx, fasqmx, gammax;
/// } comwig_;
/// 
/// #define comwig_1 comwig_

extern comwig_ comwig_1 ;
extern doublereal delta_(integer *, integer *, integer *);

/* *DECK F3L */
doublereal f3l_(integer *ia, integer *ib, integer *ic)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i1, i2, i3, it, iu;

/*                                 CHECKED TASHI 10/25/89 HUCC */
/*     DIMENSION     RT(57),FASQ(154),DEN(57),FAC(57) */

/*     GO TO 77777 */
/*     ENTRY QF3L   (RT,FASQ,DEN,FAC) */
/*     RETURN */
/* 7777 CONTINUE */

    ret_val = delta_(ia, ib, ic);
/*     write ( *, '(a,3i3d16.8)' ) */
/*    & 'ia,ib,ic,dlt',ia,ib,ic,f3l */
    if (ret_val == 0.f) {
	return ret_val;
    }
    it = (*ia + *ib + *ic + 1) / 4;
    iu = *ia + *ib + *ic + 40;
    i1 = (*ia + *ib - *ic + 3) / 4;
    i2 = (*ib + *ic - *ia + 3) / 4;
    i3 = (*ic + *ia - *ib + 3) / 4;
    ret_val = ret_val * comwig_1.fac[it - 1] / (comwig_1.fasq[iu - 1] * 
	    comwig_1.fac[i1 - 1] * comwig_1.fac[i2 - 1] * comwig_1.fac[i3 - 1]
	    );
/*     write ( *, '(a,5i3,5d16.8)' ) */
/*    & 'i1,i2,i3,iu,it',i1,i2,i3,iu,it, */
/*    & fac(it),fasq(iu),fac(i1),fac(i2),fac(i3) */
    return ret_val;
} /* f3l_ */

