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

/* *DECK DELTA */
doublereal delta_(integer *ia, integer *ib, integer *ic)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static integer i1, i2, i3;

/*                                 CHECKED TASHI 10/25/89 HUCC */
/*     DIMENSION     RT(57),FASQ(154),DEN(57),FAC(57) */
/*     SAVE RT,FASQ,DEN,FAC */

    i1 = *ia + *ib - *ic + 40;
    i2 = *ia - *ib + *ic + 40;
    i3 = -(*ia) + *ib + *ic + 40;
/*     write ( *, '(6i3,3d16.8)' ) */
/*    & ia,ib,ic,i1,i2,i3, */
/*    & fasq(i1),fasq(i2),fasq(i3) */
    ret_val = comwig_1.fasq[i1 - 1] * comwig_1.fasq[i2 - 1] * comwig_1.fasq[
	    i3 - 1];
    return ret_val;
} /* delta_ */

