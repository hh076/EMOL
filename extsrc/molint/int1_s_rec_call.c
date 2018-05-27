// ====================================================================
//
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
//
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@soc.ait.kyushu-u.ac.jp
// ====================================================================

#include <stdio.h>
#include <stdlib.h>
#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_proto_int.h"

// ====================================================================
// ====================================================================
// ====================================================================



/// addr:
/// ((0 0 0 S) 0)


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)

/* ====== (1 0 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
PSI_LI int1_s_10_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((1 1 0 S) 4)

/* ====== (1 1 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
PSI_LI int1_s_11_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_11_1 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((2 0 0 S) 4)

/* ====== (2 0 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
PSI_LI int1_s_20_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((2 0 0 S) 4)
/// ((2 1 0 S) 10)

/* ====== (2 1 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
PSI_LI int1_s_21_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_21_1 ( coef, w +     4, w +     0, w +     1, w +    10 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((1 1 0 S) 4)
/// ((2 0 0 S) 13)
/// ((2 1 0 S) 19)
/// ((2 2 0 S) 37)

/* ====== (2 2 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((2 1 0 S) (2 0 0 S) (1 1 0 S)) (2 2 0 S) 1)
PSI_LI int1_s_22_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +    13 ) ;
    int1_s_vercal_11_1 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_21_1 ( coef, w +    13, w +     0, w +     1, w +    19 ) ;
    int1_s_vercal_22_1 ( coef, w +    19, w +    13, w +     4, w +    37 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((2 0 0 S) 4)
/// ((3 0 0 S) 10)

/* ====== (3 0 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
PSI_LI int1_s_30_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_30_0 ( coef, w +     4, w +     1, w +     0, w +    10 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((2 0 0 S) 4)
/// ((3 0 0 S) 10)
/// ((3 1 0 S) 20)

/* ====== (3 1 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
PSI_LI int1_s_31_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_30_0 ( coef, w +     4, w +     1, w +     0, w +    10 ) ;
    int1_s_vercal_31_1 ( coef, w +    10, w +     0, w +     4, w +    20 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((2 0 0 S) 4)
/// ((2 1 0 S) 10)
/// ((3 0 0 S) 28)
/// ((3 1 0 S) 38)
/// ((3 2 0 S) 68)

/* ====== (3 2 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
/// (((3 1 0 S) (3 0 0 S) (2 1 0 S)) (3 2 0 S) 1)
PSI_LI int1_s_32_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_30_0 ( coef, w +     4, w +     1, w +     0, w +    28 ) ;
    int1_s_vercal_21_1 ( coef, w +     4, w +     0, w +     1, w +    10 ) ;
    int1_s_vercal_31_1 ( coef, w +    28, w +     0, w +     4, w +    38 ) ;
    int1_s_vercal_32_1 ( coef, w +    38, w +    28, w +    10, w +    68 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 S) 0)
/// ((1 0 0 S) 1)
/// ((1 1 0 S) 4)
/// ((2 0 0 S) 13)
/// ((2 1 0 S) 19)
/// ((2 2 0 S) 37)
/// ((3 0 0 S) 73)
/// ((3 1 0 S) 83)
/// ((3 2 0 S) 113)
/// ((3 3 0 S) 173)

/* ====== (3 3 0 S) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
/// (((2 1 0 S) (2 0 0 S) (1 1 0 S)) (2 2 0 S) 1)
/// (((3 1 0 S) (3 0 0 S) (2 1 0 S)) (3 2 0 S) 1)
/// (((3 2 0 S) (3 1 0 S) (2 2 0 S)) (3 3 0 S) 1)
PSI_LI int1_s_33_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     1 ) ;
    int1_s_vercal_20_0 ( coef, w +     1, w +     0, w +     0, w +    13 ) ;
    int1_s_vercal_30_0 ( coef, w +    13, w +     1, w +     0, w +    73 ) ;
    int1_s_vercal_11_1 ( coef, w +     1, w +     0, w +     0, w +     4 ) ;
    int1_s_vercal_21_1 ( coef, w +    13, w +     0, w +     1, w +    19 ) ;
    int1_s_vercal_31_1 ( coef, w +    73, w +     0, w +    13, w +    83 ) ;
    int1_s_vercal_22_1 ( coef, w +    19, w +    13, w +     4, w +    37 ) ;
    int1_s_vercal_32_1 ( coef, w +    83, w +    73, w +    19, w +   113 ) ;
    int1_s_vercal_33_1 ( coef, w +   113, w +    83, w +    37, w +   173 ) ;

    return 0 ;
}
