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
/// ((0 0 0 V) 0)


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((1 0 0 V) 2)

/* ====== (1 0 0 V) ===== */
/// rec:
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
PSI_LI int1_v_10_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     2 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((1 0 0 V) 3)
/// ((1 0 1 V) 6)
/// ((1 1 0 V) 9)

/* ====== (1 1 0 V) ===== */
/// rec:
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((1 0 0 V) (1 0 1 V) NIL NIL (0 0 0 V) (0 0 1 V)) (1 1 0 V) 1)
PSI_LI int1_v_11_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     6 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     3 ) ;
    int1_v_vercal_11_1 ( coef, w +     3, w +     6, w +     0, w +     0, w +     0, w +     1, w +     9 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((1 0 0 V) 3)
/// ((1 0 1 V) 6)
/// ((2 0 0 V) 9)

/* ====== (2 0 0 V) ===== */
/// rec:
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
PSI_LI int1_v_20_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     6 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     3 ) ;
    int1_v_vercal_20_0 ( coef, w +     3, w +     6, w +     0, w +     1, w +     0, w +     0, w +     9 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((1 0 0 V) 4)
/// ((1 0 1 V) 7)
/// ((1 0 2 V) 10)
/// ((2 0 0 V) 13)
/// ((2 0 1 V) 19)
/// ((2 1 0 V) 25)

/* ====== (2 1 0 V) ===== */
/// rec:
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((2 0 0 V) (2 0 1 V) NIL NIL (1 0 0 V) (1 0 1 V)) (2 1 0 V) 1)
PSI_LI int1_v_21_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     7 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     4 ) ;
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    10 ) ;
    int1_v_vercal_20_0 ( coef, w +     7, w +    10, w +     1, w +     2, w +     0, w +     0, w +    19 ) ;
    int1_v_vercal_20_0 ( coef, w +     4, w +     7, w +     0, w +     1, w +     0, w +     0, w +    13 ) ;
    int1_v_vercal_21_1 ( coef, w +    13, w +    19, w +     0, w +     0, w +     4, w +     7, w +    25 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((0 0 4 V) 4)
/// ((1 0 0 V) 5)
/// ((1 0 1 V) 8)
/// ((1 0 2 V) 11)
/// ((1 0 3 V) 14)
/// ((1 1 0 V) 17)
/// ((1 1 1 V) 26)
/// ((2 0 0 V) 35)
/// ((2 0 1 V) 41)
/// ((2 0 2 V) 47)
/// ((2 1 0 V) 53)
/// ((2 1 1 V) 71)
/// ((2 2 0 V) 89)

/* ====== (2 2 0 V) ===== */
/// rec:
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 3 V) (0 0 4 V) NIL NIL NIL NIL) (1 0 3 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((1 0 2 V) (1 0 3 V) (0 0 2 V) (0 0 3 V) NIL NIL) (2 0 2 V) 0)
/// (((1 0 1 V) (1 0 2 V) NIL NIL (0 0 1 V) (0 0 2 V)) (1 1 1 V) 1)
/// (((1 0 0 V) (1 0 1 V) NIL NIL (0 0 0 V) (0 0 1 V)) (1 1 0 V) 1)
/// (((2 0 1 V) (2 0 2 V) NIL NIL (1 0 1 V) (1 0 2 V)) (2 1 1 V) 1)
/// (((2 0 0 V) (2 0 1 V) NIL NIL (1 0 0 V) (1 0 1 V)) (2 1 0 V) 1)
/// (((2 1 0 V) (2 1 1 V) (2 0 0 V) (2 0 1 V) (1 1 0 V) (1 1 1 V)) (2 2 0 V) 1)
PSI_LI int1_v_22_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    11 ) ;
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     8 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     5 ) ;
    int1_v_vercal_10_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +     0, w +    14 ) ;
    int1_v_vercal_20_0 ( coef, w +     8, w +    11, w +     1, w +     2, w +     0, w +     0, w +    41 ) ;
    int1_v_vercal_20_0 ( coef, w +     5, w +     8, w +     0, w +     1, w +     0, w +     0, w +    35 ) ;
    int1_v_vercal_20_0 ( coef, w +    11, w +    14, w +     2, w +     3, w +     0, w +     0, w +    47 ) ;
    int1_v_vercal_11_1 ( coef, w +     8, w +    11, w +     0, w +     0, w +     1, w +     2, w +    26 ) ;
    int1_v_vercal_11_1 ( coef, w +     5, w +     8, w +     0, w +     0, w +     0, w +     1, w +    17 ) ;
    int1_v_vercal_21_1 ( coef, w +    41, w +    47, w +     0, w +     0, w +     8, w +    11, w +    71 ) ;
    int1_v_vercal_21_1 ( coef, w +    35, w +    41, w +     0, w +     0, w +     5, w +     8, w +    53 ) ;
    int1_v_vercal_22_1 ( coef, w +    53, w +    71, w +    35, w +    41, w +    17, w +    26, w +    89 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((1 0 0 V) 4)
/// ((1 0 1 V) 7)
/// ((1 0 2 V) 10)
/// ((2 0 0 V) 13)
/// ((2 0 1 V) 19)
/// ((3 0 0 V) 25)

/* ====== (3 0 0 V) ===== */
/// rec:
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((2 0 0 V) (2 0 1 V) (1 0 0 V) (1 0 1 V) NIL NIL) (3 0 0 V) 0)
PSI_LI int1_v_30_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     7 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     4 ) ;
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    10 ) ;
    int1_v_vercal_20_0 ( coef, w +     7, w +    10, w +     1, w +     2, w +     0, w +     0, w +    19 ) ;
    int1_v_vercal_20_0 ( coef, w +     4, w +     7, w +     0, w +     1, w +     0, w +     0, w +    13 ) ;
    int1_v_vercal_30_0 ( coef, w +    13, w +    19, w +     4, w +     7, w +     0, w +     0, w +    25 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((0 0 4 V) 4)
/// ((1 0 0 V) 5)
/// ((1 0 1 V) 8)
/// ((1 0 2 V) 11)
/// ((1 0 3 V) 14)
/// ((2 0 0 V) 17)
/// ((2 0 1 V) 23)
/// ((2 0 2 V) 29)
/// ((3 0 0 V) 35)
/// ((3 0 1 V) 45)
/// ((3 1 0 V) 55)

/* ====== (3 1 0 V) ===== */
/// rec:
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 3 V) (0 0 4 V) NIL NIL NIL NIL) (1 0 3 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((1 0 2 V) (1 0 3 V) (0 0 2 V) (0 0 3 V) NIL NIL) (2 0 2 V) 0)
/// (((2 0 1 V) (2 0 2 V) (1 0 1 V) (1 0 2 V) NIL NIL) (3 0 1 V) 0)
/// (((2 0 0 V) (2 0 1 V) (1 0 0 V) (1 0 1 V) NIL NIL) (3 0 0 V) 0)
/// (((3 0 0 V) (3 0 1 V) NIL NIL (2 0 0 V) (2 0 1 V)) (3 1 0 V) 1)
PSI_LI int1_v_31_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    11 ) ;
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     8 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     5 ) ;
    int1_v_vercal_10_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +     0, w +    14 ) ;
    int1_v_vercal_20_0 ( coef, w +     8, w +    11, w +     1, w +     2, w +     0, w +     0, w +    23 ) ;
    int1_v_vercal_20_0 ( coef, w +     5, w +     8, w +     0, w +     1, w +     0, w +     0, w +    17 ) ;
    int1_v_vercal_20_0 ( coef, w +    11, w +    14, w +     2, w +     3, w +     0, w +     0, w +    29 ) ;
    int1_v_vercal_30_0 ( coef, w +    23, w +    29, w +     8, w +    11, w +     0, w +     0, w +    45 ) ;
    int1_v_vercal_30_0 ( coef, w +    17, w +    23, w +     5, w +     8, w +     0, w +     0, w +    35 ) ;
    int1_v_vercal_31_1 ( coef, w +    35, w +    45, w +     0, w +     0, w +    17, w +    23, w +    55 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((0 0 4 V) 4)
/// ((0 0 5 V) 5)
/// ((1 0 0 V) 6)
/// ((1 0 1 V) 9)
/// ((1 0 2 V) 12)
/// ((1 0 3 V) 15)
/// ((1 0 4 V) 18)
/// ((2 0 0 V) 21)
/// ((2 0 1 V) 27)
/// ((2 0 2 V) 33)
/// ((2 0 3 V) 39)
/// ((2 1 0 V) 45)
/// ((2 1 1 V) 63)
/// ((3 0 0 V) 81)
/// ((3 0 1 V) 91)
/// ((3 0 2 V) 101)
/// ((3 1 0 V) 111)
/// ((3 1 1 V) 141)
/// ((3 2 0 V) 171)

/* ====== (3 2 0 V) ===== */
/// rec:
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 3 V) (0 0 4 V) NIL NIL NIL NIL) (1 0 3 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 4 V) (0 0 5 V) NIL NIL NIL NIL) (1 0 4 V) 0)
/// (((1 0 2 V) (1 0 3 V) (0 0 2 V) (0 0 3 V) NIL NIL) (2 0 2 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((1 0 3 V) (1 0 4 V) (0 0 3 V) (0 0 4 V) NIL NIL) (2 0 3 V) 0)
/// (((2 0 1 V) (2 0 2 V) (1 0 1 V) (1 0 2 V) NIL NIL) (3 0 1 V) 0)
/// (((2 0 0 V) (2 0 1 V) (1 0 0 V) (1 0 1 V) NIL NIL) (3 0 0 V) 0)
/// (((2 0 2 V) (2 0 3 V) (1 0 2 V) (1 0 3 V) NIL NIL) (3 0 2 V) 0)
/// (((2 0 1 V) (2 0 2 V) NIL NIL (1 0 1 V) (1 0 2 V)) (2 1 1 V) 1)
/// (((2 0 0 V) (2 0 1 V) NIL NIL (1 0 0 V) (1 0 1 V)) (2 1 0 V) 1)
/// (((3 0 1 V) (3 0 2 V) NIL NIL (2 0 1 V) (2 0 2 V)) (3 1 1 V) 1)
/// (((3 0 0 V) (3 0 1 V) NIL NIL (2 0 0 V) (2 0 1 V)) (3 1 0 V) 1)
/// (((3 1 0 V) (3 1 1 V) (3 0 0 V) (3 0 1 V) (2 1 0 V) (2 1 1 V)) (3 2 0 V) 1)
PSI_LI int1_v_32_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    12 ) ;
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +     9 ) ;
    int1_v_vercal_10_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +     0, w +    15 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     6 ) ;
    int1_v_vercal_10_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +     0, w +    18 ) ;
    int1_v_vercal_20_0 ( coef, w +    12, w +    15, w +     2, w +     3, w +     0, w +     0, w +    33 ) ;
    int1_v_vercal_20_0 ( coef, w +     9, w +    12, w +     1, w +     2, w +     0, w +     0, w +    27 ) ;
    int1_v_vercal_20_0 ( coef, w +     6, w +     9, w +     0, w +     1, w +     0, w +     0, w +    21 ) ;
    int1_v_vercal_20_0 ( coef, w +    15, w +    18, w +     3, w +     4, w +     0, w +     0, w +    39 ) ;
    int1_v_vercal_30_0 ( coef, w +    27, w +    33, w +     9, w +    12, w +     0, w +     0, w +    91 ) ;
    int1_v_vercal_30_0 ( coef, w +    21, w +    27, w +     6, w +     9, w +     0, w +     0, w +    81 ) ;
    int1_v_vercal_30_0 ( coef, w +    33, w +    39, w +    12, w +    15, w +     0, w +     0, w +   101 ) ;
    int1_v_vercal_21_1 ( coef, w +    27, w +    33, w +     0, w +     0, w +     9, w +    12, w +    63 ) ;
    int1_v_vercal_21_1 ( coef, w +    21, w +    27, w +     0, w +     0, w +     6, w +     9, w +    45 ) ;
    int1_v_vercal_31_1 ( coef, w +    91, w +   101, w +     0, w +     0, w +    27, w +    33, w +   141 ) ;
    int1_v_vercal_31_1 ( coef, w +    81, w +    91, w +     0, w +     0, w +    21, w +    27, w +   111 ) ;
    int1_v_vercal_32_1 ( coef, w +   111, w +   141, w +    81, w +    91, w +    45, w +    63, w +   171 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 V) 0)
/// ((0 0 1 V) 1)
/// ((0 0 2 V) 2)
/// ((0 0 3 V) 3)
/// ((0 0 4 V) 4)
/// ((0 0 5 V) 5)
/// ((0 0 6 V) 6)
/// ((1 0 0 V) 7)
/// ((1 0 1 V) 10)
/// ((1 0 2 V) 13)
/// ((1 0 3 V) 16)
/// ((1 0 4 V) 19)
/// ((1 0 5 V) 22)
/// ((1 1 0 V) 25)
/// ((1 1 1 V) 34)
/// ((1 1 2 V) 43)
/// ((2 0 0 V) 52)
/// ((2 0 1 V) 58)
/// ((2 0 2 V) 64)
/// ((2 0 3 V) 70)
/// ((2 0 4 V) 76)
/// ((2 1 0 V) 82)
/// ((2 1 1 V) 100)
/// ((2 1 2 V) 118)
/// ((2 2 0 V) 136)
/// ((2 2 1 V) 172)
/// ((3 0 0 V) 208)
/// ((3 0 1 V) 218)
/// ((3 0 2 V) 228)
/// ((3 0 3 V) 238)
/// ((3 1 0 V) 248)
/// ((3 1 1 V) 278)
/// ((3 1 2 V) 308)
/// ((3 2 0 V) 338)
/// ((3 2 1 V) 398)
/// ((3 3 0 V) 458)

/* ====== (3 3 0 V) ===== */
/// rec:
/// (((0 0 3 V) (0 0 4 V) NIL NIL NIL NIL) (1 0 3 V) 0)
/// (((0 0 2 V) (0 0 3 V) NIL NIL NIL NIL) (1 0 2 V) 0)
/// (((0 0 1 V) (0 0 2 V) NIL NIL NIL NIL) (1 0 1 V) 0)
/// (((0 0 4 V) (0 0 5 V) NIL NIL NIL NIL) (1 0 4 V) 0)
/// (((0 0 0 V) (0 0 1 V) NIL NIL NIL NIL) (1 0 0 V) 0)
/// (((0 0 5 V) (0 0 6 V) NIL NIL NIL NIL) (1 0 5 V) 0)
/// (((1 0 2 V) (1 0 3 V) (0 0 2 V) (0 0 3 V) NIL NIL) (2 0 2 V) 0)
/// (((1 0 1 V) (1 0 2 V) (0 0 1 V) (0 0 2 V) NIL NIL) (2 0 1 V) 0)
/// (((1 0 3 V) (1 0 4 V) (0 0 3 V) (0 0 4 V) NIL NIL) (2 0 3 V) 0)
/// (((1 0 0 V) (1 0 1 V) (0 0 0 V) (0 0 1 V) NIL NIL) (2 0 0 V) 0)
/// (((1 0 4 V) (1 0 5 V) (0 0 4 V) (0 0 5 V) NIL NIL) (2 0 4 V) 0)
/// (((2 0 2 V) (2 0 3 V) (1 0 2 V) (1 0 3 V) NIL NIL) (3 0 2 V) 0)
/// (((2 0 1 V) (2 0 2 V) (1 0 1 V) (1 0 2 V) NIL NIL) (3 0 1 V) 0)
/// (((2 0 0 V) (2 0 1 V) (1 0 0 V) (1 0 1 V) NIL NIL) (3 0 0 V) 0)
/// (((2 0 3 V) (2 0 4 V) (1 0 3 V) (1 0 4 V) NIL NIL) (3 0 3 V) 0)
/// (((1 0 2 V) (1 0 3 V) NIL NIL (0 0 2 V) (0 0 3 V)) (1 1 2 V) 1)
/// (((1 0 1 V) (1 0 2 V) NIL NIL (0 0 1 V) (0 0 2 V)) (1 1 1 V) 1)
/// (((2 0 2 V) (2 0 3 V) NIL NIL (1 0 2 V) (1 0 3 V)) (2 1 2 V) 1)
/// (((2 0 1 V) (2 0 2 V) NIL NIL (1 0 1 V) (1 0 2 V)) (2 1 1 V) 1)
/// (((1 0 0 V) (1 0 1 V) NIL NIL (0 0 0 V) (0 0 1 V)) (1 1 0 V) 1)
/// (((2 0 0 V) (2 0 1 V) NIL NIL (1 0 0 V) (1 0 1 V)) (2 1 0 V) 1)
/// (((3 0 1 V) (3 0 2 V) NIL NIL (2 0 1 V) (2 0 2 V)) (3 1 1 V) 1)
/// (((3 0 0 V) (3 0 1 V) NIL NIL (2 0 0 V) (2 0 1 V)) (3 1 0 V) 1)
/// (((3 0 2 V) (3 0 3 V) NIL NIL (2 0 2 V) (2 0 3 V)) (3 1 2 V) 1)
/// (((2 1 1 V) (2 1 2 V) (2 0 1 V) (2 0 2 V) (1 1 1 V) (1 1 2 V)) (2 2 1 V) 1)
/// (((2 1 0 V) (2 1 1 V) (2 0 0 V) (2 0 1 V) (1 1 0 V) (1 1 1 V)) (2 2 0 V) 1)
/// (((3 1 1 V) (3 1 2 V) (3 0 1 V) (3 0 2 V) (2 1 1 V) (2 1 2 V)) (3 2 1 V) 1)
/// (((3 1 0 V) (3 1 1 V) (3 0 0 V) (3 0 1 V) (2 1 0 V) (2 1 1 V)) (3 2 0 V) 1)
/// (((3 2 0 V) (3 2 1 V) (3 1 0 V) (3 1 1 V) (2 2 0 V) (2 2 1 V)) (3 3 0 V) 1)
PSI_LI int1_v_33_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_v_vercal_10_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +     0, w +    16 ) ;
    int1_v_vercal_10_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     0, w +    13 ) ;
    int1_v_vercal_10_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     0, w +    10 ) ;
    int1_v_vercal_10_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +     0, w +    19 ) ;
    int1_v_vercal_10_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     0, w +     7 ) ;
    int1_v_vercal_10_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +     0, w +    22 ) ;
    int1_v_vercal_20_0 ( coef, w +    13, w +    16, w +     2, w +     3, w +     0, w +     0, w +    64 ) ;
    int1_v_vercal_20_0 ( coef, w +    10, w +    13, w +     1, w +     2, w +     0, w +     0, w +    58 ) ;
    int1_v_vercal_20_0 ( coef, w +    16, w +    19, w +     3, w +     4, w +     0, w +     0, w +    70 ) ;
    int1_v_vercal_20_0 ( coef, w +     7, w +    10, w +     0, w +     1, w +     0, w +     0, w +    52 ) ;
    int1_v_vercal_20_0 ( coef, w +    19, w +    22, w +     4, w +     5, w +     0, w +     0, w +    76 ) ;
    int1_v_vercal_30_0 ( coef, w +    64, w +    70, w +    13, w +    16, w +     0, w +     0, w +   228 ) ;
    int1_v_vercal_30_0 ( coef, w +    58, w +    64, w +    10, w +    13, w +     0, w +     0, w +   218 ) ;
    int1_v_vercal_30_0 ( coef, w +    52, w +    58, w +     7, w +    10, w +     0, w +     0, w +   208 ) ;
    int1_v_vercal_30_0 ( coef, w +    70, w +    76, w +    16, w +    19, w +     0, w +     0, w +   238 ) ;
    int1_v_vercal_11_1 ( coef, w +    13, w +    16, w +     0, w +     0, w +     2, w +     3, w +    43 ) ;
    int1_v_vercal_11_1 ( coef, w +    10, w +    13, w +     0, w +     0, w +     1, w +     2, w +    34 ) ;
    int1_v_vercal_21_1 ( coef, w +    64, w +    70, w +     0, w +     0, w +    13, w +    16, w +   118 ) ;
    int1_v_vercal_21_1 ( coef, w +    58, w +    64, w +     0, w +     0, w +    10, w +    13, w +   100 ) ;
    int1_v_vercal_11_1 ( coef, w +     7, w +    10, w +     0, w +     0, w +     0, w +     1, w +    25 ) ;
    int1_v_vercal_21_1 ( coef, w +    52, w +    58, w +     0, w +     0, w +     7, w +    10, w +    82 ) ;
    int1_v_vercal_31_1 ( coef, w +   218, w +   228, w +     0, w +     0, w +    58, w +    64, w +   278 ) ;
    int1_v_vercal_31_1 ( coef, w +   208, w +   218, w +     0, w +     0, w +    52, w +    58, w +   248 ) ;
    int1_v_vercal_31_1 ( coef, w +   228, w +   238, w +     0, w +     0, w +    64, w +    70, w +   308 ) ;
    int1_v_vercal_22_1 ( coef, w +   100, w +   118, w +    58, w +    64, w +    34, w +    43, w +   172 ) ;
    int1_v_vercal_22_1 ( coef, w +    82, w +   100, w +    52, w +    58, w +    25, w +    34, w +   136 ) ;
    int1_v_vercal_32_1 ( coef, w +   278, w +   308, w +   218, w +   228, w +   100, w +   118, w +   398 ) ;
    int1_v_vercal_32_1 ( coef, w +   248, w +   278, w +   208, w +   218, w +    82, w +   100, w +   338 ) ;
    int1_v_vercal_33_1 ( coef, w +   338, w +   398, w +   248, w +   278, w +   136, w +   172, w +   458 ) ;

    return 0 ;
}
