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
/// ((0 0 0 T) 0)


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)

/* ====== (1 0 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
PSI_LI int1_t_10_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((1 1 0 T) 8)
/// ((1 1 0 S) 17)

/* ====== (1 1 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
/// (((1 0 0 T) NIL (0 0 0 T) (1 1 0 S) NIL) (1 1 0 T) 1)
PSI_LI int1_t_11_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_s_vercal_11_1 ( coef, w +     5, w +     0, w +     1, w +    17 ) ;
    int1_t_vercal_11_1 ( coef, w +     2, w +     0, w +     0, w +    17, w +     0, w +     8 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((2 0 0 T) 8)
/// ((2 0 0 S) 14)

/* ====== (2 0 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
PSI_LI int1_t_20_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    14 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    14, w +     1, w +     8 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((2 0 0 T) 8)
/// ((2 0 0 S) 14)
/// ((2 1 0 T) 20)
/// ((2 1 0 S) 38)

/* ====== (2 1 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((2 0 0 T) NIL (1 0 0 T) (2 1 0 S) NIL) (2 1 0 T) 1)
PSI_LI int1_t_21_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    14 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    14, w +     1, w +     8 ) ;
    int1_s_vercal_21_1 ( coef, w +    14, w +     0, w +     5, w +    38 ) ;
    int1_t_vercal_21_1 ( coef, w +     8, w +     0, w +     2, w +    38, w +     0, w +    20 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((1 1 0 T) 8)
/// ((1 1 0 S) 17)
/// ((2 0 0 T) 26)
/// ((2 0 0 S) 32)
/// ((2 1 0 T) 38)
/// ((2 1 0 S) 56)
/// ((2 2 0 T) 74)
/// ((2 2 0 S) 110)

/* ====== (2 2 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((2 1 0 S) (2 0 0 S) (1 1 0 S)) (2 2 0 S) 1)
/// (((1 0 0 T) NIL (0 0 0 T) (1 1 0 S) NIL) (1 1 0 T) 1)
/// (((2 0 0 T) NIL (1 0 0 T) (2 1 0 S) NIL) (2 1 0 T) 1)
/// (((2 1 0 T) (2 0 0 T) (1 1 0 T) (2 2 0 S) (2 0 0 S)) (2 2 0 T) 1)
PSI_LI int1_t_22_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    32 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    32, w +     1, w +    26 ) ;
    int1_s_vercal_11_1 ( coef, w +     5, w +     0, w +     1, w +    17 ) ;
    int1_s_vercal_21_1 ( coef, w +    32, w +     0, w +     5, w +    56 ) ;
    int1_s_vercal_22_1 ( coef, w +    56, w +    32, w +    17, w +   110 ) ;
    int1_t_vercal_11_1 ( coef, w +     2, w +     0, w +     0, w +    17, w +     0, w +     8 ) ;
    int1_t_vercal_21_1 ( coef, w +    26, w +     0, w +     2, w +    56, w +     0, w +    38 ) ;
    int1_t_vercal_22_1 ( coef, w +    38, w +    26, w +     8, w +   110, w +    32, w +    74 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((2 0 0 T) 8)
/// ((2 0 0 S) 14)
/// ((3 0 0 T) 20)
/// ((3 0 0 S) 30)

/* ====== (3 0 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((2 0 0 T) (1 0 0 T) NIL (3 0 0 S) (1 0 0 S)) (3 0 0 T) 0)
PSI_LI int1_t_30_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    14 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_s_vercal_30_0 ( coef, w +    14, w +     5, w +     0, w +    30 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    14, w +     1, w +     8 ) ;
    int1_t_vercal_30_0 ( coef, w +     8, w +     2, w +     0, w +    30, w +     5, w +    20 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((2 0 0 T) 8)
/// ((2 0 0 S) 14)
/// ((3 0 0 T) 20)
/// ((3 0 0 S) 30)
/// ((3 1 0 T) 40)
/// ((3 1 0 S) 70)

/* ====== (3 1 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((2 0 0 T) (1 0 0 T) NIL (3 0 0 S) (1 0 0 S)) (3 0 0 T) 0)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
/// (((3 0 0 T) NIL (2 0 0 T) (3 1 0 S) NIL) (3 1 0 T) 1)
PSI_LI int1_t_31_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    14 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_s_vercal_30_0 ( coef, w +    14, w +     5, w +     0, w +    30 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    14, w +     1, w +     8 ) ;
    int1_t_vercal_30_0 ( coef, w +     8, w +     2, w +     0, w +    30, w +     5, w +    20 ) ;
    int1_s_vercal_31_1 ( coef, w +    30, w +     0, w +    14, w +    70 ) ;
    int1_t_vercal_31_1 ( coef, w +    20, w +     0, w +     8, w +    70, w +     0, w +    40 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((2 0 0 T) 8)
/// ((2 0 0 S) 14)
/// ((2 1 0 T) 20)
/// ((2 1 0 S) 38)
/// ((3 0 0 T) 56)
/// ((3 0 0 S) 66)
/// ((3 1 0 T) 76)
/// ((3 1 0 S) 106)
/// ((3 2 0 T) 136)
/// ((3 2 0 S) 196)

/* ====== (3 2 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((2 0 0 T) (1 0 0 T) NIL (3 0 0 S) (1 0 0 S)) (3 0 0 T) 0)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
/// (((3 1 0 S) (3 0 0 S) (2 1 0 S)) (3 2 0 S) 1)
/// (((2 0 0 T) NIL (1 0 0 T) (2 1 0 S) NIL) (2 1 0 T) 1)
/// (((3 0 0 T) NIL (2 0 0 T) (3 1 0 S) NIL) (3 1 0 T) 1)
/// (((3 1 0 T) (3 0 0 T) (2 1 0 T) (3 2 0 S) (3 0 0 S)) (3 2 0 T) 1)
PSI_LI int1_t_32_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    14 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_s_vercal_30_0 ( coef, w +    14, w +     5, w +     0, w +    66 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    14, w +     1, w +     8 ) ;
    int1_t_vercal_30_0 ( coef, w +     8, w +     2, w +     0, w +    66, w +     5, w +    56 ) ;
    int1_s_vercal_21_1 ( coef, w +    14, w +     0, w +     5, w +    38 ) ;
    int1_s_vercal_31_1 ( coef, w +    66, w +     0, w +    14, w +   106 ) ;
    int1_s_vercal_32_1 ( coef, w +   106, w +    66, w +    38, w +   196 ) ;
    int1_t_vercal_21_1 ( coef, w +     8, w +     0, w +     2, w +    38, w +     0, w +    20 ) ;
    int1_t_vercal_31_1 ( coef, w +    56, w +     0, w +     8, w +   106, w +     0, w +    76 ) ;
    int1_t_vercal_32_1 ( coef, w +    76, w +    56, w +    20, w +   196, w +    66, w +   136 ) ;

    return 0 ;
}


/// addr:
/// ((0 0 0 T) 0)
/// ((0 0 0 S) 1)
/// ((1 0 0 T) 2)
/// ((1 0 0 S) 5)
/// ((1 1 0 T) 8)
/// ((1 1 0 S) 17)
/// ((2 0 0 T) 26)
/// ((2 0 0 S) 32)
/// ((2 1 0 T) 38)
/// ((2 1 0 S) 56)
/// ((2 2 0 T) 74)
/// ((2 2 0 S) 110)
/// ((3 0 0 T) 146)
/// ((3 0 0 S) 156)
/// ((3 1 0 T) 166)
/// ((3 1 0 S) 196)
/// ((3 2 0 T) 226)
/// ((3 2 0 S) 286)
/// ((3 3 0 T) 346)
/// ((3 3 0 S) 446)

/* ====== (3 3 0 T) ===== */
/// rec:
/// (((0 0 0 S) NIL NIL) (1 0 0 S) 0)
/// (((1 0 0 S) (0 0 0 S) NIL) (2 0 0 S) 0)
/// (((0 0 0 T) NIL NIL (1 0 0 S) NIL) (1 0 0 T) 0)
/// (((2 0 0 S) (1 0 0 S) NIL) (3 0 0 S) 0)
/// (((1 0 0 T) (0 0 0 T) NIL (2 0 0 S) (0 0 0 S)) (2 0 0 T) 0)
/// (((2 0 0 T) (1 0 0 T) NIL (3 0 0 S) (1 0 0 S)) (3 0 0 T) 0)
/// (((3 0 0 S) NIL (2 0 0 S)) (3 1 0 S) 1)
/// (((1 0 0 S) NIL (0 0 0 S)) (1 1 0 S) 1)
/// (((2 0 0 S) NIL (1 0 0 S)) (2 1 0 S) 1)
/// (((2 1 0 S) (2 0 0 S) (1 1 0 S)) (2 2 0 S) 1)
/// (((3 1 0 S) (3 0 0 S) (2 1 0 S)) (3 2 0 S) 1)
/// (((1 0 0 T) NIL (0 0 0 T) (1 1 0 S) NIL) (1 1 0 T) 1)
/// (((2 0 0 T) NIL (1 0 0 T) (2 1 0 S) NIL) (2 1 0 T) 1)
/// (((3 0 0 T) NIL (2 0 0 T) (3 1 0 S) NIL) (3 1 0 T) 1)
/// (((3 2 0 S) (3 1 0 S) (2 2 0 S)) (3 3 0 S) 1)
/// (((2 1 0 T) (2 0 0 T) (1 1 0 T) (2 2 0 S) (2 0 0 S)) (2 2 0 T) 1)
/// (((3 1 0 T) (3 0 0 T) (2 1 0 T) (3 2 0 S) (3 0 0 S)) (3 2 0 T) 1)
/// (((3 2 0 T) (3 1 0 T) (2 2 0 T) (3 3 0 S) (3 1 0 S)) (3 3 0 T) 1)
PSI_LI int1_t_33_ver ( PSI_LD *coef, PSI_LD *w )
{
    int1_s_vercal_10_0 ( coef, w +     1, w +     0, w +     0, w +     5 ) ;
    int1_s_vercal_20_0 ( coef, w +     5, w +     1, w +     0, w +    32 ) ;
    int1_t_vercal_10_0 ( coef, w +     0, w +     0, w +     0, w +     5, w +     0, w +     2 ) ;
    int1_s_vercal_30_0 ( coef, w +    32, w +     5, w +     0, w +   156 ) ;
    int1_t_vercal_20_0 ( coef, w +     2, w +     0, w +     0, w +    32, w +     1, w +    26 ) ;
    int1_t_vercal_30_0 ( coef, w +    26, w +     2, w +     0, w +   156, w +     5, w +   146 ) ;
    int1_s_vercal_31_1 ( coef, w +   156, w +     0, w +    32, w +   196 ) ;
    int1_s_vercal_11_1 ( coef, w +     5, w +     0, w +     1, w +    17 ) ;
    int1_s_vercal_21_1 ( coef, w +    32, w +     0, w +     5, w +    56 ) ;
    int1_s_vercal_22_1 ( coef, w +    56, w +    32, w +    17, w +   110 ) ;
    int1_s_vercal_32_1 ( coef, w +   196, w +   156, w +    56, w +   286 ) ;
    int1_t_vercal_11_1 ( coef, w +     2, w +     0, w +     0, w +    17, w +     0, w +     8 ) ;
    int1_t_vercal_21_1 ( coef, w +    26, w +     0, w +     2, w +    56, w +     0, w +    38 ) ;
    int1_t_vercal_31_1 ( coef, w +   146, w +     0, w +    26, w +   196, w +     0, w +   166 ) ;
    int1_s_vercal_33_1 ( coef, w +   286, w +   196, w +   110, w +   446 ) ;
    int1_t_vercal_22_1 ( coef, w +    38, w +    26, w +     8, w +   110, w +    32, w +    74 ) ;
    int1_t_vercal_32_1 ( coef, w +   166, w +   146, w +    38, w +   286, w +   156, w +   226 ) ;
    int1_t_vercal_33_1 ( coef, w +   226, w +   166, w +    74, w +   446, w +   196, w +   346 ) ;

    return 0 ;
}
