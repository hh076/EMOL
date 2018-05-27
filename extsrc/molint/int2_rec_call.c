// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@soc.ait.kyushu-u.ac.jp
// ====================================================================

#include <stdio.h>
#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_proto_rec.h"

/* ====== (1 0 0 0) ===== */
int eri_1000_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     2 ) ;

    return 0 ;
}

/* ====== (1 0 1 0) ===== */
int eri_1010_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     3 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_1010_2 ( coef, w +     3, w +     6, w +     0, w +     0, w +     1, w +     9 ) ;

    return 0 ;
}

/* ====== (1 1 0 0) ===== */
int eri_1100_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     3 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_2000_0 ( coef, w +     3, w +     6, w +     0, w +     1, w +     0, w +    18 ) ;

    return 0 ;
}
int eri_1100_hor ( double *coef, double *w )
{
    eri_horcal_1100_1 ( coef, w +    18, w +     3, w +     9 ) ;

    return 0 ;
}

/* ====== (1 1 1 0) ===== */
int eri_1110_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     4 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_1010_2 ( coef, w +     4, w +     7, w +     0, w +     0, w +     1, w +    13 ) ;
    eri_vercal_2000_0 ( coef, w +     4, w +     7, w +     0, w +     1, w +     0, w +    49 ) ;
    eri_vercal_2000_0 ( coef, w +     7, w +    10, w +     1, w +     2, w +     0, w +    55 ) ;
    eri_vercal_2010_2 ( coef, w +    49, w +    55, w +     0, w +     0, w +     7, w +    61 ) ;

    return 0 ;
}
int eri_1110_hor ( double *coef, double *w )
{
    eri_horcal_1110_1 ( coef, w +    61, w +    13, w +    22 ) ;

    return 0 ;
}

/* ====== (1 1 1 1) ===== */
int eri_1111_ver ( double *coef, double *w )
{
    eri_vercal_0010_2 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     5 ) ;
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    17 ) ;
    eri_vercal_1010_2 ( coef, w +     8, w +    11, w +     0, w +     0, w +     1, w +    20 ) ;
    eri_vercal_1010_2 ( coef, w +    11, w +    14, w +     0, w +     0, w +     2, w +    29 ) ;
    eri_vercal_2000_0 ( coef, w +     8, w +    11, w +     0, w +     1, w +     0, w +   146 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     1, w +     2, w +     0, w +   152 ) ;
    eri_vercal_2000_0 ( coef, w +    14, w +    17, w +     2, w +     3, w +     0, w +   158 ) ;
    eri_vercal_1020_2 ( coef, w +    20, w +    29, w +     8, w +    11, w +     5, w +   164 ) ;
    eri_vercal_2010_2 ( coef, w +   146, w +   152, w +     0, w +     0, w +    11, w +   182 ) ;
    eri_vercal_2010_2 ( coef, w +   152, w +   158, w +     0, w +     0, w +    14, w +   200 ) ;
    eri_vercal_2020_2 ( coef, w +   182, w +   200, w +   146, w +   152, w +    29, w +   272 ) ;

    return 0 ;
}
int eri_1111_hor ( double *coef, double *w )
{
    eri_horcal_1120_1 ( coef, w +   272, w +   164, w +   218 ) ;
    eri_horcal_1110_1 ( coef, w +   182, w +    20, w +    38 ) ;
    eri_horcal_1111_3 ( coef, w +   218, w +    38, w +    65 ) ;

    return 0 ;
}

/* ====== (2 0 0 0) ===== */
int eri_2000_ver ( double *coef, double *w )
{
    {
        int i ;
        for ( i = 0 ; i < 3 ; i++ ) {
            fprintf( stderr, "(20,00): (ss,ss): %6d, %28.16e\n", i, w[ i ] ) ;
        }
    }
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     3 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_2000_0 ( coef, w +     3, w +     6, w +     0, w +     1, w +     0, w +     9 ) ;

    return 0 ;
}

/* ====== (2 0 1 0) ===== */
int eri_2010_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     4 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_2000_0 ( coef, w +     4, w +     7, w +     0, w +     1, w +     0, w +    13 ) ;
    eri_vercal_2000_0 ( coef, w +     7, w +    10, w +     1, w +     2, w +     0, w +    19 ) ;
    eri_vercal_2010_2 ( coef, w +    13, w +    19, w +     0, w +     0, w +     7, w +    25 ) ;

    return 0 ;
}

/* ====== (1 1 2 0) ===== */
int eri_1120_ver ( double *coef, double *w )
{
    eri_vercal_0010_2 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     5 ) ;
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    17 ) ;
    eri_vercal_1010_2 ( coef, w +     8, w +    11, w +     0, w +     0, w +     1, w +    20 ) ;
    eri_vercal_1010_2 ( coef, w +    11, w +    14, w +     0, w +     0, w +     2, w +    29 ) ;
    eri_vercal_2000_0 ( coef, w +     8, w +    11, w +     0, w +     1, w +     0, w +    38 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     1, w +     2, w +     0, w +    44 ) ;
    eri_vercal_2000_0 ( coef, w +    14, w +    17, w +     2, w +     3, w +     0, w +    50 ) ;
    eri_vercal_1020_2 ( coef, w +    20, w +    29, w +     8, w +    11, w +     5, w +    56 ) ;
    eri_vercal_2010_2 ( coef, w +    38, w +    44, w +     0, w +     0, w +    11, w +    74 ) ;
    eri_vercal_2010_2 ( coef, w +    44, w +    50, w +     0, w +     0, w +    14, w +    92 ) ;
    eri_vercal_2020_2 ( coef, w +    74, w +    92, w +    38, w +    44, w +    29, w +   164 ) ;

    return 0 ;
}
int eri_1120_hor ( double *coef, double *w )
{
    eri_horcal_1120_1 ( coef, w +   164, w +    56, w +   110 ) ;

    return 0 ;
}

/* ====== (2 0 2 0) ===== */
int eri_2020_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     5 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_1010_2 ( coef, w +     8, w +    11, w +     0, w +     0, w +     2, w +    17 ) ;
    eri_vercal_2000_0 ( coef, w +     5, w +     8, w +     0, w +     1, w +     0, w +    26 ) ;
    eri_vercal_2000_0 ( coef, w +     8, w +    11, w +     1, w +     2, w +     0, w +    32 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     2, w +     3, w +     0, w +    38 ) ;
    eri_vercal_2010_2 ( coef, w +    26, w +    32, w +     0, w +     0, w +     8, w +    44 ) ;
    eri_vercal_2010_2 ( coef, w +    32, w +    38, w +     0, w +     0, w +    11, w +    62 ) ;
    eri_vercal_2020_2 ( coef, w +    44, w +    62, w +    26, w +    32, w +    17, w +    80 ) ;

    return 0 ;
}

/* ====== (2 1 0 0) ===== */
int eri_2100_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     4 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_2000_0 ( coef, w +     4, w +     7, w +     0, w +     1, w +     0, w +    13 ) ;
    eri_vercal_2000_0 ( coef, w +     7, w +    10, w +     1, w +     2, w +     0, w +    19 ) ;
    eri_vercal_3000_0 ( coef, w +    13, w +    19, w +     4, w +     7, w +     0, w +    43 ) ;

    return 0 ;
}
int eri_2100_hor ( double *coef, double *w )
{
    eri_horcal_2100_1 ( coef, w +    43, w +    13, w +    25 ) ;

    return 0 ;
}

/* ====== (2 1 1 0) ===== */
int eri_2110_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     5 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_2000_0 ( coef, w +     5, w +     8, w +     0, w +     1, w +     0, w +    17 ) ;
    eri_vercal_2000_0 ( coef, w +     8, w +    11, w +     1, w +     2, w +     0, w +    23 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     2, w +     3, w +     0, w +    29 ) ;
    eri_vercal_2010_2 ( coef, w +    17, w +    23, w +     0, w +     0, w +     8, w +    35 ) ;
    eri_vercal_3000_0 ( coef, w +    17, w +    23, w +     5, w +     8, w +     0, w +   107 ) ;
    eri_vercal_3000_0 ( coef, w +    23, w +    29, w +     8, w +    11, w +     0, w +   117 ) ;
    eri_vercal_3010_2 ( coef, w +   107, w +   117, w +     0, w +     0, w +    23, w +   127 ) ;

    return 0 ;
}
int eri_2110_hor ( double *coef, double *w )
{
    eri_horcal_2110_1 ( coef, w +   127, w +    35, w +    53 ) ;

    return 0 ;
}

/* ====== (2 1 1 1) ===== */
int eri_2111_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     9 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    12 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    15 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    18 ) ;
    eri_vercal_1010_2 ( coef, w +     9, w +    12, w +     0, w +     0, w +     2, w +    21 ) ;
    eri_vercal_2000_0 ( coef, w +     6, w +     9, w +     0, w +     1, w +     0, w +    30 ) ;
    eri_vercal_2000_0 ( coef, w +     9, w +    12, w +     1, w +     2, w +     0, w +    36 ) ;
    eri_vercal_2000_0 ( coef, w +    12, w +    15, w +     2, w +     3, w +     0, w +    42 ) ;
    eri_vercal_2000_0 ( coef, w +    15, w +    18, w +     3, w +     4, w +     0, w +    48 ) ;
    eri_vercal_2010_2 ( coef, w +    30, w +    36, w +     0, w +     0, w +     9, w +    54 ) ;
    eri_vercal_2010_2 ( coef, w +    36, w +    42, w +     0, w +     0, w +    12, w +    72 ) ;
    eri_vercal_2020_2 ( coef, w +    54, w +    72, w +    30, w +    36, w +    21, w +    90 ) ;
    eri_vercal_3000_0 ( coef, w +    30, w +    36, w +     6, w +     9, w +     0, w +   450 ) ;
    eri_vercal_3000_0 ( coef, w +    36, w +    42, w +     9, w +    12, w +     0, w +   460 ) ;
    eri_vercal_3000_0 ( coef, w +    42, w +    48, w +    12, w +    15, w +     0, w +   470 ) ;
    eri_vercal_3010_2 ( coef, w +   450, w +   460, w +     0, w +     0, w +    36, w +   480 ) ;
    eri_vercal_3010_2 ( coef, w +   460, w +   470, w +     0, w +     0, w +    42, w +   510 ) ;
    eri_vercal_3020_2 ( coef, w +   480, w +   510, w +   450, w +   460, w +    72, w +   540 ) ;

    return 0 ;
}
int eri_2111_hor ( double *coef, double *w )
{
    eri_horcal_2120_1 ( coef, w +   540, w +    90, w +   342 ) ;
    eri_horcal_2110_1 ( coef, w +   480, w +    54, w +   126 ) ;
    eri_horcal_2111_3 ( coef, w +   342, w +   126, w +   180 ) ;

    return 0 ;
}

/* ====== (2 1 2 0) ===== */
int eri_2120_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     9 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    12 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    15 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    18 ) ;
    eri_vercal_1010_2 ( coef, w +     9, w +    12, w +     0, w +     0, w +     2, w +    21 ) ;
    eri_vercal_2000_0 ( coef, w +     6, w +     9, w +     0, w +     1, w +     0, w +    30 ) ;
    eri_vercal_2000_0 ( coef, w +     9, w +    12, w +     1, w +     2, w +     0, w +    36 ) ;
    eri_vercal_2000_0 ( coef, w +    12, w +    15, w +     2, w +     3, w +     0, w +    42 ) ;
    eri_vercal_2000_0 ( coef, w +    15, w +    18, w +     3, w +     4, w +     0, w +    48 ) ;
    eri_vercal_2010_2 ( coef, w +    30, w +    36, w +     0, w +     0, w +     9, w +    54 ) ;
    eri_vercal_2010_2 ( coef, w +    36, w +    42, w +     0, w +     0, w +    12, w +    72 ) ;
    eri_vercal_2020_2 ( coef, w +    54, w +    72, w +    30, w +    36, w +    21, w +    90 ) ;
    eri_vercal_3000_0 ( coef, w +    30, w +    36, w +     6, w +     9, w +     0, w +   234 ) ;
    eri_vercal_3000_0 ( coef, w +    36, w +    42, w +     9, w +    12, w +     0, w +   244 ) ;
    eri_vercal_3000_0 ( coef, w +    42, w +    48, w +    12, w +    15, w +     0, w +   254 ) ;
    eri_vercal_3010_2 ( coef, w +   234, w +   244, w +     0, w +     0, w +    36, w +   264 ) ;
    eri_vercal_3010_2 ( coef, w +   244, w +   254, w +     0, w +     0, w +    42, w +   294 ) ;
    eri_vercal_3020_2 ( coef, w +   264, w +   294, w +   234, w +   244, w +    72, w +   324 ) ;

    return 0 ;
}
int eri_2120_hor ( double *coef, double *w )
{
    eri_horcal_2120_1 ( coef, w +   324, w +    90, w +   126 ) ;

    return 0 ;
}

/* ====== (2 1 2 1) ===== */
int eri_2121_ver ( double *coef, double *w )
{
    eri_vercal_0010_2 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    13 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    16 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    19 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    22 ) ;
    eri_vercal_1000_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +    25 ) ;
    eri_vercal_1010_2 ( coef, w +    13, w +    16, w +     0, w +     0, w +     2, w +    28 ) ;
    eri_vercal_1010_2 ( coef, w +    16, w +    19, w +     0, w +     0, w +     3, w +    37 ) ;
    eri_vercal_2000_0 ( coef, w +    10, w +    13, w +     0, w +     1, w +     0, w +    46 ) ;
    eri_vercal_2000_0 ( coef, w +    13, w +    16, w +     1, w +     2, w +     0, w +    52 ) ;
    eri_vercal_2000_0 ( coef, w +    16, w +    19, w +     2, w +     3, w +     0, w +    58 ) ;
    eri_vercal_2000_0 ( coef, w +    19, w +    22, w +     3, w +     4, w +     0, w +    64 ) ;
    eri_vercal_2000_0 ( coef, w +    22, w +    25, w +     4, w +     5, w +     0, w +    70 ) ;
    eri_vercal_1020_2 ( coef, w +    28, w +    37, w +    13, w +    16, w +     7, w +    76 ) ;
    eri_vercal_2010_2 ( coef, w +    46, w +    52, w +     0, w +     0, w +    13, w +    94 ) ;
    eri_vercal_2010_2 ( coef, w +    52, w +    58, w +     0, w +     0, w +    16, w +   112 ) ;
    eri_vercal_2010_2 ( coef, w +    58, w +    64, w +     0, w +     0, w +    19, w +   130 ) ;
    eri_vercal_2020_2 ( coef, w +    94, w +   112, w +    46, w +    52, w +    28, w +   148 ) ;
    eri_vercal_2020_2 ( coef, w +   112, w +   130, w +    52, w +    58, w +    37, w +   184 ) ;
    eri_vercal_3000_0 ( coef, w +    46, w +    52, w +    10, w +    13, w +     0, w +   652 ) ;
    eri_vercal_3000_0 ( coef, w +    52, w +    58, w +    13, w +    16, w +     0, w +   662 ) ;
    eri_vercal_3000_0 ( coef, w +    58, w +    64, w +    16, w +    19, w +     0, w +   672 ) ;
    eri_vercal_3000_0 ( coef, w +    64, w +    70, w +    19, w +    22, w +     0, w +   682 ) ;
    eri_vercal_3010_2 ( coef, w +   652, w +   662, w +     0, w +     0, w +    52, w +   692 ) ;
    eri_vercal_3010_2 ( coef, w +   662, w +   672, w +     0, w +     0, w +    58, w +   722 ) ;
    eri_vercal_3010_2 ( coef, w +   672, w +   682, w +     0, w +     0, w +    64, w +   752 ) ;
    eri_vercal_2030_2 ( coef, w +   148, w +   184, w +    94, w +   112, w +    76, w +   782 ) ;
    eri_vercal_3020_2 ( coef, w +   692, w +   722, w +   652, w +   662, w +   112, w +   842 ) ;
    eri_vercal_3020_2 ( coef, w +   722, w +   752, w +   662, w +   672, w +   130, w +   902 ) ;
    eri_vercal_3030_2 ( coef, w +   842, w +   902, w +   692, w +   722, w +   184, w +  1142 ) ;

    return 0 ;
}
int eri_2121_hor ( double *coef, double *w )
{
    eri_horcal_2130_1 ( coef, w +  1142, w +   782, w +   962 ) ;
    eri_horcal_2120_1 ( coef, w +   842, w +   148, w +   220 ) ;
    eri_horcal_2121_3 ( coef, w +   962, w +   220, w +   328 ) ;

    return 0 ;
}

/* ====== (2 2 0 0) ===== */
int eri_2200_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     5 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_2000_0 ( coef, w +     5, w +     8, w +     0, w +     1, w +     0, w +    17 ) ;
    eri_vercal_2000_0 ( coef, w +     8, w +    11, w +     1, w +     2, w +     0, w +    23 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     2, w +     3, w +     0, w +    29 ) ;
    eri_vercal_3000_0 ( coef, w +    17, w +    23, w +     5, w +     8, w +     0, w +    89 ) ;
    eri_vercal_3000_0 ( coef, w +    23, w +    29, w +     8, w +    11, w +     0, w +    99 ) ;
    eri_vercal_4000_0 ( coef, w +    89, w +    99, w +    17, w +    23, w +     0, w +   139 ) ;

    return 0 ;
}
int eri_2200_hor ( double *coef, double *w )
{
    eri_horcal_3100_1 ( coef, w +   139, w +    89, w +   109 ) ;
    eri_horcal_2100_1 ( coef, w +    89, w +    17, w +    35 ) ;
    eri_horcal_2200_1 ( coef, w +   109, w +    35, w +    53 ) ;

    return 0 ;
}

/* ====== (2 2 1 0) ===== */
int eri_2210_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     6 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +     9 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    12 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    15 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    18 ) ;
    eri_vercal_2000_0 ( coef, w +     6, w +     9, w +     0, w +     1, w +     0, w +    21 ) ;
    eri_vercal_2000_0 ( coef, w +     9, w +    12, w +     1, w +     2, w +     0, w +    27 ) ;
    eri_vercal_2000_0 ( coef, w +    12, w +    15, w +     2, w +     3, w +     0, w +    33 ) ;
    eri_vercal_2000_0 ( coef, w +    15, w +    18, w +     3, w +     4, w +     0, w +    39 ) ;
    eri_vercal_2010_2 ( coef, w +    21, w +    27, w +     0, w +     0, w +     9, w +    45 ) ;
    eri_vercal_3000_0 ( coef, w +    21, w +    27, w +     6, w +     9, w +     0, w +   225 ) ;
    eri_vercal_3000_0 ( coef, w +    27, w +    33, w +     9, w +    12, w +     0, w +   235 ) ;
    eri_vercal_3000_0 ( coef, w +    33, w +    39, w +    12, w +    15, w +     0, w +   245 ) ;
    eri_vercal_3010_2 ( coef, w +   225, w +   235, w +     0, w +     0, w +    27, w +   255 ) ;
    eri_vercal_4000_0 ( coef, w +   225, w +   235, w +    21, w +    27, w +     0, w +   375 ) ;
    eri_vercal_4000_0 ( coef, w +   235, w +   245, w +    27, w +    33, w +     0, w +   390 ) ;
    eri_vercal_4010_2 ( coef, w +   375, w +   390, w +     0, w +     0, w +   235, w +   405 ) ;

    return 0 ;
}
int eri_2210_hor ( double *coef, double *w )
{
    eri_horcal_3110_1 ( coef, w +   405, w +   255, w +   285 ) ;
    eri_horcal_2110_1 ( coef, w +   255, w +    45, w +    63 ) ;
    eri_horcal_2210_1 ( coef, w +   285, w +    63, w +   117 ) ;

    return 0 ;
}

/* ====== (2 2 1 1) ===== */
int eri_2211_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    13 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    16 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    19 ) ;
    eri_vercal_1000_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +    22 ) ;
    eri_vercal_1010_2 ( coef, w +    10, w +    13, w +     0, w +     0, w +     2, w +    25 ) ;
    eri_vercal_2000_0 ( coef, w +     7, w +    10, w +     0, w +     1, w +     0, w +    34 ) ;
    eri_vercal_2000_0 ( coef, w +    10, w +    13, w +     1, w +     2, w +     0, w +    40 ) ;
    eri_vercal_2000_0 ( coef, w +    13, w +    16, w +     2, w +     3, w +     0, w +    46 ) ;
    eri_vercal_2000_0 ( coef, w +    16, w +    19, w +     3, w +     4, w +     0, w +    52 ) ;
    eri_vercal_2000_0 ( coef, w +    19, w +    22, w +     4, w +     5, w +     0, w +    58 ) ;
    eri_vercal_2010_2 ( coef, w +    34, w +    40, w +     0, w +     0, w +    10, w +    64 ) ;
    eri_vercal_2010_2 ( coef, w +    40, w +    46, w +     0, w +     0, w +    13, w +    82 ) ;
    eri_vercal_2020_2 ( coef, w +    64, w +    82, w +    34, w +    40, w +    25, w +   100 ) ;
    eri_vercal_3000_0 ( coef, w +    34, w +    40, w +     7, w +    10, w +     0, w +   946 ) ;
    eri_vercal_3000_0 ( coef, w +    40, w +    46, w +    10, w +    13, w +     0, w +   956 ) ;
    eri_vercal_3000_0 ( coef, w +    46, w +    52, w +    13, w +    16, w +     0, w +   966 ) ;
    eri_vercal_3000_0 ( coef, w +    52, w +    58, w +    16, w +    19, w +     0, w +   976 ) ;
    eri_vercal_3010_2 ( coef, w +   946, w +   956, w +     0, w +     0, w +    40, w +   986 ) ;
    eri_vercal_3010_2 ( coef, w +   956, w +   966, w +     0, w +     0, w +    46, w +  1016 ) ;
    eri_vercal_3020_2 ( coef, w +   986, w +  1016, w +   946, w +   956, w +    82, w +  1046 ) ;
    eri_vercal_4000_0 ( coef, w +   946, w +   956, w +    34, w +    40, w +     0, w +  1376 ) ;
    eri_vercal_4000_0 ( coef, w +   956, w +   966, w +    40, w +    46, w +     0, w +  1391 ) ;
    eri_vercal_4000_0 ( coef, w +   966, w +   976, w +    46, w +    52, w +     0, w +  1406 ) ;
    eri_vercal_4010_2 ( coef, w +  1376, w +  1391, w +     0, w +     0, w +   956, w +  1421 ) ;
    eri_vercal_4010_2 ( coef, w +  1391, w +  1406, w +     0, w +     0, w +   966, w +  1466 ) ;
    eri_vercal_4020_2 ( coef, w +  1421, w +  1466, w +  1376, w +  1391, w +  1016, w +  1511 ) ;

    return 0 ;
}
int eri_2211_hor ( double *coef, double *w )
{
    eri_horcal_3120_1 ( coef, w +  1511, w +  1046, w +  1196 ) ;
    eri_horcal_2120_1 ( coef, w +  1046, w +   100, w +   190 ) ;
    eri_horcal_3110_1 ( coef, w +  1421, w +   986, w +  1106 ) ;
    eri_horcal_2110_1 ( coef, w +   986, w +    64, w +   136 ) ;
    eri_horcal_2220_1 ( coef, w +  1196, w +   190, w +   730 ) ;
    eri_horcal_2210_1 ( coef, w +  1106, w +   136, w +   298 ) ;
    eri_horcal_2211_3 ( coef, w +   730, w +   298, w +   406 ) ;

    return 0 ;
}

/* ====== (2 2 2 0) ===== */
int eri_2220_ver ( double *coef, double *w )
{
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +     7 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    10 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    13 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    16 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    19 ) ;
    eri_vercal_1000_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +    22 ) ;
    eri_vercal_1010_2 ( coef, w +    10, w +    13, w +     0, w +     0, w +     2, w +    25 ) ;
    eri_vercal_2000_0 ( coef, w +     7, w +    10, w +     0, w +     1, w +     0, w +    34 ) ;
    eri_vercal_2000_0 ( coef, w +    10, w +    13, w +     1, w +     2, w +     0, w +    40 ) ;
    eri_vercal_2000_0 ( coef, w +    13, w +    16, w +     2, w +     3, w +     0, w +    46 ) ;
    eri_vercal_2000_0 ( coef, w +    16, w +    19, w +     3, w +     4, w +     0, w +    52 ) ;
    eri_vercal_2000_0 ( coef, w +    19, w +    22, w +     4, w +     5, w +     0, w +    58 ) ;
    eri_vercal_2010_2 ( coef, w +    34, w +    40, w +     0, w +     0, w +    10, w +    64 ) ;
    eri_vercal_2010_2 ( coef, w +    40, w +    46, w +     0, w +     0, w +    13, w +    82 ) ;
    eri_vercal_2020_2 ( coef, w +    64, w +    82, w +    34, w +    40, w +    25, w +   100 ) ;
    eri_vercal_3000_0 ( coef, w +    34, w +    40, w +     7, w +    10, w +     0, w +   460 ) ;
    eri_vercal_3000_0 ( coef, w +    40, w +    46, w +    10, w +    13, w +     0, w +   470 ) ;
    eri_vercal_3000_0 ( coef, w +    46, w +    52, w +    13, w +    16, w +     0, w +   480 ) ;
    eri_vercal_3000_0 ( coef, w +    52, w +    58, w +    16, w +    19, w +     0, w +   490 ) ;
    eri_vercal_3010_2 ( coef, w +   460, w +   470, w +     0, w +     0, w +    40, w +   500 ) ;
    eri_vercal_3010_2 ( coef, w +   470, w +   480, w +     0, w +     0, w +    46, w +   530 ) ;
    eri_vercal_3020_2 ( coef, w +   500, w +   530, w +   460, w +   470, w +    82, w +   560 ) ;
    eri_vercal_4000_0 ( coef, w +   460, w +   470, w +    34, w +    40, w +     0, w +   800 ) ;
    eri_vercal_4000_0 ( coef, w +   470, w +   480, w +    40, w +    46, w +     0, w +   815 ) ;
    eri_vercal_4000_0 ( coef, w +   480, w +   490, w +    46, w +    52, w +     0, w +   830 ) ;
    eri_vercal_4010_2 ( coef, w +   800, w +   815, w +     0, w +     0, w +   470, w +   845 ) ;
    eri_vercal_4010_2 ( coef, w +   815, w +   830, w +     0, w +     0, w +   480, w +   890 ) ;
    eri_vercal_4020_2 ( coef, w +   845, w +   890, w +   800, w +   815, w +   530, w +   935 ) ;

    return 0 ;
}
int eri_2220_hor ( double *coef, double *w )
{
    eri_horcal_3120_1 ( coef, w +   935, w +   560, w +   620 ) ;
    eri_horcal_2120_1 ( coef, w +   560, w +   100, w +   136 ) ;
    eri_horcal_2220_1 ( coef, w +   620, w +   136, w +   244 ) ;

    return 0 ;
}

/* ====== (2 2 2 1) ===== */
int eri_2221_ver ( double *coef, double *w )
{
    eri_vercal_0010_2 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     8 ) ;
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +    11 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    14 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    17 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    20 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    23 ) ;
    eri_vercal_1000_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +    26 ) ;
    eri_vercal_1000_0 ( coef, w +     6, w +     7, w +     0, w +     0, w +     0, w +    29 ) ;
    eri_vercal_1010_2 ( coef, w +    14, w +    17, w +     0, w +     0, w +     2, w +    32 ) ;
    eri_vercal_1010_2 ( coef, w +    17, w +    20, w +     0, w +     0, w +     3, w +    41 ) ;
    eri_vercal_2000_0 ( coef, w +    11, w +    14, w +     0, w +     1, w +     0, w +    50 ) ;
    eri_vercal_2000_0 ( coef, w +    14, w +    17, w +     1, w +     2, w +     0, w +    56 ) ;
    eri_vercal_2000_0 ( coef, w +    17, w +    20, w +     2, w +     3, w +     0, w +    62 ) ;
    eri_vercal_2000_0 ( coef, w +    20, w +    23, w +     3, w +     4, w +     0, w +    68 ) ;
    eri_vercal_2000_0 ( coef, w +    23, w +    26, w +     4, w +     5, w +     0, w +    74 ) ;
    eri_vercal_2000_0 ( coef, w +    26, w +    29, w +     5, w +     6, w +     0, w +    80 ) ;
    eri_vercal_1020_2 ( coef, w +    32, w +    41, w +    14, w +    17, w +     8, w +    86 ) ;
    eri_vercal_2010_2 ( coef, w +    50, w +    56, w +     0, w +     0, w +    14, w +   104 ) ;
    eri_vercal_2010_2 ( coef, w +    56, w +    62, w +     0, w +     0, w +    17, w +   122 ) ;
    eri_vercal_2010_2 ( coef, w +    62, w +    68, w +     0, w +     0, w +    20, w +   140 ) ;
    eri_vercal_2020_2 ( coef, w +   104, w +   122, w +    50, w +    56, w +    32, w +   158 ) ;
    eri_vercal_2020_2 ( coef, w +   122, w +   140, w +    56, w +    62, w +    41, w +   194 ) ;
    eri_vercal_3000_0 ( coef, w +    50, w +    56, w +    11, w +    14, w +     0, w +  1202 ) ;
    eri_vercal_3000_0 ( coef, w +    56, w +    62, w +    14, w +    17, w +     0, w +  1212 ) ;
    eri_vercal_3000_0 ( coef, w +    62, w +    68, w +    17, w +    20, w +     0, w +  1222 ) ;
    eri_vercal_3000_0 ( coef, w +    68, w +    74, w +    20, w +    23, w +     0, w +  1232 ) ;
    eri_vercal_3000_0 ( coef, w +    74, w +    80, w +    23, w +    26, w +     0, w +  1242 ) ;
    eri_vercal_3010_2 ( coef, w +  1202, w +  1212, w +     0, w +     0, w +    56, w +  1252 ) ;
    eri_vercal_3010_2 ( coef, w +  1212, w +  1222, w +     0, w +     0, w +    62, w +  1282 ) ;
    eri_vercal_3010_2 ( coef, w +  1222, w +  1232, w +     0, w +     0, w +    68, w +  1312 ) ;
    eri_vercal_2030_2 ( coef, w +   158, w +   194, w +   104, w +   122, w +    86, w +  1342 ) ;
    eri_vercal_3020_2 ( coef, w +  1252, w +  1282, w +  1202, w +  1212, w +   122, w +  1402 ) ;
    eri_vercal_3020_2 ( coef, w +  1282, w +  1312, w +  1212, w +  1222, w +   140, w +  1462 ) ;
    eri_vercal_3030_2 ( coef, w +  1402, w +  1462, w +  1252, w +  1282, w +   194, w +  2062 ) ;
    eri_vercal_4000_0 ( coef, w +  1202, w +  1212, w +    50, w +    56, w +     0, w +  2642 ) ;
    eri_vercal_4000_0 ( coef, w +  1212, w +  1222, w +    56, w +    62, w +     0, w +  2657 ) ;
    eri_vercal_4000_0 ( coef, w +  1222, w +  1232, w +    62, w +    68, w +     0, w +  2672 ) ;
    eri_vercal_4000_0 ( coef, w +  1232, w +  1242, w +    68, w +    74, w +     0, w +  2687 ) ;
    eri_vercal_4010_2 ( coef, w +  2642, w +  2657, w +     0, w +     0, w +  1212, w +  2702 ) ;
    eri_vercal_4010_2 ( coef, w +  2657, w +  2672, w +     0, w +     0, w +  1222, w +  2747 ) ;
    eri_vercal_4010_2 ( coef, w +  2672, w +  2687, w +     0, w +     0, w +  1232, w +  2792 ) ;
    eri_vercal_4020_2 ( coef, w +  2702, w +  2747, w +  2642, w +  2657, w +  1282, w +  2837 ) ;
    eri_vercal_4020_2 ( coef, w +  2747, w +  2792, w +  2657, w +  2672, w +  1312, w +  2927 ) ;
    eri_vercal_4030_2 ( coef, w +  2837, w +  2927, w +  2702, w +  2747, w +  1462, w +  3017 ) ;

    return 0 ;
}
int eri_2221_hor ( double *coef, double *w )
{
    eri_horcal_3130_1 ( coef, w +  3017, w +  2062, w +  2342 ) ;
    eri_horcal_2130_1 ( coef, w +  2062, w +  1342, w +  1522 ) ;
    eri_horcal_3120_1 ( coef, w +  2837, w +  1402, w +  2162 ) ;
    eri_horcal_2120_1 ( coef, w +  1402, w +   158, w +   230 ) ;
    eri_horcal_2230_1 ( coef, w +  2342, w +  1522, w +  1702 ) ;
    eri_horcal_2220_1 ( coef, w +  2162, w +   230, w +   338 ) ;
    eri_horcal_2221_3 ( coef, w +  1702, w +   338, w +   554 ) ;

    return 0 ;
}

/* ====== (2 2 2 2) ===== */
int eri_2222_ver ( double *coef, double *w )
{
    eri_vercal_0010_2 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +     9 ) ;
    eri_vercal_0010_2 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    12 ) ;
    eri_vercal_1000_0 ( coef, w +     0, w +     1, w +     0, w +     0, w +     0, w +    15 ) ;
    eri_vercal_1000_0 ( coef, w +     1, w +     2, w +     0, w +     0, w +     0, w +    18 ) ;
    eri_vercal_1000_0 ( coef, w +     2, w +     3, w +     0, w +     0, w +     0, w +    21 ) ;
    eri_vercal_1000_0 ( coef, w +     3, w +     4, w +     0, w +     0, w +     0, w +    24 ) ;
    eri_vercal_1000_0 ( coef, w +     4, w +     5, w +     0, w +     0, w +     0, w +    27 ) ;
    eri_vercal_1000_0 ( coef, w +     5, w +     6, w +     0, w +     0, w +     0, w +    30 ) ;
    eri_vercal_1000_0 ( coef, w +     6, w +     7, w +     0, w +     0, w +     0, w +    33 ) ;
    eri_vercal_1000_0 ( coef, w +     7, w +     8, w +     0, w +     0, w +     0, w +    36 ) ;
    eri_vercal_1010_2 ( coef, w +    18, w +    21, w +     0, w +     0, w +     2, w +    39 ) ;
    eri_vercal_1010_2 ( coef, w +    21, w +    24, w +     0, w +     0, w +     3, w +    48 ) ;
    eri_vercal_1010_2 ( coef, w +    24, w +    27, w +     0, w +     0, w +     4, w +    57 ) ;
    eri_vercal_0020_2 ( coef, w +     9, w +    12, w +     2, w +     3, w +     0, w +    66 ) ;
    eri_vercal_2000_0 ( coef, w +    15, w +    18, w +     0, w +     1, w +     0, w +    72 ) ;
    eri_vercal_2000_0 ( coef, w +    18, w +    21, w +     1, w +     2, w +     0, w +    78 ) ;
    eri_vercal_2000_0 ( coef, w +    21, w +    24, w +     2, w +     3, w +     0, w +    84 ) ;
    eri_vercal_2000_0 ( coef, w +    24, w +    27, w +     3, w +     4, w +     0, w +    90 ) ;
    eri_vercal_2000_0 ( coef, w +    27, w +    30, w +     4, w +     5, w +     0, w +    96 ) ;
    eri_vercal_2000_0 ( coef, w +    30, w +    33, w +     5, w +     6, w +     0, w +   102 ) ;
    eri_vercal_2000_0 ( coef, w +    33, w +    36, w +     6, w +     7, w +     0, w +   108 ) ;
    eri_vercal_1020_2 ( coef, w +    39, w +    48, w +    18, w +    21, w +     9, w +   114 ) ;
    eri_vercal_1020_2 ( coef, w +    48, w +    57, w +    21, w +    24, w +    12, w +   132 ) ;
    eri_vercal_2010_2 ( coef, w +    72, w +    78, w +     0, w +     0, w +    18, w +   150 ) ;
    eri_vercal_2010_2 ( coef, w +    78, w +    84, w +     0, w +     0, w +    21, w +   168 ) ;
    eri_vercal_2010_2 ( coef, w +    84, w +    90, w +     0, w +     0, w +    24, w +   186 ) ;
    eri_vercal_2010_2 ( coef, w +    90, w +    96, w +     0, w +     0, w +    27, w +   204 ) ;
    eri_vercal_2020_2 ( coef, w +   150, w +   168, w +    72, w +    78, w +    39, w +   222 ) ;
    eri_vercal_2020_2 ( coef, w +   168, w +   186, w +    78, w +    84, w +    48, w +   258 ) ;
    eri_vercal_2020_2 ( coef, w +   186, w +   204, w +    84, w +    90, w +    57, w +   294 ) ;
    eri_vercal_3000_0 ( coef, w +    72, w +    78, w +    15, w +    18, w +     0, w +  2598 ) ;
    eri_vercal_3000_0 ( coef, w +    78, w +    84, w +    18, w +    21, w +     0, w +  2608 ) ;
    eri_vercal_3000_0 ( coef, w +    84, w +    90, w +    21, w +    24, w +     0, w +  2618 ) ;
    eri_vercal_3000_0 ( coef, w +    90, w +    96, w +    24, w +    27, w +     0, w +  2628 ) ;
    eri_vercal_3000_0 ( coef, w +    96, w +   102, w +    27, w +    30, w +     0, w +  2638 ) ;
    eri_vercal_3000_0 ( coef, w +   102, w +   108, w +    30, w +    33, w +     0, w +  2648 ) ;
    eri_vercal_1030_2 ( coef, w +   114, w +   132, w +    39, w +    48, w +    66, w +  2658 ) ;
    eri_vercal_3010_2 ( coef, w +  2598, w +  2608, w +     0, w +     0, w +    78, w +  2688 ) ;
    eri_vercal_3010_2 ( coef, w +  2608, w +  2618, w +     0, w +     0, w +    84, w +  2718 ) ;
    eri_vercal_3010_2 ( coef, w +  2618, w +  2628, w +     0, w +     0, w +    90, w +  2748 ) ;
    eri_vercal_3010_2 ( coef, w +  2628, w +  2638, w +     0, w +     0, w +    96, w +  2778 ) ;
    eri_vercal_2030_2 ( coef, w +   222, w +   258, w +   150, w +   168, w +   114, w +  2808 ) ;
    eri_vercal_2030_2 ( coef, w +   258, w +   294, w +   168, w +   186, w +   132, w +  2868 ) ;
    eri_vercal_3020_2 ( coef, w +  2688, w +  2718, w +  2598, w +  2608, w +   168, w +  2928 ) ;
    eri_vercal_3020_2 ( coef, w +  2718, w +  2748, w +  2608, w +  2618, w +   186, w +  2988 ) ;
    eri_vercal_3020_2 ( coef, w +  2748, w +  2778, w +  2618, w +  2628, w +   204, w +  3048 ) ;
    eri_vercal_3030_2 ( coef, w +  2928, w +  2988, w +  2688, w +  2718, w +   258, w +  3648 ) ;
    eri_vercal_3030_2 ( coef, w +  2988, w +  3048, w +  2718, w +  2748, w +   294, w +  3748 ) ;
    eri_vercal_4000_0 ( coef, w +  2598, w +  2608, w +    72, w +    78, w +     0, w +  5408 ) ;
    eri_vercal_4000_0 ( coef, w +  2608, w +  2618, w +    78, w +    84, w +     0, w +  5423 ) ;
    eri_vercal_4000_0 ( coef, w +  2618, w +  2628, w +    84, w +    90, w +     0, w +  5438 ) ;
    eri_vercal_4000_0 ( coef, w +  2628, w +  2638, w +    90, w +    96, w +     0, w +  5453 ) ;
    eri_vercal_4000_0 ( coef, w +  2638, w +  2648, w +    96, w +   102, w +     0, w +  5468 ) ;
    eri_vercal_4010_2 ( coef, w +  5408, w +  5423, w +     0, w +     0, w +  2608, w +  5483 ) ;
    eri_vercal_4010_2 ( coef, w +  5423, w +  5438, w +     0, w +     0, w +  2618, w +  5528 ) ;
    eri_vercal_4010_2 ( coef, w +  5438, w +  5453, w +     0, w +     0, w +  2628, w +  5573 ) ;
    eri_vercal_4010_2 ( coef, w +  5453, w +  5468, w +     0, w +     0, w +  2638, w +  5618 ) ;
    eri_vercal_2040_2 ( coef, w +  2808, w +  2868, w +   222, w +   258, w +  2658, w +  5663 ) ;
    eri_vercal_4020_2 ( coef, w +  5483, w +  5528, w +  5408, w +  5423, w +  2718, w +  5753 ) ;
    eri_vercal_4020_2 ( coef, w +  5528, w +  5573, w +  5423, w +  5438, w +  2748, w +  5843 ) ;
    eri_vercal_4020_2 ( coef, w +  5573, w +  5618, w +  5438, w +  5453, w +  2778, w +  5933 ) ;
    eri_vercal_3040_2 ( coef, w +  3648, w +  3748, w +  2928, w +  2988, w +  2868, w +  6833 ) ;
    eri_vercal_4030_2 ( coef, w +  5753, w +  5843, w +  5483, w +  5528, w +  2988, w +  6983 ) ;
    eri_vercal_4030_2 ( coef, w +  5843, w +  5933, w +  5528, w +  5573, w +  3048, w +  7133 ) ;
    eri_vercal_4040_2 ( coef, w +  6983, w +  7133, w +  5753, w +  5843, w +  3748, w +  7733 ) ;

    return 0 ;
}
int eri_2222_hor ( double *coef, double *w )
{
    eri_horcal_3140_1 ( coef, w +  7733, w +  6833, w +  7283 ) ;
    eri_horcal_2140_1 ( coef, w +  6833, w +  5663, w +  6023 ) ;
    eri_horcal_3130_1 ( coef, w +  6983, w +  3648, w +  5108 ) ;
    eri_horcal_2130_1 ( coef, w +  3648, w +  2808, w +  3108 ) ;
    eri_horcal_3120_1 ( coef, w +  5753, w +  2928, w +  3848 ) ;
    eri_horcal_2120_1 ( coef, w +  2928, w +   222, w +   330 ) ;
    eri_horcal_2240_1 ( coef, w +  7283, w +  6023, w +  6293 ) ;
    eri_horcal_2230_1 ( coef, w +  5108, w +  3108, w +  3288 ) ;
    eri_horcal_2220_1 ( coef, w +  3848, w +   330, w +   438 ) ;
    eri_horcal_2231_3 ( coef, w +  6293, w +  3288, w +  4028 ) ;
    eri_horcal_2221_3 ( coef, w +  3288, w +   438, w +   654 ) ;
    eri_horcal_2222_3 ( coef, w +  4028, w +   654, w +  1302 ) ;

    return 0 ;
}
