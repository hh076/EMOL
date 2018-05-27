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
#include "PSI_proto_rec.h"

// ====================================================================
// ====================================================================
// ====================================================================

int eri0_rec ( int inttype,
              double *coef_ver, double *coef_hor, double *w, double **eri_prim ) 
{
    switch ( inttype ) {
    case  1: { *eri_prim = w +      0 ;
               eri_0000_ver ( coef_ver, w ) ; } break ;
    case  2: { *eri_prim = w +      2 ;
               eri_1000_ver ( coef_ver, w ) ; } break ;
    case  3: { *eri_prim = w +      9 ;
               eri_1010_ver ( coef_ver, w ) ; } break ;
    case  4: { *eri_prim = w +      9 ;
               eri_1100_ver ( coef_ver, w ) ;
               eri_1100_hor ( coef_hor, w ) ; } break ;
    case  5: { *eri_prim = w +     22 ;
               eri_1110_ver ( coef_ver, w ) ;
               eri_1110_hor ( coef_hor, w ) ; } break ;
    case  6: { *eri_prim = w +     65 ;
               eri_1111_ver ( coef_ver, w ) ;
               eri_1111_hor ( coef_hor, w ) ; } break ;
    case  7: { *eri_prim = w +      9 ;
               eri_2000_ver ( coef_ver, w ) ; } break ;
    case  8: { *eri_prim = w +     25 ;
               eri_2010_ver ( coef_ver, w ) ; } break ;
    case  9: { *eri_prim = w +    110 ;
               eri_1120_ver ( coef_ver, w ) ;
               eri_1120_hor ( coef_hor, w ) ; } break ;
    case 10: { *eri_prim = w +     80 ;
               eri_2020_ver ( coef_ver, w ) ; } break ;
    case 11: { *eri_prim = w +     25 ;
               eri_2100_ver ( coef_ver, w ) ;
               eri_2100_hor ( coef_hor, w ) ; } break ;
    case 12: { *eri_prim = w +     53 ;
               eri_2110_ver ( coef_ver, w ) ;
               eri_2110_hor ( coef_hor, w ) ; } break ;
    case 13: { *eri_prim = w +    180 ;
               eri_2111_ver ( coef_ver, w ) ;
               eri_2111_hor ( coef_hor, w ) ; } break ;
    case 14: { *eri_prim = w +    126 ;
               eri_2120_ver ( coef_ver, w ) ;
               eri_2120_hor ( coef_hor, w ) ; } break ;
    case 15: { *eri_prim = w +    328 ;
               eri_2121_ver ( coef_ver, w ) ;
               eri_2121_hor ( coef_hor, w ) ; } break ;
    case 16: { *eri_prim = w +     53 ;
               eri_2200_ver ( coef_ver, w ) ;
               eri_2200_hor ( coef_hor, w ) ; } break ;
    case 17: { *eri_prim = w +    117 ;
               eri_2210_ver ( coef_ver, w ) ;
               eri_2210_hor ( coef_hor, w ) ; } break ;
    case 18: { *eri_prim = w +    406 ;
               eri_2211_ver ( coef_ver, w ) ;
               eri_2211_hor ( coef_hor, w ) ; } break ;
    case 19: { *eri_prim = w +    244 ;
               eri_2220_ver ( coef_ver, w ) ;
               eri_2220_hor ( coef_hor, w ) ; } break ;
    case 20: { *eri_prim = w +    554 ;
               eri_2221_ver ( coef_ver, w ) ;
               eri_2221_hor ( coef_hor, w ) ; } break ;
    case 21: { *eri_prim = w +   1302 ;
               eri_2222_ver ( coef_ver, w ) ;
               eri_2222_hor ( coef_hor, w ) ; } break ;
    default: { fprintf ( stderr, "ERI: Not supported: L > 2. \n" ) ; exit ( 1 ) ; }
    }

    return 0 ;
}

PSI_LI eri_0000_ver ( PSI_LD *coef_ver, PSI_LD *w )
{
    return 0 ;
}
