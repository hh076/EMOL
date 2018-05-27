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

PSI_LI int1_s_rec ( PSI_LI inttype,
              PSI_LD *coef_ver, PSI_LD *w, PSI_LD **s_prim ) 
{
    switch ( inttype ) {
    case  1: { *s_prim = w +      0 ;
               /* int1_s_00_ver ( coef_ver, w ) ; */ } break ;
    case  2: { *s_prim = w +      1 ;
               int1_s_10_ver ( coef_ver, w ) ; } break ;
    case  3: { *s_prim = w +      4 ;
               int1_s_11_ver ( coef_ver, w ) ; } break ;
    case  4: { *s_prim = w +      4 ;
               int1_s_20_ver ( coef_ver, w ) ; } break ;
    case  5: { *s_prim = w +     10 ;
               int1_s_21_ver ( coef_ver, w ) ; } break ;
    case  6: { *s_prim = w +     37 ;
               int1_s_22_ver ( coef_ver, w ) ; } break ;
    case  7: { *s_prim = w +     10 ;
               int1_s_30_ver ( coef_ver, w ) ; } break ;
    case  8: { *s_prim = w +     20 ;
               int1_s_31_ver ( coef_ver, w ) ; } break ;
    case  9: { *s_prim = w +     68 ;
               int1_s_32_ver ( coef_ver, w ) ; } break ;
    case 10: { *s_prim = w +    173 ;
               int1_s_33_ver ( coef_ver, w ) ; } break ;
    default: { fprintf ( stderr, "Not supported: L > 3. \n" ) ; exit ( 1 ) ; }
    }

    return 0 ;
}
