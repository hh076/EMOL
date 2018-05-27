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

PSI_LI int1_v_rec ( PSI_LI inttype,
              PSI_LD *coef_ver, PSI_LD *w, PSI_LD **v_prim ) 
{
    switch ( inttype ) {
    case  1: { *v_prim = w +      0 ;
               /* int1_v_00_ver ( coef_ver, w ) ; */ } break ;
    case  2: { *v_prim = w +      2 ;
               int1_v_10_ver ( coef_ver, w ) ; } break ;
    case  3: { *v_prim = w +      9 ;
               int1_v_11_ver ( coef_ver, w ) ; } break ;
    case  4: { *v_prim = w +      9 ;
               int1_v_20_ver ( coef_ver, w ) ; } break ;
    case  5: { *v_prim = w +     25 ;
               int1_v_21_ver ( coef_ver, w ) ; } break ;
    case  6: { *v_prim = w +     89 ;
               int1_v_22_ver ( coef_ver, w ) ; } break ;
    case  7: { *v_prim = w +     25 ;
               int1_v_30_ver ( coef_ver, w ) ; } break ;
    case  8: { *v_prim = w +     55 ;
               int1_v_31_ver ( coef_ver, w ) ; } break ;
    case  9: { *v_prim = w +    171 ;
               int1_v_32_ver ( coef_ver, w ) ; } break ;
    case 10: { *v_prim = w +    458 ;
               int1_v_33_ver ( coef_ver, w ) ; } break ;
    default: { fprintf ( stderr, "Not supported: L > 3. \n" ) ; exit ( 1 ) ; }
    }

    return 0 ;
}
