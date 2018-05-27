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

PSI_LI int1_t_rec ( PSI_LI inttype,
              PSI_LD *coef_ver, PSI_LD *w, PSI_LD **t_prim ) 
{
    switch ( inttype ) {
    case  1: { *t_prim = w +      0 ;
               /* int1_t_00_ver ( coef_ver, w ) ; */ } break ;
    case  2: { *t_prim = w +      2 ;
               int1_t_10_ver ( coef_ver, w ) ; } break ;
    case  3: { *t_prim = w +      8 ;
               int1_t_11_ver ( coef_ver, w ) ; } break ;
    case  4: { *t_prim = w +      8 ;
               int1_t_20_ver ( coef_ver, w ) ; } break ;
    case  5: { *t_prim = w +     20 ;
               int1_t_21_ver ( coef_ver, w ) ; } break ;
    case  6: { *t_prim = w +     74 ;
               int1_t_22_ver ( coef_ver, w ) ; } break ;
    case  7: { *t_prim = w +     20 ;
               int1_t_30_ver ( coef_ver, w ) ; } break ;
    case  8: { *t_prim = w +     40 ;
               int1_t_31_ver ( coef_ver, w ) ; } break ;
    case  9: { *t_prim = w +    136 ;
               int1_t_32_ver ( coef_ver, w ) ; } break ;
    case 10: { *t_prim = w +    346 ;
               int1_t_33_ver ( coef_ver, w ) ; } break ;
    default: { fprintf ( stderr, "Not supported: L > 3. \n" ) ; exit ( 1 ) ; }
    }

    return 0 ;
}
