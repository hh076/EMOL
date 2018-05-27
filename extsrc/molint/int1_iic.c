// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================
#include <math.h>
#include <stdio.h>

#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_proto_int1.h"
#include "PSI_fmt.h"

#define EXP        exp
#define SQRT       sqrt

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
///// Initial Integral Calculations ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

PSI_LI int1_s_iic ( PSI_LD *coefs, PSI_LD *s_ss )
{
    PSI_LD *zeta  = coefs + _ADDR_ZETA_  ;
    PSI_LD *xi    = coefs + _ADDR_XI_    ;
    PSI_LD *ABAB  = coefs + _ADDR_ABAB_  ;

    PSI_LD CP     = PI / *zeta ; 
    *s_ss = CP * SQRT ( CP ) * EXP ( - (*xi) * (*ABAB) ) ;
///fprintf ( stdout, "iic: %26.16e%26.16e%26.16e%26.16e\n",
///                  *s_ss, *zeta, *xi, *ABAB ) ;
///fflush ( stdout ) ;

    return 0 ;
}

PSI_LI int1_t_iic ( PSI_LD *coefs, PSI_LD *t_ss )
{
    PSI_LD *zeta  = coefs + _ADDR_ZETA_  ;
    PSI_LD *xi    = coefs + _ADDR_XI_    ;
    PSI_LD *ABAB  = coefs + _ADDR_ABAB_  ;

    PSI_LD CP     = PI / *zeta ; 
    PSI_LD s_ss   = CP * SQRT ( CP ) * EXP ( - (*xi) * (*ABAB) ) ;
    *t_ss         = (*xi) * ( 3 - 2 * (*xi) * (*ABAB) ) * s_ss ;

    return 0 ;
}


PSI_LI int1_v_iic ( PSI_LI M_value, PSI_LD Z_atom,
                    PSI_LD *coefs,  PSI_LD *v_ss )
{
    PSI_LD *zeta  = coefs + _ADDR_ZETA_  ;
    PSI_LD *xi    = coefs + _ADDR_XI_    ;
    PSI_LD *ABAB  = coefs + _ADDR_ABAB_  ;
    PSI_LD *PCPC  = coefs + _ADDR_PCPC_  ;

    PSI_LD CP     = PI / *zeta ; 
    PSI_LD s_ss   = CP * SQRT ( CP ) * EXP ( - (*xi) * (*ABAB) ) ;
    PSI_LD U      = (*zeta) * (*PCPC) ;
    PSI_LD C      = - 2 * Z_atom * SQRT ( (*zeta) / PI ) * s_ss ; 
    
    fmt( v_ss, M_value, U, C ) ;

    return 0 ;
}
