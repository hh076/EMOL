// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================
#include <math.h>

#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_proto_int1.h"
#include "PSI_fmt.h"

#define EXP        exp
#define SQRT       sqrt

////////////////////////////////////////////////////////////////////////////
///// Almost All Coefficients for One Integral Calculations ////////////////
////////////////////////////////////////////////////////////////////////////
PSI_LI int1_coefs_stv ( PSI_LI flg_La_gt_Lb,
                        PSI_LD expa, PSI_LD expb,
                        PSI_LD *Ca,  PSI_LD *Cb,  PSI_LD *coefs,
                        PSI_LI *MAXCOMP, PSI_LI *ADDR_FACT )
{
    PSI_LI i ;
    PSI_LD gamm_ab, gamm_aba, gamm_abb ;
    PSI_LD *PA, *PB, *P, *g_2z, *zeta, *zeta2,
           *xi, *xi2, *za, *zb, *xi_za, *xi_zb, *ca, *cb ;

    PA    = coefs + _ADDR_PA_    ;
    PB    = coefs + _ADDR_PB_    ;
/// PC    = coefs + _ADDR_PC_    ;
    P     = coefs + _ADDR_P_     ;
/// ABAB  = coefs + _ADDR_ABAB_  ;
/// PCPC  = coefs + _ADDR_PCPC_  ;
    g_2z  = coefs + _ADDR_G_2Z_  ;
    zeta  = coefs + _ADDR_ZETA_  ;
    zeta2 = coefs + _ADDR_2ZETA_ ;
    xi    = coefs + _ADDR_XI_    ;
    xi2   = coefs + _ADDR_2XI_   ;
    za    = coefs + _ADDR_ZA_   ;
    zb    = coefs + _ADDR_ZB_   ;
    xi_za = coefs + _ADDR_XI_ZA_ ;
    xi_zb = coefs + _ADDR_XI_ZB_ ;

    if ( flg_La_gt_Lb ) {
        *za  = expa ;
        *zb  = expb ;
         ca  = Ca   ;
         cb  = Cb   ;
    } else {
        *za  = expb ;
        *zb  = expa ;
         ca  = Cb   ;
         cb  = Ca   ;
    }
    *zeta    = *za + *zb ;
    *zeta2   = *zeta * 2 ;
    gamm_ab  = 1 / *zeta ;
    gamm_aba = *za * gamm_ab ;
    gamm_abb = *zb * gamm_ab ;
    *g_2z    = gamm_ab * HALF ;
    *xi      = *zb * gamm_aba ;
    *xi2     = *xi * 2;
    *xi_za   = *zb / *zeta ;
    *xi_zb   = *za / *zeta ;
    for ( i = 0 ; i < 3 ; i++ ) {
        P  [ i ] = ( (*za) * ca [ i ] + (*zb) * cb [ i ] ) * gamm_ab ;
        PA [ i ] = P [ i ] - ca [ i ] ;
        PB [ i ] = P [ i ] - cb [ i ] ;
    }
    if ( flg_La_gt_Lb ) {
        ADDR_FACT [ 0 ] = MAXCOMP [ 0 ] ;
        ADDR_FACT [ 1 ] = MAXCOMP [ 1 ] ;
    } else {
        ADDR_FACT [ 0 ] = MAXCOMP [ 1 ] ;
        ADDR_FACT [ 1 ] = MAXCOMP [ 2 ] ;
    }

    return 0 ;
}

