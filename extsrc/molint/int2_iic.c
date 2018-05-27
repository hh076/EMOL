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
#include "PSI_proto_int2.h"
#include "PSI_fmt.h"

#define EXP        exp
#define SQRT       sqrt

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///// SSSS Calculation /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
PSI_LI eri0_ssss ( PSI_LI M_value, PSI_LD a, PSI_LD b, PSI_LD c, PSI_LD d,
                   PSI_LD *abcdac, PSI_LD *coef_ver, PSI_LD *ssss )
{
    PSI_LD Cv ;
    PSI_LD Const, ABAB, CDCD, ACAC, BAAC, ABCD, CDAC ;
    PSI_LD Gab, Gcd, Ovchab, Ovchcd, Kab, Kcd, Gabcd, den ;
    PSI_LD rho, Gab_a, Gab_b, Gcd_c, Gcd_d, PM, PMC, QM, QMC, T ;
    PSI_LD *C_ab_1, *C_ab_2, *C_ab_3, *C_ab_4, *C_ab_5 ;
    PSI_LD *C_cd_1, *C_cd_2, *C_cd_3, *C_cd_4, *C_cd_5 ;
    PSI_LD *AB, *CD, *AC ;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
    C_ab_1   = coef_ver +  0 ;
    C_ab_2   = coef_ver +  3 ;
    C_ab_3   = coef_ver +  6 ;
    C_ab_4   = coef_ver +  7 ;
    C_ab_5   = coef_ver +  8 ;
    C_cd_1   = coef_ver +  9 ;
    C_cd_2   = coef_ver + 12 ;
    C_cd_3   = coef_ver + 15 ;
    C_cd_4   = coef_ver + 16 ;
    C_cd_5   = coef_ver + 17 ;
    AB       = abcdac   +  0 ;
    CD       = abcdac   +  6 ;
    AC       = abcdac   + 12 ;
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
    Const =  2 * PI * PI * sqrt ( PI ) ;
    ABAB  =  AB [ 0 ] * AB [ 0 ] + AB [ 1 ] * AB [ 1 ]
           + AB [ 2 ] * AB [ 2 ] ;
    CDCD  =  CD [ 0 ] * CD [ 0 ] + CD [ 1 ] * CD [ 1 ]
           + CD [ 2 ] * CD [ 2 ] ;
    ACAC  =  AC [ 0 ] * AC [ 0 ] + AC [ 1 ] * AC [ 1 ]
           + AC [ 2 ] * AC [ 2 ] ;
    BAAC  = - (
             AB [ 0 ] * AC [ 0 ] + AB [ 1 ] * AC [ 1 ]
           + AB [ 2 ] * AC [ 2 ]
              ) ;
    ABCD  =  AB [ 0 ] * CD [ 0 ] + AB [ 1 ] * CD [ 1 ]
           + AB [ 2 ] * CD [ 2 ] ;
    CDAC  =  CD [ 0 ] * AC [ 0 ] + CD [ 1 ] * AC [ 1 ]
           + CD [ 2 ] * AC [ 2 ] ;
    Gab    = 1 / ( a + b ) ;
    Gcd    = 1 / ( c + d ) ;
    Ovchab = a * b * Gab * ABAB ;
    Ovchcd = c * d * Gcd * CDCD ;
    Kab    = Gab * EXP ( - Ovchab ) ;
///fprintf ( stdout, "0: a, b, Gab, ABAB, Ov: %12.4e%12.4e%12.4e%12.4e%12.4e\n",
///                   a, b, Gab, ABAB, Ovchab ) ;
///fprintf ( stdout, "0: Kab: %12.4e, %12.4e, %12.4e\n",
///                   Kab, EXP ( - Ovchab ), Gab * EXP ( - Ovchab ) ) ;
    Kcd    = Gcd * EXP ( - Ovchcd ) ;
    Gabcd  = 1 / ( a + b + c + d ) ;
    den    = SQRT ( Gabcd ) ;

    rho    = 1 / ( Gab + Gcd ) ;
    Gab_a  = Gab * a ;
    Gab_b  = Gab * b ;
    Gcd_c  = Gcd * c ;
    Gcd_d  = Gcd * d ;
    PM     = Gab_b * Gab_b * ABAB + ACAC
               + 2 * Gab_b * BAAC ; 
    PMC    = -2 * Gab_b * ABCD ;
    QM     = Gcd_d * Gcd_d * CDCD + 2 * Gcd_d * CDAC ;
    QMC    = Gcd_d ;
    T      = rho * ( PM + QM + PMC * QMC ) ;
    Cv     = Const * den * Kab * Kcd ;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
    C_ab_1 [ 0 ] = - Gab_b * AB [ 0 ] ;
    C_ab_1 [ 1 ] = - Gab_b * AB [ 1 ] ;
    C_ab_1 [ 2 ] = - Gab_b * AB [ 2 ] ;

    C_ab_2 [ 0 ] = rho * ( Gab * Gab_b * AB [ 0 ] - Gab * Gcd_d * CD [ 0 ]
                         - Gab * AC [ 0 ] ) ;
    C_ab_2 [ 1 ] = rho * ( Gab * Gab_b * AB [ 1 ] - Gab * Gcd_d * CD [ 1 ] 
                         - Gab * AC [ 1 ] ) ;
    C_ab_2 [ 2 ] = rho * ( Gab * Gab_b * AB [ 2 ] - Gab * Gcd_d * CD [ 2 ]
                         - Gab * AC [ 2 ] ) ;

    C_ab_3 [ 0 ] =   HALF * Gab ;
    C_ab_4 [ 0 ] = - HALF * Gab * Gab * rho ;
    C_ab_5 [ 0 ] =   HALF * Gab * Gcd * rho ;

    C_cd_1 [ 0 ] = - Gcd_d * CD [ 0 ] ;
    C_cd_1 [ 1 ] = - Gcd_d * CD [ 1 ] ;
    C_cd_1 [ 2 ] = - Gcd_d * CD [ 2 ] ;

    C_cd_2 [ 0 ] = rho * ( Gcd * Gcd_d * CD [ 0 ] - Gcd * Gab_b * AB [ 0 ]
                         + Gcd * AC [ 0 ] ) ;
    C_cd_2 [ 1 ] = rho * ( Gcd * Gcd_d * CD [ 1 ] - Gcd * Gab_b * AB [ 1 ]
                         + Gcd * AC [ 1 ] ) ;
    C_cd_2 [ 2 ] = rho * ( Gcd * Gcd_d * CD [ 2 ] - Gcd * Gab_b * AB [ 2 ]
                         + Gcd * AC [ 2 ] ) ;

    C_cd_3 [ 0 ] =   HALF * Gcd ;
    C_cd_4 [ 0 ] = - HALF * Gcd * Gcd * rho ;
    C_cd_5 [ 0 ] =   HALF * Gab * Gcd * rho ;

    fmt( ssss, M_value, T, Cv ) ;

///fprintf ( stdout, "3: Kab: %12.4e\n", Kab ) ;
///fprintf ( stdout, "1: a, b, Gab, ABAB, Ov: %12.4e%12.4e%12.4e%12.4e%12.4e\n",
///                   a, b, Gab, ABAB, Ovchab ) ;
///fprintf ( stdout, "T, Cv, den, Kab, Kcd: %12.4e%12.4e%12.4e%12.4e%12.4e\n",
///                   T, Cv, den, Kab, Kcd ) ;
///fprintf ( stdout, "%12.4e%12.4e%12.4e%12.4e%12.4e\n",
///                   ssss[0], ssss[1], ssss[2], ssss[3], ssss[4] ) ;
/*
    if ( M_value == 2 ) {
        int i ;
        for ( i = 0 ; i <= M_value ; i++ ) {
           fprintf ( stdout, "ssss: %5d: %23.16e\n", i, ssss [ i ] ) ;
        }
    }
*/

////////////////////////////////////////////////////////
////// for debug ///////////////////////////////////////////
////////////////////////////////////////////////////////
///    for ( i = 0 ; i < 3 ; i++ ) {
///        D_AB [ i ] =   AB [ i ] ;
///        D_BA [ i ] = - AB [ i ] ;
///        D_CD [ i ] =   CD [ i ] ;
///        D_DC [ i ] = - CD [ i ] ;
///        D_AC [ i ] =   AC [ i ] ;
///        D_CA [ i ] = - AC [ i ] ;
///    }
///    //for debug end:
///    fprintf ( stdout, "ABAB, CDCD, ACAC, ABAC, CDAC, ABCD:\n\t%16.8e %16.8e %16.8e\n\t%16.8e %16.8e %16.8e\n", ABAB, CDCD, ACAC, -BAAC, CDAC, ABCD ) ;
///    fprintf ( stdout, "b, d, Gab, Gcd =%16.8e,%16.8e: %16.8e,%16.8e\n", b, d, Gab, Gcd ) ;
///    fprintf ( stdout, "den, CKab, Kcd =%16.8e,%16.8e,%16.8e\n",
///              den, Const*Kab*D_COEa*D_COEb, Kcd*D_COEc*D_COEd ) ;
///
///    D_AUX_0_011_000 =       Gab_b ;
///    D_AUX_1_012_000 = rho * Gab * Gab_b ;
///    D_AUX_1_001_011 = rho * Gab * Gcd_d ;
///    D_AUX_1_001_000 = rho * Gab ;
///    D_AUX_0_001_000 =       Gab ;
///    D_AUX_1_002_000 = rho * Gab * Gab ;
///
///    D_AUX_0_000_011 =       Gcd_d ;
///    D_AUX_1_000_012 = rho * Gcd * Gcd_d ;
///    D_AUX_1_011_001 = rho * Gcd * Gab_b ;
///    D_AUX_1_000_001 = rho * Gcd ;
///    D_AUX_0_000_001 =       Gcd ;
///    D_AUX_1_000_002 = rho * Gcd * Gcd ;
///
///    D_AUX_1_001_001 = rho * Gab * Gcd ;
///
///    fprintf ( stdout, "rho, PM, QM, PMc, QMc =\n\t%16.8e %16.8e %16.8e\n\t%16.8e %16.8e\n", rho, PM, QM, PMC, QMC ) ;
///    fprintf ( stdout, "C, T =%16.8e,  %16.8e\n", Cv*D_COEa*D_COEb*D_COEc*D_COEd, T ) ;
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
    return 0 ;
}
