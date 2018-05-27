#include <stdio.h>
#include <stdlib.h>

#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_fmt.h"
#include "PSI_proto_int2.h"

static PSI_LI NUMB_ATOM ;
static PSI_LI NUMB_SHELL ;
static PSI_LI NUMB_PRIM ;
static PSI_LI *SHEL_LQN ;
static PSI_LI *SHEL_TEM ;
static PSI_LI *SHEL_ATM ;
static PSI_LI *SHEL_ADD ;
static PSI_LD *ATOM_CHARG ;
static PSI_LD *ATOM_XYZ ;
static PSI_LD *PRIM_EXP ;
static PSI_LD *PRIM_COE ;
static PSI_LD *ERI_WORK ;
static PSI_LD TOL ;
static PSI_LI LEVEL_DBG_PRINT ;
static PSI_LI ERI_MAX_COMP  [ 4 ] ;
static PSI_LI ERI_ADDR_FACT [ 4 ] ;

static PSI_LD *ZETA ;
static PSI_LD *DKPS ;
static PSI_LD *XIZA ;
static PSI_LD *ZETA2 ;
static PSI_LD *ETA2 ;
static PSI_LD *ZE2 ;
static PSI_LD *RZ ;
static PSI_LD *RE ;
static PSI_LD *PA ;
static PSI_LD *QC ;
static PSI_LD *WP ;
static PSI_LD *WQ ;
static PSI_LD *PC ;

///////////////////////////////////////////////////////////////////////////////
////// initialization & finalization //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSI_LI eri0_init (
       PSI_LI _NUMB_ATOM, PSI_LI _NUMB_SHELL, PSI_LI _NUMB_PRIM,
       PSI_LI *_SHEL_LQN, PSI_LI *_SHEL_TEM, PSI_LI *_SHEL_ATM,
       PSI_LI *_SHEL_ADD,
       PSI_LD *_ATOM_CHARG, PSI_LD *_ATOM_XYZ,
       PSI_LD *_PRIM_EXP, PSI_LD *_PRIM_COE,
       PSI_LD _TOL, PSI_LI _LEVEL_DBG_PRINT )
{
    PSI_LI i, j, max_lqn, max_val ;
    NUMB_ATOM  = _NUMB_ATOM ;
    NUMB_SHELL = _NUMB_SHELL ;
    NUMB_PRIM  = _NUMB_PRIM ;
    TOL        = _TOL ;
    LEVEL_DBG_PRINT
               = _LEVEL_DBG_PRINT ;

    ///printf ( "%d, %d, %d\n", NUMB_ATOM, NUMB_SHELL, NUMB_PRIM ) ;

    SHEL_LQN   = ( int * ) malloc ( NUMB_SHELL * sizeof ( PSI_LI ) ) ;
    SHEL_TEM   = ( int * ) malloc ( NUMB_SHELL * sizeof ( PSI_LI ) ) ;
    SHEL_ATM   = ( int * ) malloc ( NUMB_SHELL * sizeof ( PSI_LI ) ) ;
    SHEL_ADD   = ( int * ) malloc ( NUMB_SHELL * sizeof ( PSI_LI ) ) ;

    ATOM_CHARG = ( double * ) malloc ( NUMB_ATOM * sizeof ( PSI_LD ) ) ;
    ATOM_XYZ   = ( double * ) malloc ( NUMB_ATOM * sizeof ( PSI_LD ) * 3 ) ;
    PRIM_EXP   = ( double * ) malloc ( NUMB_PRIM * sizeof ( PSI_LD ) ) ;
    PRIM_COE   = ( double * ) malloc ( NUMB_PRIM * sizeof ( PSI_LD ) ) ;

    ERI_WORK   = ( double * ) malloc ( _N_ERI_WORK_ * sizeof ( PSI_LD ) ) ;

    for ( i = 0 ; i < NUMB_SHELL ; i++ ) {
        SHEL_LQN [ i ] = _SHEL_LQN [ i ] ;
        SHEL_TEM [ i ] = _SHEL_TEM [ i ] ;
        SHEL_ATM [ i ] = _SHEL_ATM [ i ] ;
        SHEL_ADD [ i ] = _SHEL_ADD [ i ] ;
    }

    for ( i = 0 ; i < NUMB_ATOM ; i++ ) {
        ATOM_CHARG [ i ] = _ATOM_CHARG [ i ] ;
        for ( j = 0 ; j < 3 ; j++ ) {
            ATOM_XYZ [ i * 3 + j ] = _ATOM_XYZ [ i * 3 + j ] ;
        }
    }

    for ( i = 0 ; i < NUMB_PRIM ; i++ ) {
        PRIM_EXP [ i ] = _PRIM_EXP [ i ] ;
        PRIM_COE [ i ] = _PRIM_COE [ i ] ;
    }

    max_lqn = 0 ;
    for ( i = 0 ; i < NUMB_SHELL ; i++ ) {
        if ( SHEL_LQN [ i ] > max_lqn ) {
            max_lqn = SHEL_LQN [ i ] ;
        }
    }
    max_val            = ( ( max_lqn + 1 ) * ( max_lqn + 2 ) ) / 2 ;
    if ( max_lqn == 1 ) {
        max_val = 4 ;
    }
    ERI_MAX_COMP [ 0 ] = max_val * max_val * max_val ;
    ERI_MAX_COMP [ 1 ] = max_val * max_val ;
    ERI_MAX_COMP [ 2 ] = max_val ;
    ERI_MAX_COMP [ 3 ] = 1 ;

    fmt_initialize ( _MAX_L_ ) ;
    
    return 0 ;
}

PSI_LI eri0_fin ( )
{
    free ( SHEL_LQN ) ;
    free ( SHEL_TEM ) ;
    free ( SHEL_ATM ) ;
    free ( SHEL_ADD ) ;
    free ( ATOM_XYZ ) ;
    free ( PRIM_EXP ) ;
    free ( PRIM_COE ) ;
    free ( ERI_WORK ) ;

    fmt_finalize ( ) ;
    
    return 0 ;
}

/*
PSI_LI eri0_coefs( PSI_LI as, PSI_LI bs, PSI_LI cs, PSI_LI ds )
{
    int a, b, c, d, ijps, ijcs ;
    tem_a    = ncomp_lqn( SHEL_TEM[ as ] ) ;
    tem_b    = ncomp_lqn( SHEL_TEM[ bs ] ) ;
    tem_c    = ncomp_lqn( SHEL_TEM[ cs ] ) ;
    tem_d    = ncomp_lqn( SHEL_TEM[ ds ] ) ;
    nelem = tem_a * tem_b * tem_c * tem_d ;
///
    ZETA  = ( PSI_LD * )malloc( tem_a * tem_b                     * sizeof( PSI_LD ) ) ;
    DKPS  = ( PSI_LD * )malloc( tem_a * tem_b                     * sizeof( PSI_LD ) ) ;
    XIZA  = ( PSI_LD * )malloc( tem_a * tem_b                     * sizeof( PSI_LD ) ) ;
    ZETA2 = ( PSI_LD * )malloc( tem_a * tem_b                     * sizeof( PSI_LD ) ) ;
    ETA2  = ( PSI_LD * )malloc( tem_c * tem_d                     * sizeof( PSI_LD ) ) ;
    ZE2   = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d     * sizeof( PSI_LD ) ) ;
    RZ    = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d     * sizeof( PSI_LD ) ) ;
    RE    = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d     * sizeof( PSI_LD ) ) ;
    PA    = ( PSI_LD * )malloc( tem_a * tem_b * 3                 * sizeof( PSI_LD ) ) ;
    QC    = ( PSI_LD * )malloc( tem_a * tem_b * 3                 * sizeof( PSI_LD ) ) ;
    WP    = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d * 3 * sizeof( PSI_LD ) ) ;
    WQ    = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d * 3 * sizeof( PSI_LD ) ) ;
    PC    = ( PSI_LD * )malloc( tem_a * tem_b * tem_c * tem_d * 3 * sizeof( PSI_LD ) ) ;
///
    abps = 0 ;
    ap0 = SHEL_ADD[ as ] ;
    bp0 = SHEL_ADD[ bs ] ;
    for ( a = 0 ; a < tem_a ; a++ ) {
        ap = ap0 + a ;
        for ( b = 0 ; b < tem_b ; b++ ) {
            bp = bp0 + b ;
            ZETA[ abps ]  = 1/( PRIM_EXP[ ap ] + PRIM_EXP[ bp ] ) ;
            DKPS[ abps ]  = sqrt( 2.0e0 ) * pi54 / ZETA[ abps ]
                          * exp( - PRIM_EXP[ ap ] * PRIM_EXP[ bp ] / ( PRIM_EXP[ ap ] + PRIM_EXP[ bp ] ) * ABAB )
            XIZA[ abps ]  = PRIM_EXP[ bp ] / ZETA[ abps ]
            ZETA2[ abps ] = ZETA * 2 ;
            abps++ ;
        }
    }
    cdps = 0 ;
    cp0 = SHEL_ADD[ cs ] ;
    dp0 = SHEL_ADD[ ds ] ;
    for ( c = 0 ; c < tem_c ; c++ ) {
        cp = cp0 + c ;
        for ( d = 0 ; d < tem_d ; d++ ) {
            dp = dp0 + d ;
            ZETA2[ abps ] = 2 * 1 / ( PRIM_EXP[ cp ] + PRIM_EXP[ dp ] ) ;;
            ZETA[ abps ]  = 1/( PRIM_EXP[ ap ] + PRIM_EXP[ bp ] ) ;
            DKPS[ abps ]  = sqrt( 2.0e0 ) * pi54 / ZETA[ abps ]
                          * exp( - PRIM_EXP[ ap ] * PRIM_EXP[ bp ] / ( PRIM_EXP[ ap ] + PRIM_EXP[ bp ] ) * ABAB )
            XIZA[ abps ]  = PRIM_EXP[ bp ] / ZETA[ abps ]
            abps++ ;
        }
    }

    ips = 0 ;
    for ( a = 0 ; a < tem_a ; a++ ) {
        ap = ap0 + a ;
        for ( b = 0 ; b < tem_b ; b++ ) {
            bp = bp0 + b ;
            for ( c = 0 ; c < tem_c ; c++ ) {
                cp = cp0 + c ;
                for ( d = 0 ; d < tem_d ; d++ ) {
                    dp = dp0 + d ;
                    RHO[ ips ] = 1/( zeta[ abps ] * eta[ cdps ] ) ;
                    ips++ ;

}
*/
///////////////////////////////////////////////////////////////////////////////
////// Contracted Shell ERI calculation ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSI_LI eri0_calc (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI *_inttype, PSI_LI *_nsize_eris,
       PSI_LD *_eris )
{
    int Ash, Bsh, Csh, Dsh ;
    int m_value, ii ;
    int Li, Lj, Lk, Ll, Lij, Lkl ;
    int La, Lb, Lc, Ld ;
    int flg_ex_ij, flg_ex_kl, flg_ex_ijkl ;
    int ipa, ipb, ipc, ipd, ieris ;
    int inttype, nsize_eris ;
    int NSIZE_ERIS_BYTYPE [ ] = {   0,
            1,    3,    9,    9,   27,   81,    6,   18,   54,   36,
           18,   54,  162,  108,  324,   36,  108,  324,  216,  648,
         1296,   10,   30,   90,   60,   180, 360,  100,   30,   90,
          270,  180,  540, 1080,  300,   900,  60,  180,  540,  360,
         1080, 2160,  600, 1800, 3600,   100, 300,  900,  600, 1800,
         3600, 1000, 3000, 6000, 10000 } ;
    double EXPa, EXPb, EXPc, EXPd, COEa, COEb, COEc, COEd ;
    double *Ci, *Cj, *Ck, *Cl, Cij [ 3 ], Ckl [ 3 ] ;
    double *expo_a, *expo_b, *expo_c, *expo_d,
           *coef_a, *coef_b, *coef_c, *coef_d ;
    double coef_rec [ _MAX_COEFREC_ ] ;
    double *coef_ver, *coef_hor, *eris_prim ;
    double eris [ _N_ERI_ ] ;

    coef_ver = coef_rec +  0 ;
    coef_hor = coef_rec + 18 ;

    Li = SHEL_LQN [ Ish ] ; Lj = SHEL_LQN [ Jsh ] ;
    Lk = SHEL_LQN [ Ksh ] ; Ll = SHEL_LQN [ Lsh ] ;
    Ci = ATOM_XYZ + 3 * SHEL_ATM [ Ish ] ;
    Cj = ATOM_XYZ + 3 * SHEL_ATM [ Jsh ] ;
    Ck = ATOM_XYZ + 3 * SHEL_ATM [ Ksh ] ;
    Cl = ATOM_XYZ + 3 * SHEL_ATM [ Lsh ] ;

    for ( ii = 0 ; ii < 3 ; ii++ ) {
        Cij [ ii ] = Ci [ ii ] - Cj [ ii ] ;
        Ckl [ ii ] = Ck [ ii ] - Cl [ ii ] ;
    }

    if ( Li >= Lj ) {
        flg_ex_ij = 0 ;
        Lij = (Li*(Li+1))/2 + Lj+1 ;
    } else {
        flg_ex_ij = 1 ;
        Lij = (Lj*(Lj+1))/2 + Li+1 ;
    }
    if ( Lk >= Ll ) {
        flg_ex_kl = 0 ;
        Lkl = (Lk*(Lk+1))/2 + Ll+1 ;
    } else {
        flg_ex_kl = 1 ;
        Lkl = (Ll*(Ll+1))/2 + Lk+1 ;
    }
    if ( Lij >= Lkl ) {
        flg_ex_ijkl = 0 ;
        inttype = (Lij*(Lij-1))/2 + Lkl ;
    } else {
        flg_ex_ijkl = 1 ;
        inttype = (Lkl*(Lkl-1))/2 + Lij ;
    }

    /// ==============================================
    m_value = Li + Lj + Lk + Ll ;
    nsize_eris = NSIZE_ERIS_BYTYPE [ inttype ] ;

    /// =========================================================
    /// exchange shell IDs to canonical(?) integral order =======
    /// =========================================================
    if ( _FLG_PACKED_ ) {
        ERI_MAX_COMP [ 3 ] = 1 ;
        ERI_MAX_COMP [ 2 ] = ERI_MAX_COMP [ 3 ] * ( ( ( Ll + 1 ) * ( Ll + 2 ) ) / 2 ) ;
        ERI_MAX_COMP [ 1 ] = ERI_MAX_COMP [ 2 ] * ( ( ( Lk + 1 ) * ( Lk + 2 ) ) / 2 ) ;
        ERI_MAX_COMP [ 0 ] = ERI_MAX_COMP [ 1 ] * ( ( ( Lj + 1 ) * ( Lj + 2 ) ) / 2 ) ;
    }
    eri0_atomic ( Ish, Jsh, Ksh, Lsh, 
                  inttype,
                  flg_ex_ij, flg_ex_kl, flg_ex_ijkl,
                  Ci, Cj, Ck, Cl, Cij, Ckl,
                  &Ash, &Bsh, &Csh, &Dsh,
                  ERI_MAX_COMP, ERI_ADDR_FACT, coef_hor ) ;
/*
    if ( inttype == 9 ) {
        int i ;
        fprintf ( stdout, "%4d%4d%4d%4d -> %4d%4d%4d%4d\n",
                          Ish, Jsh, Ksh, Lsh, Ash, Bsh, Csh, Dsh ) ;
        fprintf ( stdout, "%4d%4d%4d%4d: %4d%4d:%4d -> %d, %d, %d\n",
                          Li,Lj,Lk,Ll, Lij,Lkl, inttype,
                          flg_ex_ij, flg_ex_kl, flg_ex_ijkl ) ;
        fprintf ( stdout, "coef_hor:\n" ) ;
        for ( i = 0 ; i < 18 ; i++ ) {
            fprintf ( stdout, "%4d %16.6e\n", i, coef_hor [ i ] ) ;
        }
    }
*/
    La = SHEL_LQN [ Ash ] ; Lb = SHEL_LQN [ Bsh ] ;
    Lc = SHEL_LQN [ Csh ] ; Ld = SHEL_LQN [ Dsh ] ;

    for ( ieris = 0 ; ieris < nsize_eris ; ieris++ ) {
        eris [ ieris ] = 0 ;
    }

    /// ===============================================
    /// CONTRACTION 4 FOLD LOOP...            =========
    /// ===============================================
    expo_a = PRIM_EXP + SHEL_ADD [ Ash ] ;
    expo_b = PRIM_EXP + SHEL_ADD [ Bsh ] ;
    expo_c = PRIM_EXP + SHEL_ADD [ Csh ] ;
    expo_d = PRIM_EXP + SHEL_ADD [ Dsh ] ;
    coef_a = PRIM_COE + SHEL_ADD [ Ash ] ;
    coef_b = PRIM_COE + SHEL_ADD [ Bsh ] ;
    coef_c = PRIM_COE + SHEL_ADD [ Csh ] ;
    coef_d = PRIM_COE + SHEL_ADD [ Dsh ] ;
   
    for ( ipb = 0 ; ipb < SHEL_TEM [ Bsh ] ; ipb++ ) {
     EXPb = expo_b [ ipb ] ;
     COEb = coef_b [ ipb ] ;
     for ( ipa = 0 ; ipa < SHEL_TEM [ Ash ] ; ipa++ ) {
      EXPa = expo_a [ ipa ] ;
      COEa = coef_a [ ipa ] ;
      for ( ipd = 0 ; ipd < SHEL_TEM [ Dsh ] ; ipd++ ) {
       EXPd = expo_d [ ipd ] ;
       COEd = coef_d [ ipd ] ;
       for ( ipc = 0 ; ipc < SHEL_TEM [ Csh ] ; ipc++ ) {
        EXPc = expo_c [ ipc ] ;
        COEc = coef_c [ ipc ] ;
        //------------------------------------
        // Initial Integral Calculation ( IIC )
        // and coefficients for RC
        eri0_ssss ( m_value, EXPa, EXPb, EXPc, EXPd, coef_hor,
                   coef_ver, ERI_WORK ) ;

        // Recursive Calculation ( RC )
        eri0_rec  ( inttype, coef_ver, coef_hor, ERI_WORK, &eris_prim ) ;

        for ( ieris = 0 ; ieris < nsize_eris ; ieris++ ) {
            eris [ ieris ] += eris_prim [ ieris ]
                            * COEa * COEb * COEc * COEd ;
        }
       } // end ipc //
      } // end ipd //
     } // end ipa //
    } // end ipb //

    eri0_save ( inttype, La, Lb, Lc, Ld, ERI_ADDR_FACT, eris, _eris ) ;
    *_inttype    = inttype ;
    *_nsize_eris = nsize_eris ;

    ///eri0_print ( Ish, Jsh, Ksh, Lsh, inttype, nsize_eris, _eris ) ;

    return 0 ;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
PSI_LI eri0_print (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI inttype, PSI_LI nsize_eris,
       PSI_LD *eris )
{
    PSI_LI Ash, Bsh, Csh, Dsh ;
    PSI_LI i, a, b, c, d, la, lb, lc, ld ;
    PSI_LI Li, Lj, Lk, Ll ;
    PSI_LI ncomp_a, ncomp_b, ncomp_c, ncomp_d, ieris ;
    PSI_LI offset_a, offset_b, offset_c, offset_d ;

    Li = SHEL_LQN [ Ish ] ;
    Lj = SHEL_LQN [ Jsh ] ;
    Lk = SHEL_LQN [ Ksh ] ;
    Ll = SHEL_LQN [ Lsh ] ;

    eri0_exchange_to_atomic_order ( Ish, Jsh, Ksh, Lsh, Li, Lj, Lk, Ll,
                                    &Ash, &Bsh, &Csh, &Dsh ) ;
    la = SHEL_LQN [ Ash ] ;
    lb = SHEL_LQN [ Bsh ] ;
    lc = SHEL_LQN [ Csh ] ;
    ld = SHEL_LQN [ Dsh ] ;

    fprintf ( stdout, "(%3d%3d,%3d%3d ) -> (%3d%3d,%3d%3d ) : (%2d%2d,%2d%2d ) : %d\n",
                      Ish, Jsh, Ksh, Lsh, Ash, Bsh, Csh, Dsh,
                      la, lb, lc, ld, inttype ) ;

    ncomp_a = ( ( la + 1 ) * ( la + 2 ) ) / 2 ;
    ncomp_b = ( ( lb + 1 ) * ( lb + 2 ) ) / 2 ;
    ncomp_c = ( ( lc + 1 ) * ( lc + 2 ) ) / 2 ;
    ncomp_d = ( ( ld + 1 ) * ( ld + 2 ) ) / 2 ;

    offset_a = 0 ;
    for ( i = 0 ; i < Ash ; i++ ) {
        offset_a += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }
    offset_b = 0 ;
    for ( i = 0 ; i < Bsh ; i++ ) {
        offset_b += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }
    offset_c = 0 ;
    for ( i = 0 ; i < Csh ; i++ ) {
        offset_c += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }
    offset_d = 0 ;
    for ( i = 0 ; i < Dsh ; i++ ) {
        offset_d += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }

    ieris = 0 ;
    for ( a = 0 ; a < ncomp_a ; a++ ) {
    for ( b = 0 ; b < ncomp_b ; b++ ) {
    for ( c = 0 ; c < ncomp_c ; c++ ) {
    for ( d = 0 ; d < ncomp_d ; d++ ) {
        fprintf ( stdout, "%5d%5d,%5d%5d:%5d %24.16e\n",
                            offset_a + a, offset_b + b, offset_c + c, offset_d + d,
                            ieris, eris [ ieris ] ) ;
        ieris++ ;
    }
    }
    }
    }
  
    return 0 ;
}
