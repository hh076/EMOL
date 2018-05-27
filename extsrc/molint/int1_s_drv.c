#include <stdio.h>
#include <stdlib.h>

#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_fmt.h"
#include "PSI_proto_int1.h"
#include "PSI_proto_rec_s.h"

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
static PSI_LD *DINT_WORK ;
static PSI_LD TOL ;
static PSI_LI LEVEL_DBG_PRINT ;
static PSI_LI MAXCOMP   [ 2 ] ;
static PSI_LI ADDR_FACT [ 2 ] ;
static PSI_LI NSIZE_INT1_BYTYPE [ ] = { 0,
    1,   3,   9,   6,  18,
   36,  10,  30,  60, 100 } ;

///////////////////////////////////////////////////////////////////////////////
////// initialization & finalization //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSI_LI int1_s_init (
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

    DINT_WORK   = ( double * ) malloc ( _N_SINT_WORK_ * sizeof ( PSI_LD ) ) ;

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
    MAXCOMP [ 0 ] = max_val ;
    MAXCOMP [ 1 ] = 1       ;

    ///fmt_initialize ( _MAX_L_ ) ;
    
    return 0 ;
}

PSI_LI int1_s_fin ( )
{
    free ( SHEL_LQN ) ;
    free ( SHEL_TEM ) ;
    free ( SHEL_ATM ) ;
    free ( SHEL_ADD ) ;
    free ( ATOM_XYZ ) ;
    free ( PRIM_EXP ) ;
    free ( PRIM_COE ) ;
    free ( DINT_WORK ) ;

    ///fmt_finalize ( ) ;
    
    return 0 ;
}

///////////////////////////////////////////////////////////////////////////////
////// Contracted Shell INT1 calculation //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSI_LI int1_s_calc (
       PSI_LI Ish, PSI_LI Jsh,
       PSI_LI *_inttype, PSI_LI *_nsize_int, PSI_LD *_DINT_cont )
{
    PSI_LI m_value, inttype, nsize_int, a, b, ap, bp, i, inuc ;
    PSI_LI La, Lb, La_ex, Lb_ex, flg_La_gt_Lb, flg_a_and_b ;
    PSI_LD ABAB, c_ab, expa, expb, coea, coeb ;
    PSI_LD *Ca, *Cb, *Cc, *expo_a,  *expo_b, *coef_a, *coef_b ;
    PSI_LD *dint_prim ; 
    PSI_LD coefs_stv [ _N_COEFS1_ ] ;
    PSI_LD DINT_cont [ _N_STV_ ] ; 

    a = Ish ; b = Jsh ;

    La     = SHEL_LQN [ a ] ;
    Lb     = SHEL_LQN [ b ] ;
    Ca     = ATOM_XYZ + 3 * SHEL_ATM [ a ] ;
    Cb     = ATOM_XYZ + 3 * SHEL_ATM [ b ] ;
    expo_a = PRIM_EXP + SHEL_ADD [ a ] ;
    expo_b = PRIM_EXP + SHEL_ADD [ b ] ;
    coef_a = PRIM_COE + SHEL_ADD [ a ] ;
    coef_b = PRIM_COE + SHEL_ADD [ b ] ;
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////
    flg_a_and_b = ( a == b ) ;
    m_value     = La + Lb + 1 ;
    if ( La > Lb ) {
        flg_La_gt_Lb = 1 ;
        inttype      = ( La*(La+1)/2 + Lb+1 ) ;
        La_ex        = La ;
        Lb_ex        = Lb ;
    } else {
        flg_La_gt_Lb = 0 ; 
        inttype      = ( Lb*(Lb+1)/2 + La+1 ) ;
        La_ex        = Lb ;
        Lb_ex        = La ;
    }
    nsize_int   = NSIZE_INT1_BYTYPE [ inttype ] ;

    for ( i = 0 ; i < nsize_int ; i++ ) {
        DINT_cont [ i ] = 0 ;
    }
    for ( i = 0, ABAB = 0 ; i < 3 ; i++ ) {
        ABAB   += ( Ca [ i ] - Cb [ i ] ) * ( Ca [ i ] - Cb [ i ] ) ;
    }
    coefs_stv [ _ADDR_ABAB_ ] = ABAB ;

///fprintf ( stdout, "ish, jsh: %d, %d: %16.8e\n", Ish, Jsh, ABAB ) ;
    /// ===============================================
    /// CONTRACTION LOOP ==============================
    /// ===============================================
    for ( ap = 0 ; ap < SHEL_TEM [ a ] ; ap++ ) {
      expa = expo_a [ ap ] ;
      coea = coef_a [ ap ] ;
      for ( bp = 0 ; bp < SHEL_TEM [ b ] ; bp++ ) {
        expb = expo_b [ bp ] ;
        coeb = coef_b [ bp ] ;

        //------------------------------------
        // Preparing coefficients (PA <=> PB in the case: Lb > La)
        //------------------------------------
        int1_coefs_stv ( flg_La_gt_Lb,
                         expa, expb, Ca, Cb,
                         coefs_stv, MAXCOMP, ADDR_FACT ) ;

        //------------------------------------
        // Initial Integral Calculations ( IICs ) for S or T
        //------------------------------------
        int1_s_iic ( coefs_stv, DINT_WORK ) ;

        //------------------------------------
        // Recursive Calculations ( RCs ) for S and T
        //------------------------------------
        int1_s_rec ( inttype, coefs_stv, DINT_WORK, &dint_prim ) ; 

        //------------------------------------
        // Contraction
        //------------------------------------
        c_ab = coea * coeb ;
        for ( i = 0 ; i < nsize_int ; i++ ) {
          DINT_cont [ i ] += dint_prim [ i ] * c_ab ;
        }
      } // bp
    } // ap
///
///fprintf ( stdout, "end contraction:\n", Ish, Jsh ) ;
    int1_save ( inttype, La_ex, Lb_ex, ADDR_FACT,
                DINT_cont, _DINT_cont ) ;
    *_inttype   = inttype ;
    *_nsize_int = nsize_int ;

    return ( PSI_LI ) 0 ;
}

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//////////////////////////////////////////////////////
PSI_LI int1_s_print (
       PSI_LI Ash, PSI_LI Bsh,
       PSI_LI inttype, PSI_LI nsize_int,
       PSI_LD *Sint )
{
    PSI_LI i, a, b, la, lb ;
    PSI_LI ncomp_a, ncomp_b, offset_a, offset_b, iint ;

    la = SHEL_LQN [ Ash ] ;
    lb = SHEL_LQN [ Bsh ] ;

    fprintf ( stdout, "(%3d,%3d ) : (%2d,%2d ) : %d\n",
              Ash, Bsh, la, lb, inttype ) ;

    ncomp_a = ( ( la + 1 ) * ( la + 2 ) ) / 2 ;
    ncomp_b = ( ( lb + 1 ) * ( lb + 2 ) ) / 2 ;

    offset_a = 0 ;
    for ( i = 0 ; i < Ash ; i++ ) {
        offset_a += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }
    offset_b = 0 ;
    for ( i = 0 ; i < Bsh ; i++ ) {
        offset_b += ( ( SHEL_LQN [ i ] + 1 ) * ( SHEL_LQN [ i ] + 2 ) ) / 2 ;
    }

    iint = 0 ;
    for ( a = 0 ; a < ncomp_a ; a++ ) {
    for ( b = 0 ; b < ncomp_b ; b++ ) {
        fprintf ( stdout, "%5d%5d:%5d %24.16e\n",
                  offset_a + a, offset_b + b, iint, Sint [ iint ] ) ;
        iint++ ;
    }
    }
  
    return 0 ;
}
