// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#include "PSI_stan.h"
#include "PSI_def_int.h"
#include "PSI_proto_int2.h"

static PSI_LI ERI_FORCE_EXCHANGE [] = {
   0, 0, 0, 0, 0,  0, 0, 0, 0, 3,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0, 0, 0, 0, 3,  0, 3, 3, 0, 0,
   0, 0, 0, 0, 3,  0, 0, 0, 0, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0, 3, 3, 3, 0,  0, 0, 0,
   0, 0, 0, 0, 0,  3, 3, 0, 0, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 3, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0, 0, 0, 0, 0,  0, 0, 0, 0, 0,
   0 } ;

static PSI_LI ASH = 0, BSH = 0, CSH = 0, DSH = 0 ;

PSI_LD eri0_atomic (
                PSI_LI i, PSI_LI j, PSI_LI k, PSI_LI l,
                PSI_LI inttype,
                PSI_LI flg_ex_ij, PSI_LI flg_ex_kl,
                PSI_LI flg_ex_ijkl,
                PSI_LD* Ci, PSI_LD* Cj, PSI_LD* Ck, PSI_LD* Cl,
                PSI_LD* Cij, PSI_LD* Ckl,
                PSI_LI* a, PSI_LI* b, PSI_LI* c, PSI_LI* d,
                PSI_LI *max_comp, PSI_LI *addr_fact,
                PSI_LD *coef_rec )
{
    PSI_LI force_ex ;
    PSI_LD *Ca, *Cb, *Cc, *Cd ;
    
    force_ex =  ERI_FORCE_EXCHANGE [ inttype ] ;
    if ( ( flg_ex_ijkl <= 0 && ( force_ex <= 0 ) ) ||
         ( flg_ex_ijkl >  0 && ( force_ex >  0 ) ) ) {
        if ( flg_ex_ij <= 0 ) {
            // (ij,**)
            *a              = i ;
            *b              = j ;
            coef_rec [  0 ] =   Cij [ 0 ] ;
            coef_rec [  1 ] =   Cij [ 1 ] ;
            coef_rec [  2 ] =   Cij [ 2 ] ;
            coef_rec [  3 ] = - Cij [ 0 ] ;
            coef_rec [  4 ] = - Cij [ 1 ] ;
            coef_rec [  5 ] = - Cij [ 2 ] ;
            Ca              =   Ci ;
            Cb              =   Cj ;
            addr_fact [ 0 ] = max_comp [ 0 ] ;
            addr_fact [ 1 ] = max_comp [ 1 ] ;
        } else {
            // (ji,**)
            *a              = j ;
            *b              = i ;
            coef_rec [  0 ] = - Cij [ 0 ] ;
            coef_rec [  1 ] = - Cij [ 1 ] ;
            coef_rec [  2 ] = - Cij [ 2 ] ;
            coef_rec [  3 ] =   Cij [ 0 ] ;
            coef_rec [  4 ] =   Cij [ 1 ] ;
            coef_rec [  5 ] =   Cij [ 2 ] ;
            Ca              =   Cj ;
            Cb              =   Ci ;
            addr_fact [ 0 ] = max_comp [ 1 ] ;
            addr_fact [ 1 ] = max_comp [ 0 ] ;
        }
        if ( flg_ex_kl <= 0 ) {
            // (**,kl)
            *c              = k ;
            *d              = l ;
            coef_rec [  6 ] =   Ckl [ 0 ] ;
            coef_rec [  7 ] =   Ckl [ 1 ] ;
            coef_rec [  8 ] =   Ckl [ 2 ] ;
            coef_rec [  9 ] = - Ckl [ 0 ] ;
            coef_rec [ 10 ] = - Ckl [ 1 ] ;
            coef_rec [ 11 ] = - Ckl [ 2 ] ;
            Cc              =   Ck ;
            Cd              =   Cl ;
            addr_fact [ 2 ] = max_comp [ 2 ] ;
            addr_fact [ 3 ] = max_comp [ 3 ] ;
        } else {
            // (**,lk)
            *c              = l ;
            *d              = k ;
            coef_rec [  6 ] = - Ckl [ 0 ] ;
            coef_rec [  7 ] = - Ckl [ 1 ] ;
            coef_rec [  8 ] = - Ckl [ 2 ] ;
            coef_rec [  9 ] =   Ckl [ 0 ] ;
            coef_rec [ 10 ] =   Ckl [ 1 ] ;
            coef_rec [ 11 ] =   Ckl [ 2 ] ;
            Cc              =   Cl ;
            Cd              =   Ck ;
            addr_fact [ 2 ] = max_comp [ 3 ] ;
            addr_fact [ 3 ] = max_comp [ 2 ] ;
        }
    } else {
        if ( flg_ex_kl <= 0 ) {
            // (kl,**)
            *a              = k ;
            *b              = l ;
            coef_rec [  0 ] =   Ckl [ 0 ] ;
            coef_rec [  1 ] =   Ckl [ 1 ] ;
            coef_rec [  2 ] =   Ckl [ 2 ] ;
            coef_rec [  3 ] = - Ckl [ 0 ] ;
            coef_rec [  4 ] = - Ckl [ 1 ] ;
            coef_rec [  5 ] = - Ckl [ 2 ] ;
            Ca              =   Ck ;
            Cb              =   Cl ;
            addr_fact [ 0 ] = max_comp [ 2 ] ;
            addr_fact [ 1 ] = max_comp [ 3 ] ;
        } else {
            // (lk,**)
            *a              = l ;
            *b              = k ;
            coef_rec [  0 ] = - Ckl [ 0 ] ;
            coef_rec [  1 ] = - Ckl [ 1 ] ;
            coef_rec [  2 ] = - Ckl [ 2 ] ;
            coef_rec [  3 ] =   Ckl [ 0 ] ;
            coef_rec [  4 ] =   Ckl [ 1 ] ;
            coef_rec [  5 ] =   Ckl [ 2 ] ;
            Ca              =   Cl ;
            Cb              =   Ck ;
            addr_fact [ 0 ] = max_comp [ 3 ] ;
            addr_fact [ 1 ] = max_comp [ 2 ] ;
        }
        if ( flg_ex_ij <= 0 ) {
            // (**,ij)
            *c              = i ;
            *d              = j ;
            coef_rec [  6 ] =   Cij [ 0 ] ;
            coef_rec [  7 ] =   Cij [ 1 ] ;
            coef_rec [  8 ] =   Cij [ 2 ] ;
            coef_rec [  9 ] = - Cij [ 0 ] ;
            coef_rec [ 10 ] = - Cij [ 1 ] ;
            coef_rec [ 11 ] = - Cij [ 2 ] ;
            Cc              =   Ci ;
            Cd              =   Cj ;
            addr_fact [ 2 ] = max_comp [ 0 ] ;
            addr_fact [ 3 ] = max_comp [ 1 ] ;
        } else {
            // (**,ji)
            *c              = j ;
            *d              = i ;
            coef_rec [  6 ] = - Cij [ 0 ] ;
            coef_rec [  7 ] = - Cij [ 1 ] ;
            coef_rec [  8 ] = - Cij [ 2 ] ;
            coef_rec [  9 ] =   Cij [ 0 ] ;
            coef_rec [ 10 ] =   Cij [ 1 ] ;
            coef_rec [ 11 ] =   Cij [ 2 ] ;
            Cc              =   Cj ;
            Cd              =   Ci ;
            addr_fact [ 2 ] = max_comp [ 1 ] ;
            addr_fact [ 3 ] = max_comp [ 0 ] ;
        }
    }
    coef_rec [ 12 ] = Ca [ 0 ] - Cc [ 0 ] ;
    coef_rec [ 13 ] = Ca [ 1 ] - Cc [ 1 ] ; 
    coef_rec [ 14 ] = Ca [ 2 ] - Cc [ 2 ] ; 
    coef_rec [ 15 ] = - coef_rec [ 12 ] ;
    coef_rec [ 16 ] = - coef_rec [ 13 ] ;
    coef_rec [ 17 ] = - coef_rec [ 14 ] ;

    ASH = *a ; BSH = *b ; CSH = *c ; DSH = *d ;

    return 0 ;
}

///////////////////////////////////////////////////////////////////////////////
////// Exchange to Atomic Orger ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
PSI_LI eri0_exchange_to_atomic_order (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI Li,  PSI_LI Lj,  PSI_LI Lk,  PSI_LI Ll,
       PSI_LI *a, PSI_LI *b, PSI_LI *c, PSI_LI *d )
{

    *a = ASH, *b = BSH, *c = CSH, *d = DSH ;
/*
    PSI_LI Lij, Lkl ;
    PSI_LI flg_ex_ij, flg_ex_kl, flg_ex_ijkl ;

    *a = Ish ;
    *b = Jsh ;
    *c = Ksh ;
    *d = Lsh ;

    if ( Li >= Lj ) {
        flg_ex_ij = 0 ;
        Lij = ( Li*(Li+1)/2 + Lj+1 ) ;
    } else {
        flg_ex_ij = 1 ;
        Lij = ( Lj*(Lj+1)/2 + Li+1 ) ;
    }
    if ( Lk >= Ll ) {
       flg_ex_kl = 0 ;
       Lkl = Lk*(Lk+1)/2 + Ll+1 ;
    } else {
       flg_ex_kl = 1 ;
       Lkl = Ll*(Ll+1)/2 + Lk+1 ;
    }
    if ( Lij >= Lkl ) {
        flg_ex_ijkl = 0 ;
    } else {
        flg_ex_ijkl = 1 ;
    }

    if ( flg_ex_ijkl <= 0 ) {
        if ( flg_ex_ij <= 0 ) {
            // (ij,**)
            *a              = Ish ;
            *b              = Jsh ;
        } else {
            // (ji,**)
            *a              = Ksh ;
            *b              = Lsh ;
        }
        if ( flg_ex_kl <= 0 ) {
            // (**,kl)
            *c              = Ksh ;
            *d              = Lsh ;
        } else {
            // (**,lk)
            *c              = Lsh ;
            *d              = Ksh ;
        }
    } else {
        if ( flg_ex_kl <= 0 ) {
            // (kl,**)
            *a              = Ksh ;
            *b              = Lsh ;
        } else {
            // (lk,**)
            *a              = Lsh ;
            *b              = Ksh ;
        }
        if ( flg_ex_ij <= 0 ) {
            // (**,ij)
            *c              = Ish ;
            *d              = Jsh ;
        } else {
            // (**,ji)
            *c              = Jsh ;
            *d              = Ish ;
        }
    }
*/
    return 0 ;
}
