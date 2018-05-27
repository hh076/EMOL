// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_INT2_
#define _INCLUDE_PSI_PROTO_INT2_

#include "PSI_stan.h"
#include "PSI_proto_rec.h"

///////////////////////////////////
///////////////////////////////////
///////// DRIVER //////////////////
///////////////////////////////////
///////////////////////////////////
PSI_LI eri0_init (
       PSI_LI _NUMB_ATOM, PSI_LI _NUMB_SHELL, PSI_LI _NUMB_PRIM,
       PSI_LI *_SHEL_LQN, PSI_LI *_SHEL_TEM, PSI_LI *_SHEL_ATM,
       PSI_LI *_SHEL_ADD,
       PSI_LD *_ATOM_CHARG, PSI_LD *_ATOM_XYZ,
       PSI_LD *_PRIM_EXP,   PSI_LD *_PRIM_COE,
       PSI_LD _TOL,         PSI_LI _LEVEL_DBG_PRINT ) ;

PSI_LI eri0_calc (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI *_inttype, PSI_LI *_nsize_eris,
       PSI_LD *_eris ) ;

PSI_LI eri0_print (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI inttype, PSI_LI nsize_eris,
       PSI_LD *eris ) ;

///////////////////////////////////
///////// atomic //////////////////
///////////////////////////////////

PSI_LD eri0_atomic (
       PSI_LI i, PSI_LI j, PSI_LI k, PSI_LI l,
       PSI_LI inttype,
       PSI_LI flg_ex_ij, PSI_LI flg_ex_kl,
       PSI_LI flg_ex_ijkl,
       PSI_LD* Ci, PSI_LD* Cj, PSI_LD* Ck, PSI_LD* Cl,
       PSI_LD* Cij, PSI_LD* Ckl,
       PSI_LI* a, PSI_LI* b, PSI_LI* c, PSI_LI* d,
       PSI_LI *max_comp, PSI_LI *addr_fact,
       PSI_LD *abcdac ) ;
PSI_LI eri0_exchange_to_atomic_order (
       PSI_LI Ish, PSI_LI Jsh, PSI_LI Ksh, PSI_LI Lsh,
       PSI_LI Li, PSI_LI Lj, PSI_LI Lk, PSI_LI Ll,
       PSI_LI *a, PSI_LI *b, PSI_LI *c, PSI_LI *d ) ;

///////////////////////////////////
///////// save ////////////////////
///////////////////////////////////
PSI_LI eri0_save ( PSI_LI inttype, PSI_LI la, PSI_LI lb, PSI_LI lc, PSI_LI ld,
                   PSI_LI *addr_fact, PSI_LD *eri_from, PSI_LD *eri_to ) ;

///////////////////////////////////
///////////////////////////////////
///////////////////////////////////
///
#endif
