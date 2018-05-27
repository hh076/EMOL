// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_HGMAT_
#define _INCLUDE_PSI_PROTO_HGMAT_

#include "PSI_stan.h"

PSI_LI hgmat_init (
       PSI_LI _NSHELL,      PSI_LI _NAO,
       PSI_LI _NSIZE_FOCK,  PSI_LI _NCSPLQN,   PSI_LI _NCSP,
       PSI_LI *_SHEL_LQN,   PSI_LI *_SHEL_ATM, PSI_LI *_SHEL_INI,
       PSI_LI *_CSP_ADD,    PSI_LI *_CSP_TEM,  PSI_LI *_CSP_AB,
       PSI_LD _THR_SCHWARZ, PSI_LD _THR_ERI,   PSI_LI _LEVEL_PRINT ) ;
PSI_LI gen_hmat ( PSI_LD *Smat, PSI_LD *Tmat, PSI_LD *Hmat ) ;
PSI_LI gen_gmat ( PSI_LD *D,    PSI_LD *G ) ;

PSI_LI copy_to_sthmat (
    PSI_LI ISH,     PSI_LI JSH,
    PSI_LI *INIS,   PSI_LI *LQNS,
    PSI_LD *S_cont, PSI_LD *T_cont, PSI_LD *V_cont,
    PSI_LD *Smat, PSI_LD *Tmat, PSI_LD *Hmat ) ;

PSI_LI eri0_partial_fock (
        PSI_LI _ish, PSI_LI _jsh, PSI_LI _ksh, PSI_LI _lsh,
        PSI_LI *inis, PSI_LI *lqns, PSI_LD *eris, PSI_LD *D, PSI_LD *G ) ;
///
#endif
