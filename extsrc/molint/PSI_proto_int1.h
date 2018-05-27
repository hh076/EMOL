// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_INT1_
#define _INCLUDE_PSI_PROTO_INT1_

#include "PSI_stan.h"

////////////////////////////////////////////
///////// ONEINT, TWOINT INIT, FIN /////////
////////////////////////////////////////////
PSI_LI int1_s_init (
           PSI_LI _NUMB_ATOM, PSI_LI _NUMB_SHELL, PSI_LI _NUMB_PRIM,
           PSI_LI *_SHEL_LQN, PSI_LI *_SHEL_TEM, PSI_LI *_SHEL_ATM,
           PSI_LI *_SHEL_ADD,
           PSI_LD *_ATOM_CHARG, PSI_LD *_ATOM_XYZ,
           PSI_LD *_PRIM_EXP,   PSI_LD *_PRIM_COE,
           PSI_LD _TOL,         PSI_LI _LEVEL_DBG_PRINT ) ;
PSI_LI int1_t_init (
           PSI_LI _NUMB_ATOM, PSI_LI _NUMB_SHELL, PSI_LI _NUMB_PRIM,
           PSI_LI *_SHEL_LQN, PSI_LI *_SHEL_TEM, PSI_LI *_SHEL_ATM,
           PSI_LI *_SHEL_ADD,
           PSI_LD *_ATOM_CHARG, PSI_LD *_ATOM_XYZ,
           PSI_LD *_PRIM_EXP,   PSI_LD *_PRIM_COE,
           PSI_LD _TOL,         PSI_LI _LEVEL_DBG_PRINT ) ;
PSI_LI int1_v_init (
           PSI_LI _NUMB_ATOM, PSI_LI _NUMB_SHELL, PSI_LI _NUMB_PRIM,
           PSI_LI *_SHEL_LQN, PSI_LI *_SHEL_TEM, PSI_LI *_SHEL_ATM,
           PSI_LI *_SHEL_ADD,
           PSI_LD *_ATOM_CHARG, PSI_LD *_ATOM_XYZ,
           PSI_LD *_PRIM_EXP,   PSI_LD *_PRIM_COE,
           PSI_LD _TOL,         PSI_LI _LEVEL_DBG_PRINT ) ;

PSI_LI int1_s_fin ( ) ;
PSI_LI int1_t_fin ( ) ;
PSI_LI int1_v_fin ( ) ;

PSI_LI int1_s_print (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI inttype, PSI_LI nsize_int,
           PSI_LD *DINT_cont ) ;
PSI_LI int1_t_print (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI inttype, PSI_LI nsize_int,
           PSI_LD *DINT_cont ) ;
PSI_LI int1_v_print (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI inttype, PSI_LI nsize_int,
           PSI_LD *DINT_cont ) ;

PSI_LI int1_s_calc (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI *inttype, PSI_LI *nsize_int,
           PSI_LD *DINT_cont ) ;
PSI_LI int1_t_calc (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI *inttype, PSI_LI *nsize_int,
           PSI_LD *DINT_cont ) ;
PSI_LI int1_v_calc (
           PSI_LI Ish, PSI_LI Jsh,
           PSI_LI *inttype, PSI_LI *nsize_int,
           PSI_LD *DINT_cont ) ;


PSI_LI int1_save ( int inttype, int la, int lb, int *addr_fact,
                       double *dint_from, double *dint_to ) ;

///////////////////////////////////
///////// ONEINT iic, rec /////////
///////////////////////////////////
PSI_LI int1_coefs_st ( PSI_LI flg_La_gt_Lb,
                       PSI_LD expa, PSI_LD expb,
                       PSI_LD *Ca,  PSI_LD *Cb,  PSI_LD *coefs,
                       PSI_LI *MAXCOMP, PSI_LI *ADDR_FACT ) ;
PSI_LI int1_coefs_v  ( PSI_LI flg_La_gt_Lb,
                       PSI_LD expa, PSI_LD expb, PSI_LD *Ca,
                       PSI_LD *Cb,  PSI_LD *Cc,  PSI_LD *coefs,
                       PSI_LI *MAXCOMP, PSI_LI *ADDR_FACT ) ;

PSI_LI int1_s_iic ( PSI_LD *coefs_st,  PSI_LD *DINT_WORK ) ;
PSI_LI int1_s_rec ( PSI_LI inttype,    PSI_LD *coefs_st,
                      PSI_LD *DINT_WORK, PSI_LD **dint_prim ) ; 

PSI_LI int1_t_iic ( PSI_LD *coefs_st,  PSI_LD *DINT_WORK ) ;
PSI_LI int1_t_rec ( PSI_LI inttype,    PSI_LD *coefs_st,
                      PSI_LD *DINT_WORK, PSI_LD **dint_prim ) ; 

PSI_LI int1_v_iic ( PSI_LI m_value,    PSI_LD Z_atom,
                      PSI_LD *coefs_v,   PSI_LD *DINT_WORK ) ;
PSI_LI int1_v_rec ( PSI_LI inttype,    PSI_LD *coefs_st,
                      PSI_LD *DINT_WORK, PSI_LD **dint_prim ) ; 

///////////////////////////////////
///////////////////////////////////
///////////////////////////////////
///
#endif
