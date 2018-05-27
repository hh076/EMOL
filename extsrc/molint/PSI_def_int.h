/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      RHF Part of OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@c.csce.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_DEF_INT_
#define _INCLUDE_PSI_DEF_INT_

/// constant
#define BOHR_RADIUS    0.52917724924e0
#define HALF           0.5e0
#define PI             3.14159265358979323846264338327950288e0
//
#define _FLG_PACKED_   1     // 0 => integrals are saved as GAMESS format,
                             // 1 =>                     as packed format

#define _FOCK_THRZERO_ 1.0e-30
//

/// rec allocate parameter
#define _MAX_L_        3     // (limitation: max(L) = 3(f) )
#define _MAX_COEFREC_  (18 + 18)
//
#define _N_SS_         (2 * (_MAX_L_) + 1)
#define _N_SSSS_       (4 * (_MAX_L_) + 1)
#define _N_STV_        100   // _MAX_L_ = 1,2,3,4 =>  16,  36,  100,  225
#define _N_ERI_        10000 // _MAX_L_ = 1,2,3,4 => 256,1296,10000,50625 
#define _N_INTER_ERI_  84149 // _MAX_L_ = 1,2,3,4 => 308,7958,84149
#define _N_ERI_WORK_   ((_N_INTER_ERI_) + _N_ERI_)

#define _N_INTER_SINT_ 10000 // _MAX_L_ = 1,2,3,4 => 308,7958,84149
#define _N_SINT_WORK_  ((_N_INTER_SINT_) + _N_STV_)
#define _N_INTER_TINT_ 10000 // _MAX_L_ = 1,2,3,4 => 308,7958,84149
#define _N_TINT_WORK_  ((_N_INTER_TINT_) + _N_STV_)
#define _N_INTER_VINT_ 10000 // _MAX_L_ = 1,2,3,4 => 308,7958,84149
#define _N_VINT_WORK_  ((_N_INTER_VINT_) + _N_STV_)

/// one electron: coefficients memory map
#define _ADDR_PA_         0
#define _ADDR_PB_     ((_ADDR_PA_)    + 3)
#define _ADDR_PC_     ((_ADDR_PB_)    + 3)
#define _ADDR_P_      ((_ADDR_PC_)    + 3)
#define _ADDR_ABAB_   ((_ADDR_P_)     + 3)
#define _ADDR_PCPC_   ((_ADDR_ABAB_)  + 1)
#define _ADDR_G_2Z_   ((_ADDR_PCPC_)  + 1)
#define _ADDR_ZETA_   ((_ADDR_G_2Z_)  + 1)
#define _ADDR_2ZETA_  ((_ADDR_ZETA_)  + 1)
#define _ADDR_XI_     ((_ADDR_2ZETA_) + 1)
#define _ADDR_2XI_    ((_ADDR_XI_)    + 1)
#define _ADDR_ZA_     ((_ADDR_2XI_)   + 1)
#define _ADDR_ZB_     ((_ADDR_ZA_)    + 1)
#define _ADDR_XI_ZA_  ((_ADDR_ZB_)    + 1)
#define _ADDR_XI_ZB_  ((_ADDR_XI_ZA_) + 1)
#define _ADDR_AB_     ((_ADDR_XI_ZB_) + 1)
#define _ADDR_BA_     ((_ADDR_AB_)    + 3)
#define _ADDR_AC_     ((_ADDR_BA_)    + 3)
#define _ADDR_CA_     ((_ADDR_AC_)    + 3)
#define _ADDR_BC_     ((_ADDR_CA_)    + 3)
#define _ADDR_CB_     ((_ADDR_BA_)    + 3)

#define _N_COEFS1_    ((_ADDR_CB_)    + 3)

///#define _DEBUG_INT2SAVE_   ///active for debugging int2_save.c

#endif ///_INCLUDE_PSI_DEF_INT_
