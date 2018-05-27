// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_VINT_
#define _INCLUDE_PSI_PROTO_VINT_

PSI_LI int1_v_rec ( PSI_LI inttype, PSI_LD *coef_ver, PSI_LD *w, PSI_LD **v_prim ) ;
PSI_LI int1_v_10_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_11_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_20_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_21_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_22_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_30_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_31_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_32_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_33_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_v_vercal_10_0 ( PSI_LD *coef,
                        PSI_LD *V_00_0, PSI_LD *V_00_1, PSI_LD *DMY_2,
                        PSI_LD *DMY_3, PSI_LD *DMY_4, PSI_LD *DMY_5, PSI_LD *V_10 ) ;
PSI_LI int1_v_vercal_20_0 ( PSI_LD *coef,
                        PSI_LD *V_10_0, PSI_LD *V_10_1, PSI_LD *V_00_2,
                        PSI_LD *V_00_3, PSI_LD *DMY_4, PSI_LD *DMY_5, PSI_LD *V_20 ) ;
PSI_LI int1_v_vercal_30_0 ( PSI_LD *coef,
                        PSI_LD *V_20_0, PSI_LD *V_20_1, PSI_LD *V_10_2,
                        PSI_LD *V_10_3, PSI_LD *DMY_4, PSI_LD *DMY_5, PSI_LD *V_30 ) ;
PSI_LI int1_v_vercal_11_1 ( PSI_LD *coef,
                        PSI_LD *V_10_0, PSI_LD *V_10_1, PSI_LD *DMY_2,
                        PSI_LD *DMY_3, PSI_LD *V_00_4, PSI_LD *V_00_5, PSI_LD *V_11 ) ;
PSI_LI int1_v_vercal_21_1 ( PSI_LD *coef,
                        PSI_LD *V_20_0, PSI_LD *V_20_1, PSI_LD *DMY_2,
                        PSI_LD *DMY_3, PSI_LD *V_10_4, PSI_LD *V_10_5, PSI_LD *V_21 ) ;
PSI_LI int1_v_vercal_31_1 ( PSI_LD *coef,
                        PSI_LD *V_30_0, PSI_LD *V_30_1, PSI_LD *DMY_2,
                        PSI_LD *DMY_3, PSI_LD *V_20_4, PSI_LD *V_20_5, PSI_LD *V_31 ) ;
PSI_LI int1_v_vercal_22_1 ( PSI_LD *coef,
                        PSI_LD *V_21_0, PSI_LD *V_21_1, PSI_LD *V_20_2,
                        PSI_LD *V_20_3, PSI_LD *V_11_4, PSI_LD *V_11_5, PSI_LD *V_22 ) ;
PSI_LI int1_v_vercal_32_1 ( PSI_LD *coef,
                        PSI_LD *V_31_0, PSI_LD *V_31_1, PSI_LD *V_30_2,
                        PSI_LD *V_30_3, PSI_LD *V_21_4, PSI_LD *V_21_5, PSI_LD *V_32 ) ;
PSI_LI int1_v_vercal_33_1 ( PSI_LD *coef,
                        PSI_LD *V_32_0, PSI_LD *V_32_1, PSI_LD *V_31_2,
                        PSI_LD *V_31_3, PSI_LD *V_22_4, PSI_LD *V_22_5, PSI_LD *V_33 ) ;

#endif
