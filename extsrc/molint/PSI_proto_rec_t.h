// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_TINT_
#define _INCLUDE_PSI_PROTO_TINT_

PSI_LI int1_t_rec ( PSI_LI inttype, PSI_LD *coef_ver, PSI_LD *w, PSI_LD **t_prim ) ;
PSI_LI int1_t_10_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_11_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_20_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_21_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_22_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_30_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_31_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_32_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_33_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_t_vercal_10_0 ( PSI_LD *coef,
                        PSI_LD *T_00_0, PSI_LD *DMY_1, PSI_LD *DMY_2,
                        PSI_LD *S_10_3, PSI_LD *DMY_4, PSI_LD *T_10 ) ;
PSI_LI int1_t_vercal_20_0 ( PSI_LD *coef,
                        PSI_LD *T_10_0, PSI_LD *T_00_1, PSI_LD *DMY_2,
                        PSI_LD *S_20_3, PSI_LD *S_00_4, PSI_LD *T_20 ) ;
PSI_LI int1_t_vercal_30_0 ( PSI_LD *coef,
                        PSI_LD *T_20_0, PSI_LD *T_10_1, PSI_LD *DMY_2,
                        PSI_LD *S_30_3, PSI_LD *S_10_4, PSI_LD *T_30 ) ;
PSI_LI int1_t_vercal_11_1 ( PSI_LD *coef,
                        PSI_LD *T_10_0, PSI_LD *DMY_1, PSI_LD *T_00_2,
                        PSI_LD *S_11_3, PSI_LD *DMY_4, PSI_LD *T_11 ) ;
PSI_LI int1_t_vercal_21_1 ( PSI_LD *coef,
                        PSI_LD *T_20_0, PSI_LD *DMY_1, PSI_LD *T_10_2,
                        PSI_LD *S_21_3, PSI_LD *DMY_4, PSI_LD *T_21 ) ;
PSI_LI int1_t_vercal_31_1 ( PSI_LD *coef,
                        PSI_LD *T_30_0, PSI_LD *DMY_1, PSI_LD *T_20_2,
                        PSI_LD *S_31_3, PSI_LD *DMY_4, PSI_LD *T_31 ) ;
PSI_LI int1_t_vercal_22_1 ( PSI_LD *coef,
                        PSI_LD *T_21_0, PSI_LD *T_20_1, PSI_LD *T_11_2,
                        PSI_LD *S_22_3, PSI_LD *S_20_4, PSI_LD *T_22 ) ;
PSI_LI int1_t_vercal_32_1 ( PSI_LD *coef,
                        PSI_LD *T_31_0, PSI_LD *T_30_1, PSI_LD *T_21_2,
                        PSI_LD *S_32_3, PSI_LD *S_30_4, PSI_LD *T_32 ) ;
PSI_LI int1_t_vercal_33_1 ( PSI_LD *coef,
                        PSI_LD *T_32_0, PSI_LD *T_31_1, PSI_LD *T_22_2,
                        PSI_LD *S_33_3, PSI_LD *S_31_4, PSI_LD *T_33 ) ;

#endif
