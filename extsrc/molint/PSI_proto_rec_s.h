// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_SINT_
#define _INCLUDE_PSI_PROTO_SINT_

PSI_LI int1_s_rec ( PSI_LI inttype, PSI_LD *coef_ver, PSI_LD *w, PSI_LD **s_prim ) ;
PSI_LI int1_s_10_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_11_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_20_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_21_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_22_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_30_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_31_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_32_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_33_ver ( PSI_LD *coef, PSI_LD *w ) ;
PSI_LI int1_s_vercal_10_0 ( PSI_LD *coef,
                        PSI_LD *S_00_0, PSI_LD *DMY_1, PSI_LD *DMY_2,
                        PSI_LD *S_10 ) ;
PSI_LI int1_s_vercal_20_0 ( PSI_LD *coef,
                        PSI_LD *S_10_0, PSI_LD *S_00_1, PSI_LD *DMY_2,
                        PSI_LD *S_20 ) ;
PSI_LI int1_s_vercal_30_0 ( PSI_LD *coef,
                        PSI_LD *S_20_0, PSI_LD *S_10_1, PSI_LD *DMY_2,
                        PSI_LD *S_30 ) ;
PSI_LI int1_s_vercal_11_1 ( PSI_LD *coef,
                        PSI_LD *S_10_0, PSI_LD *DMY_1, PSI_LD *S_00_2,
                        PSI_LD *S_11 ) ;
PSI_LI int1_s_vercal_21_1 ( PSI_LD *coef,
                        PSI_LD *S_20_0, PSI_LD *DMY_1, PSI_LD *S_10_2,
                        PSI_LD *S_21 ) ;
PSI_LI int1_s_vercal_31_1 ( PSI_LD *coef,
                        PSI_LD *S_30_0, PSI_LD *DMY_1, PSI_LD *S_20_2,
                        PSI_LD *S_31 ) ;
PSI_LI int1_s_vercal_22_1 ( PSI_LD *coef,
                        PSI_LD *S_21_0, PSI_LD *S_20_1, PSI_LD *S_11_2,
                        PSI_LD *S_22 ) ;
PSI_LI int1_s_vercal_32_1 ( PSI_LD *coef,
                        PSI_LD *S_31_0, PSI_LD *S_30_1, PSI_LD *S_21_2,
                        PSI_LD *S_32 ) ;
PSI_LI int1_s_vercal_33_1 ( PSI_LD *coef,
                        PSI_LD *S_32_0, PSI_LD *S_31_1, PSI_LD *S_22_2,
                        PSI_LD *S_33 ) ;

#endif
