// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _MATOPE_H_
#define _MATOPE_H_

extern void dcopy2(int, double* , double*);
extern void dscale(int n, double alpha, double* x);
extern void init_darray(int, double*, double);
extern void init_iarray(int, int*, int);
extern double dmaxabs(int n, double *s);
extern int imaxval(int n, int* s);
extern double dot_product(int, double*, double*);

extern void scale_diag(int, double*, double);

extern double max_diag_diff(int, double*, double);
extern int cholesky_dec(int, double*, double*);
extern void copy_diag_(int, double* , double*);
extern int eigen(int, double*, double*, double*, double*, double*);
extern int leq1_(int, double*, double*, int*);
extern int leq_(int, double*, double*, int*);
extern void dspmm(int, double*, double*, double*);
extern void dspsp(int, double*, double*, double*);
extern void sym_mul2_(int, double*, double*, double*, double*);
extern void daxpy2_(int, double, double*, double*);
extern void dspmv2_(int, double, double*, double*, double, double*);
extern void dbl_swap_elements_(int i, int j, double* x);
extern void int_swap_elements_(int i, int j, int* x);

// 行列の保存形式を変換する関数

extern void conv_OSP_OS(int, double*, double*);
extern void conv_OS_OSP(int, double*, double*);
extern void conv_SSP2_OS(int, double*, double*, int*);
extern void conv_SSP2_OSP(int, double*, double*, int*);
extern void conv_OS_SSP(int, double*, double*, int*);
extern void conv_OSP_SSP(int, double*, double*, int*);
extern void conv_SSP_OSP(int, double*, double*, int*);
extern void conv_inplace_SSP2_OSP(int, double*, double*, int*);
extern void conv_inplace_SSP_OSP(int, double*, double*, int*);
extern void conv_inplace_OSP_SSP(int, double*, double*, int*);

extern int diis_alloc( int numb_ao ) ;
extern int diis_update(double S[], const double D[], double F[]) ;

#endif
