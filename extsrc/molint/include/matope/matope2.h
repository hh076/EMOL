// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _FMO_MAT_H_
#define _FMO_MAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void fmo_mat2_unpack_matrix(int n, const double* P, double* S);
extern void fmo_mat2_pack_matrix(int n, const double* S, double* P);

extern int fmo_mat2_initialize(int maxn);
extern int fmo_mat2_perror(const char* s);
extern size_t fmo_mat2_temp_memory_size();

extern void fmo_mat2_dcopy(int n, double src[], double dest[]);
extern void fmo_mat2_daxpy(int n, double alpha, double x[], double y[]);
extern double fmo_mat2_ddot(int n, double x[], double y[]);
extern void fmo_mat2_dscale(int n, double alpha, double x[]);

extern void fmo_mat2_scale_diag(int n, double alpha, double X[]);
extern int fmo_mat2_isum(int n, int ix[]);
extern void fmo_mat2_transpose(int n, const double S[], double St[]);
extern void fmo_mat2_init_dvec(int n, const double alpha, double S[]);

extern int fmo_mat2_chodec(int n, double S[]);

extern int fmo_mat2_gsep2sep_square(int n, double A[], double U[]);
extern int fmo_mat2_gsep2sep_packed(int n, double AP[], double UP[]);
extern int fmo_mat2_gsep2sep_packed2(int n, double AP[], double U[]);

extern int fmo_mat2_gsep2sep2_square(int n, double A[], double U[]);
extern int fmo_mat2_gsep2sep2_packed2(int n, double AP[], double U[]);

extern int fmo_mat2_leq_square(int n, int m, double A[], double B[]);
extern int fmo_mat2_sep_square(int n, double A[], double e[]);
extern int fmo_mat2_solv_leq_u(int n, int m, double A[], double B[]);
extern int fmo_mat2_gsep_square(int n, double A[], double U[], double e[]);

extern int fmo_mat2_sym_mul(int n, double S[], double D[], double SDS[]);
extern int fmo_mat2_sym_update(int n, const double A[], const double B[],
	const double C[], double alpha, double E[]);

extern int fmo_mat2_asym_update(int n, const double FP[],
	const double DP[], const double U[], double E[]);

extern int fmo_mat2_dspsp(int n, const double AP[], const double BP[],
	double C[]);
#endif
