// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _LAPACK_
#define _LAPACK_

/* LAPACK driver routines */
// 一般化固有値問題(圧縮形式)
extern void dspgv_(int* ITYPE, char* JOBZ, char* UPLO, int* N,
	double* AP, double* BP, double* W, double* Z, int* LDZ,
	double* WORK, int* INFO);
// 対称固有値問題(圧縮形式)
extern void dspev_(char* JOBZ, char* UPLO, int* N, double* AP,
	double* W, double* Z, int* LDZ, double* WORK, int* INFO);
// 対称連立一次方程式(圧縮形式)
extern void dspsv_(char* UPLO, int* N, int* NRHS, double* AP,
	int* IPIV, double* B, int* LDB, int* INFO);
// 対称連立一次方程式(非圧縮形式)
extern void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV,
	double* B, int* LDB, int* INFO);

// standard eigenproblem for symmetric matrix
extern void dsyev_(char* JOBZ, char* UPLO, int* N, double* A,
	int* LDA, double* W, double* WORK, int* LWORK, int* INFO);

/* LAPACK routines needed to obtain eigen-values and vectores */
extern void dpptrf_(char* UPLO, int* N, double* AP, int* INFO);
extern void dspgst_(int* ITYPE, char* UPLO, int* N, double* AP,
	double* BP, int* INFO);
extern void dsptrd_(char* UPLO, int* N, double* AP, double* D,
	double* E, double* TAU, int* INFO);
extern void dopmtr_(char* SIDE, char* UPLO, char* TRANS, int* M,
	int* N, double* AP, double* TAU, double* C, int* LDC,
	double* WORK, int* INFO);

extern void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);
extern void dsygst_(int* ITYPE, char* UPLO, int* N, double* A,
	int* LDA, double* B, int* LDB, int* INFO);
extern void dsytrd_(char* UPLO, int* N, double* A, int* LDA, double* D,
	double* E, double* TAU, double* WORK, int* LWORK, int* INFO);
extern void dormtr_(char* SIDE, char* UPLO, char* TRANS, int* M,
	int* N, double* A, int* LDA, double* TAU, double* C, int* LDC,
	double* WORK, int* LWORK, int* INFO);

extern void dstebz_(char* RANGE, char* ORDER, int* N,
	double* VL, double* VU, int* IL, int* IU, double* ABSTOL,
	double* D, double* E, int* M, int* NSPLIT, double* W,
	int* IBLOCK, int* ISPLIT, double* WORK, int* IWORK, int* INFO);
extern void dsteqr_(char* COMPZ, int* N, double* D, double* E,
	double* Z, int* LDZ, double* WORK, int* INFO);
extern void dpteqr_(char* COMPZ, int* N, double* D, double* E,
	double* Z, int* LDZ, double* WORK, int* INFO);
extern void dstein_(int* N, double* D, double* E, int* M, double* W,
	int* IBLOCK, int* ISPLIT, double* Z, int* LDZ, double* WORK,
	int* IWORK, int* IFAIL, int* INFO);
extern void dsterf_(int *N, double* D, double* E, int* INFO);
extern void dstedc_(char* COMPZ, int* N, double* D, double* E,
	double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK,
	int* LIWORK, int* INFO);

extern void dopgtr_(char* UPLO, int* N, double* AP, double* TAU,
	double* Q, int* LDQ, double* WORK, int* INFO);
extern void dorgtr_(char* UPLO, int* N, double* A, int* LDA,
	double* TAU, double* WORK, int* LWORK, int* INFO);
extern double dlamch_(char*);

extern int ilaenv_(int* ISPEC, char* NAME, char* OPTS, int* N1,
	int* N2, int* N3, int* N4);

/* BLAS routines */
/* LEVEL 1 */
extern void dcopy_(int* N, double* X, int* INCX, double* Y, int* INCY);
extern int idamax_(int* N, double* X, int* INCX);
extern double ddot_(int* N, double* X, int* INCX, double* Y, int* INCY);
extern void dscal_(int* N, double* ALPHA, double* X, int* INCX);
extern void daxpy_(int* N, double* ALPHA, double* X, int* INCX,
	double* Y, int* INCY);

/* LEVEL 3 */
extern void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
	double* ALPHA, double* A, int* LDA, double* B, int* LDB,
	double* BETA, double* C, int* LDC);
extern void dtpsv_(char* UPLO, char* TRANS, char* DIAG, int* N,
	double* AP, double* X, int* INCX);
extern void dtrsv_(char* UPLO, char* TRANS, char* DIAG, int* N,
	double* A, int* LDA, double* X, int* INCX);

/* LEVEL 2 */
extern void dspmv_(char* UPLO, int* N, double* ALPHA, double* AP,
	double* X, int* INCX, double* BETA, double* Y, int* INCY);
extern void dgemv_(char* TRANS, int* M, int* N, double* ALPHA,
	double* A, int* LDA, double* X, int* INCX, double* BETA,
	double* Y, int* INCY);
#endif
