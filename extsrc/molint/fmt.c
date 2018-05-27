// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef int fmt_int;
typedef double fmt_dbl;

static fmt_dbl FMT_inv2 = 1.0 / 2.0;
static fmt_dbl FMT_inv6 = 1.0 / 6.0;
static fmt_int FMT_fmt_m;
static fmt_int FMT_fmt_inv_d = 512;	// 2^9
static fmt_int FMT_fmt_max_m;
static fmt_dbl FMT_fmt_t;
static fmt_int FMT_fmt_n_step;
static fmt_dbl *FMT_fmt_table;
static fmt_dbl FMT_fmt_step_size;
static fmt_dbl FMT_fmt_inv_step_size;
static fmt_dbl FMT_pi_div2;
static fmt_int FMT_fmt_max_m1;

static fmt_int initialized = 0;
static fmt_int MAXLQN = -1;

static void Free(void* p);

/* ---------------------------------------------------------------
    不完全Γ関数（誤差関数）計算の初期化ルーチン
------------------------------------------------------------------ */
void fmt_initialize(fmt_int max_lqn) {
    fmt_dbl thr_zero = 1.0e-17;
    fmt_int i,j,m,nu;
    fmt_dbl eps,t,expt,t2,term,func;
    fmt_int it0;

    if (initialized == 1) return;
    if (max_lqn <= MAXLQN) return;

    Free(FMT_fmt_table);

    // 基本的な変数の定義
    FMT_fmt_m = 4 * max_lqn + 2;
    FMT_fmt_max_m = FMT_fmt_m + 3;
    FMT_fmt_t = 2 * FMT_fmt_m + 36;
    FMT_fmt_n_step = FMT_fmt_t * FMT_fmt_inv_d;
    FMT_fmt_step_size = FMT_fmt_t / FMT_fmt_n_step;
    FMT_fmt_inv_step_size= 1.0 / FMT_fmt_step_size;
    FMT_fmt_max_m1 = FMT_fmt_max_m + 1;
    FMT_pi_div2 = 2.0 * atan(1.0);

    // 誤差関数計算のテーブルのためのメモリ確保
    FMT_fmt_table = (fmt_dbl*)malloc(sizeof(fmt_dbl)*(FMT_fmt_max_m+1)*(FMT_fmt_n_step+1));

    // 計算開始
	// T=0の場合の計算 Fm(0) = 1/(2m+1)
    for (m=0; m<=FMT_fmt_max_m; m++) FMT_fmt_table[m] = 1.0 / (2*m+1);

	// T>0の場合の計算 Fm(T) (T>0)
    m=FMT_fmt_max_m;
    for (j=1, it0=FMT_fmt_max_m1; j<=FMT_fmt_n_step; j++, it0 += FMT_fmt_max_m1){
	t = FMT_fmt_step_size * j;
	expt = exp(-t);
	nu = 2*m+1;
	t2 = 2.0 * t;
	eps = (expt/t2) * thr_zero;
	term = 1.0 / nu;
	func = term;
	i = nu;
	while (1) {
	    i += 2;
	    term *= t2 / i;
	    func += term;
	    if (term <= eps) break;
	}
	FMT_fmt_table[it0 + m] = expt * func;
	for (i=m-1; i>=0; i--){
	    nu -= 2;
	    FMT_fmt_table[it0 + i] = (expt + t2*FMT_fmt_table[it0 + (i+1)]) / nu;
	}
    }
    MAXLQN = max_lqn;
    initialized = 1;
}

/* ---------------------------------------------------------------
    不完全Γ関数（誤差関数）計算ルーチン
------------------------------------------------------------------ */
void fmt(fmt_dbl *f, fmt_int m, fmt_dbl t, fmt_dbl coef){
    fmt_int i,ts;
    fmt_dbl d,t_inv,nu;
    fmt_int ts0;
    // main
    if (t <= (2*m+36)) {
	ts = 0.5 + t * FMT_fmt_inv_step_size;
	d = ts * FMT_fmt_step_size - t;
	ts0 = ts * FMT_fmt_max_m1;
	for (i=0; i<=m; i++){
	    f[i] = coef * (((FMT_fmt_table[ts0 + i + 3] * FMT_inv6   * d
			   + FMT_fmt_table[ts0 + i + 2] * FMT_inv2 ) * d
			   + FMT_fmt_table[ts0 + i + 1]            ) * d
			   + FMT_fmt_table[ts0 + i + 0]            );
	}
    } else {
	t_inv = FMT_inv2 / t;
	f[0] = coef * sqrt(FMT_pi_div2*t_inv);
	nu = 1.0;
	for (i=1; i<=m; i++){
	    f[i] = t_inv * nu * f[i-1];
	    nu += 2.0;
	}
    }
}


/* ---------------------------------------------------------------
    不完全Γ関数（誤差関数）計算の後処理ルーチン
------------------------------------------------------------------ */
void fmt_finalize(){
    if (initialized == 0) return;
    Free(FMT_fmt_table);
    initialized = 0;
    MAXLQN = -1;
}

static void Free(void* p) {
    if (p != NULL) {
	free(p);
	p = NULL;
    }
}
