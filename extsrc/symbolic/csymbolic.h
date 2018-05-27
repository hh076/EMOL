#ifndef __INCLUDE_CSYMBOLIC__
#define __INCLUDE_CSYMBOLIC__

#include "ruby.h"
#include "drt.h"
#include "slib_list.h"
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

#define _MAX_DEPTH_              4
#define _MAX_NSIZE_MULTILOOP_0_  (1)
#define _MAX_NSIZE_MULTILOOP_1_  (64)
#define _MAX_NSIZE_MULTILOOP_2_  (64*64/2)
#define _MAX_NSIZE_MULTILOOP_3_  (64*64*64/4)
#define _MAX_NSIZE_MULTILOOP_4_  (64*64*64*64/8)
#define _MAX_NSIZE_EXPR_         (1024*1024*128)

typedef struct {
    int    ij ;
    int    pqrs ;
    double coef ;
} csymbolic_anexpression_t ;

typedef struct {
    int    depth ;
    int    nmin ;
    int    nmax ;
    int    nloop ;
    int    *res ;
} multiloop_t ;

typedef struct {
    int    norb_symbol ;
    int    norb_actual ;
    int    ncsf_symbol ;
    int    ncsf_actual ;
    int    size_ij ;
    int    size_addr ;
    int    size_pqrs ;
    int    size_coef ;
    int    n_1el_symbol ;
    int    n_1el ;
    int    n_ext ;
    int    mo_val_range_last ;
    int    mo_ext_range_first ;
    int    mo_ext_range_last ;
///
    int    **tri ;
    int    **csf_symbol ;
    int    **csf_actual ;
    char   **csf_symbolc ;
    char   **csf_actualc ;
    int    *ij ;
    int    *addr ;
    int    *pqrs ;
    double *coef ;
///
    int    nexpr_actual ;
    csymbolic_anexpression_t *expr_actual ;
///
    multiloop_t *mos ;
///
    edml_drt_t drt ;
///
} csymbolic_t ;

extern void wrap_CSymbolic_free( csymbolic_t *ptr ) ;
extern VALUE wrap_CSymbolic_allocate( VALUE self ) ;
extern VALUE wrap_CSymbolic_initialize( VALUE self,
    VALUE _norb_symbol, VALUE _norb_actual, VALUE _ncsf_symbol, VALUE _ncsf_actual,
    VALUE _size_ij, VALUE _size_addr, VALUE _size_pqrs, VALUE _size_coef,
    VALUE _n_1el_symbol, VALUE _n_1el, VALUE _n_ext, 
    VALUE _mo_val_range_last, VALUE _mo_ext_range_first, VALUE _mo_ext_range_last ) ;
extern VALUE wrap_CSymbolic_set_csf_symbol( VALUE self, VALUE _csf_symbol ) ;
extern VALUE wrap_CSymbolic_set_csf_actual( VALUE self, VALUE _csf_actual ) ;
extern VALUE wrap_CSymbolic_set_ij( VALUE self, VALUE _ij ) ;
extern VALUE wrap_CSymbolic_set_addr( VALUE self, VALUE _addr ) ;
extern VALUE wrap_CSymbolic_set_pqrs( VALUE self, VALUE _pqrs ) ;
extern VALUE wrap_CSymbolic_set_coef( VALUE self, VALUE _coef ) ;
extern VALUE wrap_CSymbolic_set_drt( VALUE self, VALUE _nelec, VALUE _spin2, VALUE _n_core, VALUE _n_val, VALUE _n_ext ) ;
extern void Init_CSymbolic ( ) ;

extern int main_actual_expression( csymbolic_t *p ) ;
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
#endif /// __INCLUDE_CSYMBOLIC__

