#ifndef _INCLUDE_CEXPRSELECT_
#define _INCLUDE_CEXPRSELECT_

#define _SZ_PQRS_ 5
const char *filename = "Data_expr_c" ;

typedef struct {
    int    ij ;
    int    pqrs ;
    double coef ;
} anexpr_t ;

typedef struct {
    int addr ;
    int size ;
} addr_size_t ;

typedef struct {
    char       *filename ;
    int        norb ;
    int        ncsf ;
    char       **csf ;
    int        nexpr ;
    anexpr_t   *expr ;
    addr_size_t *as ;
} cexprselect_t ;

#endif
