#include <stdio.h>
#include <stdlib.h>

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
    int        norb ;
    int        ncsf ;
    char       **csf ;
    int        nexpr ;
    anexpr_t   *expr ;
    addr_size_t *as ;
} expr_t ;

int **tri ;

int get_expr_csfpair( expr_t expr, int ij, int *nexpr, int *pqrs, double *coef )
{ 
    int *ids, iexpr ;
    int n1ints = expr.norb * ( expr.norb + 1 ) / 2 ;
    int addr   = expr.as[ ij ].addr ;
    *nexpr     = expr.as[ ij ].size ;
    if ( *nexpr > 0 ) {
        for ( iexpr = 0 ; iexpr < *nexpr ; iexpr++ ) {
            ids           = pqrs + _SZ_PQRS_ * iexpr ;
            ids [     4 ] = expr.expr[ addr + iexpr ].pqrs ;
            coef[ iexpr ] = expr.expr[ addr + iexpr ].coef ;
            if ( ids[ 4 ] > n1ints ) {
                ids[ 0 ] = ids[ 1 ] = ids[ 2 ] = ids[ 3 ] = 0 ;
            ///get_pqrs( pqrs[ i ].id, ids + 0, ids + 1, ids + 2, ids + 3 ) ;
            } else {
                ids[ 0 ] = ids[ 1 ] = 0 ;
                ids[ 2 ] = ids[ 3 ] = -1 ;
            }
        }
    }
    return 0 ;
}

int main()
{
    int max_tri = 1024 ;
    int i, j, ij, n, ival, max_csfpair, prev, iexpr ;
    char buf[ BUFSIZ ] ;
    expr_t expr ;
    FILE *fp = fopen( filename, "r" ) ;
///
    int    nexpr_curr ;
    int    *pqrs_curr = ( int    * ) malloc( max_tri * _SZ_PQRS_ * sizeof( int    ) ) ;
    double *coef_curr = ( double * ) malloc( max_tri *     sizeof( double ) ) ;
///    tri       = ( int ** ) malloc( ( max_tri * ( max_tri + 1 ) ) / 2 * sizeof( int * ) ) ;
///    tri       = ( int ** ) malloc( ( max_tri * ( max_tri + 1 ) ) / 2 * sizeof( int * ) ) ;
///    for ( i = 0 ; i < ( max_tri * ( max_tri + 1 ) ) / 2 ; i++ ) {
///        tri[ i ] = ( int * ) malloc( 2 * sizeof( int ) ) ;
///    }
///    for ( n = 0, i = 0 ; i < max_tri ; i++ ) {
///        for ( j = 0 ; j <= i ; j++, n++ ) {
///            tri[ n ][ 0 ] = i ;
///            tri[ n ][ 1 ] = j ;
///        }
///    }
///
    fgets( buf, sizeof( buf ), fp ) ;
    sscanf( buf, "%d%d%d", &(expr.norb), &(expr.ncsf), &(expr.nexpr) ) ;
///
    expr.csf = ( char ** ) malloc( expr.ncsf * sizeof( char * ) ) ;
    for ( i = 0 ; i < expr.ncsf ; i++ ) {
        expr.csf[ i ] = ( char * ) malloc( expr.norb * sizeof( char ) ) ;
    }
    expr.expr = ( anexpr_t * ) malloc( expr.nexpr * sizeof( anexpr_t ) ) ;
    max_csfpair = ( expr.ncsf * ( expr.ncsf + 1 ) ) / 2 ;
    expr.as = ( addr_size_t * ) malloc( max_csfpair * sizeof( addr_size_t ) ) ;
///
    for ( i = 0 ; i < expr.ncsf ; i++ ) {
        fgets( buf, sizeof( buf ), fp ) ;
    }
    for ( i = 0 ; i < expr.nexpr ; i++ ) {
        fgets( buf, sizeof( buf ), fp ) ;
        sscanf( buf, "%d%d%lf", &(expr.expr[ i ].ij), &(expr.expr[ i ].pqrs), &(expr.expr[ i ].coef) ) ;
    }
    for ( i = 0 ; i < max_csfpair ; i++ ) {
        expr.as[ i ].addr = -1 ;
        expr.as[ i ].size =  0 ;
    }

    n = 0 ;
    ival = 0 ;
    prev = -1 ;
    for ( i = 0 ; i < expr.nexpr ; i++ ) {
        int ij_curr = expr.expr[ i ].ij ;
        if ( ival != ij_curr ) {
            expr.as[ ij_curr ].addr = i ;
            expr.as[ ij_curr ].size = i - prev ;
            ival = ij_curr ;
            prev = i ;
            n++ ;
        }
    }

///    for ( i = 0 ; i < expr.nexpr ; i++ ) {
///        fprintf( stdout, "%10d%10d%10d %23.16e\n", i, expr.expr[ i ].ij, expr.expr[ i ].pqrs, expr.expr[ i ].coef ) ;
///    }
    for ( ij = 0 ; ij < max_csfpair ; ij++ ) {
        nexpr_curr = 0 ;
        get_expr_csfpair( expr, ij, &nexpr_curr, pqrs_curr, coef_curr ) ;
        for ( iexpr = 0 ; iexpr < nexpr_curr ; iexpr++ ) {
            fprintf( stdout, "%10d%4d : (%10d :%4d%4d%4d%4d ), %23.16e\n",
                     ij, iexpr,
                     pqrs_curr[ _SZ_PQRS_ * iexpr + 4 ], pqrs_curr[ _SZ_PQRS_ * iexpr + 0 ], pqrs_curr[ _SZ_PQRS_ * iexpr + 1 ],
                                                         pqrs_curr[ _SZ_PQRS_ * iexpr + 2 ], pqrs_curr[ _SZ_PQRS_ * iexpr + 3 ],
                     coef_curr[ iexpr ] ) ;
        }
    }
/*
    for ( i = 0 ; i < max_csfpair ; i++ ) {
        int addr = expr.as[ i ].addr ;
        int size = expr.as[ i ].size ;
        ///fprintf( stdout, "%10d%10d%10d\n", i, expr.as[ i ].addr, expr.as[ i ].size ) ;
        
        if ( addr >= 0 ) {
            for ( iexpr = 0 ; iexpr < size ; iexpr++ ) {
                int    ij   = expr.expr[ addr + iexpr ].ij ;
                int    pqrs = expr.expr[ addr + iexpr ].pqrs ;
                double coef = expr.expr[ addr + iexpr ].coef ;
                fprintf( stdout, "%10d%10d%10d: %10d%10d %23.16e\n", i, addr, iexpr, ij, pqrs, coef ) ;
            }
        }
    }
*/
    free( pqrs_curr ) ; 
    free( coef_curr ) ; 
    free( expr.csf ) ; 
    free( expr.expr ) ; 
    free( expr.as ) ; 
    fclose( fp ) ;
    return 0 ;
}
