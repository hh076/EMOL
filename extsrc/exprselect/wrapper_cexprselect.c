#include <stdio.h>
#include <stdlib.h>
#include "ruby.h"
#include "cexprselect.h"

char *sep = " ," ;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
static int readfile( char *filename, cexprselect_t *p )
{
    int i, j, n, ival, prev, max_csfpair ;
    char buf[ BUFSIZ ], *buf_tok ;
    FILE *fp = fopen( filename, "r" ) ;
///
    fgets( buf, sizeof( buf ), fp ) ;
    sscanf( buf, "%d%d%d", &(p->norb), &(p->ncsf), &(p->nexpr) ) ;
///
    p->csf = ( char ** ) ruby_xmalloc( p->ncsf * sizeof( char * ) ) ;
    for ( i = 0 ; i < p->ncsf ; i++ ) {
        p->csf[ i ] = ( char * ) ruby_xmalloc( p->norb * sizeof( char ) ) ;
    }
    p->expr = ( anexpr_t * ) ruby_xmalloc( p->nexpr * sizeof( anexpr_t ) ) ;
    max_csfpair = ( p->ncsf * ( p->ncsf + 1 ) ) / 2 ;
    p->as = ( addr_size_t * ) ruby_xmalloc( max_csfpair * sizeof( addr_size_t ) ) ;
///
    for ( i = 0 ; i < p->ncsf ; i++ ) {
        fgets( buf, sizeof( buf ), fp ) ;
        buf_tok = strtok( buf, sep ) ;
        j = 0 ;
        p->csf[ i ][ j++ ] = buf_tok[ 0 ] ;
        for ( ; j < p->norb ; j++ ) {
            buf_tok = strtok( NULL, sep ) ;
            p->csf[ i ][ j ] = buf_tok[ 0 ] ;
        }
    }
    for ( i = 0 ; i < p->nexpr ; i++ ) {
        fgets( buf, sizeof( buf ), fp ) ;
        sscanf( buf, "%d%d%lf", &(p->expr[ i ].ij), &(p->expr[ i ].pqrs), &(p->expr[ i ].coef) ) ;
    }
    for ( i = 0 ; i < max_csfpair ; i++ ) {
        p->as[ i ].addr = -1 ;
        p->as[ i ].size =  0 ;
    }

    n = 0 ;
    ival = 0 ;
    prev = -1 ;
    for ( i = 0 ; i < p->nexpr ; i++ ) {
        int ij_curr = p->expr[ i ].ij ;
        if ( ival != ij_curr ) {
            p->as[ ij_curr ].addr = i ;
            p->as[ ij_curr ].size = i - prev ;
            ival = ij_curr ;
            prev = i ;
            n++ ;
        }
    }
    fclose( fp ) ; 
    return 0 ;
}

static int get_expr_csfpair( cexprselect_t *p, int ij, int *nexpr, int *pqrs, double *coef )
{ 
    int *ids, iexpr ;
    int n1ints = p->norb * ( p->norb + 1 ) / 2 ;
    int addr   = p->as[ ij ].addr ;
    *nexpr     = p->as[ ij ].size ;
    ///fprintf( stderr, "addr, *nexpr: %d, %d\n", addr, *nexpr ) ;
    if ( *nexpr > 0 ) {
        for ( iexpr = 0 ; iexpr < *nexpr ; iexpr++ ) {
            int offset = _SZ_PQRS_ * iexpr ;
            ///fprintf( stderr, "iexpr: %d, %16p, %16p, %16p, %16p\n",
            ///         iexpr, &(pqrs[ offset + 4 ]), &(coef[ iexpr ]), &(p->expr[ addr + iexpr ].pqrs), &(p->expr[ addr + iexpr ].coef) ) ;
            pqrs[ offset + 4 ] = p->expr[ addr + iexpr ].pqrs ;
            coef[ iexpr ]      = p->expr[ addr + iexpr ].coef ;
            ///fprintf( stderr, "pqrs, coef: %d, %e\n", pqrs[ offset + 4 ], coef[ iexpr ] ) ;
            if ( pqrs[ offset + 4 ] > n1ints ) {
                 pqrs[ offset + 0 ] = pqrs[ offset + 1 ] = pqrs[ offset + 2 ] = pqrs[ offset + 3 ] = 0 ;
            ///get_pqrs( pqrs[ i ].id, ids + 0, ids + 1, ids + 2, ids + 3 ) ;
            } else {
                pqrs[ offset + 0 ] = pqrs[ offset + 1 ] = 0 ;
                pqrs[ offset + 2 ] = pqrs[ offset + 3 ] = -1 ;
            }
        }
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void wrap_CExprselect_free( cexprselect_t *p )
{
    ruby_xfree( p->csf ) ; 
    ruby_xfree( p->expr ) ; 
    ruby_xfree( p->as ) ; 
	ruby_xfree( p ) ;
}

VALUE wrap_CExprselect_allocate( VALUE self )
{
	void *p = ruby_xmalloc( sizeof( cexprselect_t ) ) ;
	return Data_Wrap_Struct( self, NULL, wrap_CExprselect_free, p ) ;
}

VALUE wrap_CExprselect_initialize( VALUE self )
{
    return Qnil ;
}

VALUE wrap_CExprselect_readfile( VALUE self, VALUE _filename )
{
	cexprselect_t *p;
	Data_Get_Struct( self, cexprselect_t, p ) ;
    char *filename = StringValuePtr( _filename ) ;
///
    readfile( filename, p ) ;
///
    return Qnil ;
}

VALUE wrap_CExprselect_get_expr( VALUE self, VALUE _ij, VALUE _pqrs, VALUE _coef )
{
    int    i, ij, max_tri = 1024 ;
    int    nexpr, *pqrs ;
    double *coef ;
	cexprselect_t *p;
	Data_Get_Struct( self, cexprselect_t, p ) ;
    pqrs = ( int    * ) ( _SZ_PQRS_ * max_tri * sizeof( int    ) ) ;
    coef = ( double * ) (             max_tri * sizeof( double ) ) ;
    ij   = FIX2INT( _ij ) ;
    get_expr_csfpair( p, ij, &nexpr, pqrs, coef ) ;
    for ( i = 0 ; i < nexpr ; i++ ) {
        rb_ary_push( _pqrs, INT2NUM(      pqrs[ i ] ) ) ;
        rb_ary_push( _coef, rb_float_new( coef[ i ] ) ) ;
    }
    free( pqrs ) ;
    free( coef ) ;
    return Qnil ;
}

VALUE wrap_CExprselect_cshow( VALUE self )
{
    int    j, icsf, ij, iexpr, max_csfpair, max_tri = 1024 ;
    int    nexpr, *pqrs ;
    double *coef ;
	cexprselect_t *p;
	Data_Get_Struct( self, cexprselect_t, p ) ;
    FILE *fp = fopen( "Dump_e", "w" ) ;
///
    fprintf( fp, "%d %d %d\n", p->ncsf, p->norb, p->nexpr ) ;
    for ( icsf = 0 ; icsf < p->ncsf ; icsf++ ) {
        for ( j = 0 ; j < p->norb ; j++ ) {
            fprintf( fp, "%2c", p->csf[ icsf ][ j ] ) ;
        }
        fprintf( fp, "\n" ) ;
    }
///
    pqrs = ( int    * ) ruby_xmalloc( _SZ_PQRS_ * max_tri * sizeof( int    ) ) ;
    coef = ( double * ) ruby_xmalloc(             max_tri * sizeof( double ) ) ;
    max_csfpair = ( p->ncsf * ( p->ncsf + 1 ) ) / 2 ;
    fprintf( fp, "max_csfpair: %d\n", max_csfpair ) ;
    fflush( fp ) ;
    for ( ij = 0 ; ij < max_csfpair ; ij++ ) {
        get_expr_csfpair( p, ij, &nexpr, pqrs, coef ) ;
        for ( iexpr = 0 ; iexpr < nexpr ; iexpr++ ) {
            fprintf( fp, "%10d %10d %10d %23.16e\n", ij, iexpr, pqrs[ _SZ_PQRS_ * iexpr + 4 ], coef[ iexpr ] ) ;
            fflush( fp ) ;
        }
    }
    ruby_xfree( pqrs ) ;
    ruby_xfree( coef ) ;
///
    fclose( fp ) ;
    return Qnil ;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Init_CExprselect ( )
{
    VALUE c = rb_define_class( "CExprselect",     rb_cObject ) ;
///
    rb_define_alloc_func    ( c,                  wrap_CExprselect_allocate ) ;
    rb_define_private_method( c, "initialize",    wrap_CExprselect_initialize,   0 ) ;
    rb_define_method(         c, "readfile",      wrap_CExprselect_readfile,     1 ) ;
    rb_define_method(         c, "get_expr",      wrap_CExprselect_get_expr,     3 ) ;
    rb_define_method(         c, "cshow",         wrap_CExprselect_cshow,        0 ) ;
}
