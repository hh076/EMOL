#include <stdio.h>
#include <stdlib.h>
#include "ruby.h"
#include "csymbolic.h"

const int MAX_NLOOP[] = { _MAX_NSIZE_MULTILOOP_0_, _MAX_NSIZE_MULTILOOP_1_,
                          _MAX_NSIZE_MULTILOOP_2_, _MAX_NSIZE_MULTILOOP_3_,
                          _MAX_NSIZE_MULTILOOP_4_ } ;

static int gen_loop( multiloop_t *pmloop, int k, int *val )
{
    int i, kk, nin ;
    if ( k == pmloop->depth ) {
        if ( pmloop->nloop >= MAX_NLOOP[ pmloop->depth ] ) {
            fprintf( stderr, "Error: gen_loop: depth, *nloop >= max_nloop : %4d, %d >= %d\n", pmloop->depth, pmloop->nloop, MAX_NLOOP[ pmloop->depth ] ) ;
            exit( 1 ) ;
        }
        for ( i = 0 ; i < pmloop->depth ; i++ ) {
            pmloop->res[ _MAX_DEPTH_ * pmloop->nloop + i ] = val[ i ] ;
        }
        (pmloop->nloop)++ ;
    } else {
        if ( k == 0 ) {
            nin = pmloop->nmin ;
        } else {
            nin = val[ k - 1 ] + 1 ;
        }
        for ( kk = nin ; kk <= pmloop->nmax ; kk++ ) {
            val[ k ] = kk ;
            gen_loop( pmloop, k + 1, val ) ;
        }
    }
    return 1 ;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void wrap_CSymbolic_free( csymbolic_t *p )
{
    int i ;
    for ( i = 0 ; i < p->ncsf_symbol ; i++ ) {
        ruby_xfree( p->csf_symbol [ i ] ) ;
        ruby_xfree( p->csf_symbolc[ i ] ) ;
    }
    for ( i = 0 ; i < p->ncsf_actual ; i++ ) {
        ruby_xfree( p->csf_actual [ i ] ) ;
        ruby_xfree( p->csf_actualc[ i ] ) ;
    }
    for ( i = 0 ; i < ( p->ncsf_symbol * ( p->ncsf_symbol + 1 ) / 2 ) ; i++ ) {
        ruby_xfree( p->tri[ i ] ) ;
    }
    ruby_xfree( p->csf_symbol ) ;
    ruby_xfree( p->csf_actual ) ;
    ruby_xfree( p->csf_symbolc ) ;
    ruby_xfree( p->csf_actualc ) ;
    ruby_xfree( p->ij ) ;
    ruby_xfree( p->addr ) ;
    ruby_xfree( p->pqrs ) ;
    ruby_xfree( p->coef ) ;
    ruby_xfree( p->tri ) ;
    if ( p->expr_actual != ( csymbolic_anexpression_t * ) NULL ) {
        ruby_xfree( p->expr_actual ) ;
    }
    for ( i = 0 ; i <= _MAX_DEPTH_ ; i++ ) {
        ruby_xfree( p->mos[ i ].res ) ;
    }
    ruby_xfree( p->mos ) ;
	ruby_xfree( p ) ;
}

VALUE wrap_CSymbolic_allocate( VALUE self )
{
	void *p = ruby_xmalloc( sizeof( csymbolic_t ) ) ;
	return Data_Wrap_Struct( self, NULL, wrap_CSymbolic_free, p ) ;
}

VALUE wrap_CSymbolic_initialize( VALUE self,
    VALUE _norb_symbol, VALUE _norb_actual, VALUE _ncsf_symbol, VALUE _ncsf_actual,
    VALUE _size_ij, VALUE _size_addr, VALUE _size_pqrs, VALUE _size_coef,
    VALUE _n_1el_symbol, VALUE _n_1el, VALUE _n_ext, 
    VALUE _mo_val_range_last, VALUE _mo_ext_range_first, VALUE _mo_ext_range_last )

{
    int i, j, id, depth, max_nloop ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;

    p->norb_symbol        = FIX2INT( _norb_symbol ) ;
    p->norb_actual        = FIX2INT( _norb_actual ) ;
    p->ncsf_symbol        = FIX2INT( _ncsf_symbol ) ;
    p->ncsf_actual        = FIX2INT( _ncsf_actual ) ;
    p->size_ij            = FIX2INT( _size_ij ) ;
    p->size_addr          = FIX2INT( _size_addr ) ;
    p->size_pqrs          = FIX2INT( _size_pqrs ) ;
    p->size_coef          = FIX2INT( _size_coef ) ;
    p->n_1el_symbol       = FIX2INT( _n_1el_symbol ) ;
    p->n_1el              = FIX2INT( _n_1el ) ;
    p->n_ext              = FIX2INT( _n_ext ) ;
    p->mo_val_range_last  = FIX2INT( _mo_val_range_last ) ;
    p->mo_ext_range_first = FIX2INT( _mo_ext_range_first ) ;
    p->mo_ext_range_last  = FIX2INT( _mo_ext_range_last ) ;

    p->csf_symbol  = ( int    ** ) ruby_xmalloc( p->ncsf_symbol * sizeof( int *  ) ) ;
    p->csf_actual  = ( int    ** ) ruby_xmalloc( p->ncsf_actual * sizeof( int *  ) ) ;
    p->csf_symbolc = ( char   ** ) ruby_xmalloc( p->ncsf_symbol * sizeof( int *  ) ) ;
    p->csf_actualc = ( char   ** ) ruby_xmalloc( p->ncsf_actual * sizeof( int *  ) ) ;
    for ( i = 0 ; i < p->ncsf_symbol ; i++ ) {
        p->csf_symbol [ i ] = ( int  * ) ruby_xmalloc( p->norb_symbol * sizeof( int  ) ) ;
        p->csf_symbolc[ i ] = ( char * ) ruby_xmalloc( p->norb_symbol * sizeof( char ) ) ;
    }
    for ( i = 0 ; i < p->ncsf_actual ; i++ ) {
        p->csf_actual [ i ] = ( int  * ) ruby_xmalloc( p->norb_actual * sizeof( int  ) ) ;
        p->csf_actualc[ i ] = ( char * ) ruby_xmalloc( p->norb_actual * sizeof( char ) ) ;
    }
    p->ij         = ( int    *  ) ruby_xmalloc( p->size_ij     * sizeof( int    ) ) ;
    p->addr       = ( int    *  ) ruby_xmalloc( p->size_addr   * sizeof( int    ) ) ;
    p->pqrs       = ( int    *  ) ruby_xmalloc( p->size_pqrs   * sizeof( int    ) ) ;
    p->coef       = ( double *  ) ruby_xmalloc( p->size_coef   * sizeof( double ) ) ;
    p->tri        = ( int    ** ) ruby_xmalloc( ( p->ncsf_symbol * ( p->ncsf_symbol + 1 ) / 2 ) * sizeof( int * ) ) ;
    for ( i = 0 ; i < ( p->ncsf_symbol * ( p->ncsf_symbol + 1 ) / 2 ) ; i++ ) {
        p->tri[ i ] = ( int  *  ) ruby_xmalloc( 2              * sizeof( int    ) ) ;
    }
    for ( i = 1, id = 0 ; i <= p->ncsf_symbol ; i++ ) {
        for ( j = 1 ; j <= i ; j++, id++ ) {
            p->tri[ id ][ 0 ] = i ;
            p->tri[ id ][ 1 ] = j ;
        }
    }
    p->expr_actual  = ( csymbolic_anexpression_t * ) ruby_xmalloc( _MAX_NSIZE_EXPR_ * sizeof( csymbolic_anexpression_t ) ) ;;
    p->nexpr_actual = 0 ;
    p->mos          = ( multiloop_t * ) ruby_xmalloc( ( _MAX_DEPTH_ + 1 ) * sizeof( multiloop_t ) ) ;
    for ( depth = 0 ; depth <= _MAX_DEPTH_ ; depth++ ) {
        p->mos[ depth ].res = ( int * ) ruby_xmalloc(  MAX_NLOOP[ depth ] * sizeof( multiloop_t ) ) ;
    }
    for( depth = 0 ; depth <= _MAX_DEPTH_ ; depth++ ) {
        int val[ _MAX_DEPTH_+1 ] ;
        for ( i = 0 ; i < depth ; i++ ) { val[ i ] = 0 ; }
        p->mos[ depth ].depth = depth ;
        p->mos[ depth ].nmin  = p->mo_ext_range_first  ;
        p->mos[ depth ].nmax  = p->mo_ext_range_last   ;
        p->mos[ depth ].nloop = 0 ;
        gen_loop( &(p->mos[ depth ]), 0, val ) ;
        ///fprintf( stdout, "depth, nloop: %4d, %4d\n", p->mos[ depth ].depth, p->mos[ depth ].nloop ) ;
        ///for ( i = 0 ; i < p->mos[ depth ].nloop ; i++ ) {
        ///    fprintf( stdout, "%6d :", i ) ; for ( j = 0 ; j < depth ; j++ ) { fprintf( stdout, "%4d", p->mos[ depth ].res[ _MAX_DEPTH_ * i + j ] ) ; } fprintf( stdout, "\n" ) ;
        ///}
    }
///
/*
    fprintf( stdout, "c: p->norb_symbol        : %4d\n", p->norb_symbol        ) ;
    fprintf( stdout, "c: p->norb_actual        : %4d\n", p->norb_actual        ) ;
    fprintf( stdout, "c: p->ncsf_symbol        : %4d\n", p->ncsf_symbol        ) ;
    fprintf( stdout, "c: p->ncsf_actual        : %4d\n", p->ncsf_actual        ) ;
    fprintf( stdout, "c: p->size_ij            : %4d\n", p->size_ij            ) ;
    fprintf( stdout, "c: p->size_addr          : %4d\n", p->size_addr          ) ;
    fprintf( stdout, "c: p->size_pqrs          : %4d\n", p->size_pqrs          ) ;
    fprintf( stdout, "c: p->size_coef          : %4d\n", p->size_coef          ) ;
    fprintf( stdout, "c: p->n_1el_symbol       : %4d\n", p->n_1el_symbol       ) ;
    fprintf( stdout, "c: p->n_1el              : %4d\n", p->n_1el              ) ;
    fprintf( stdout, "c: p->n_ext              : %4d\n", p->n_ext              ) ;
    fprintf( stdout, "c: p->mo_val_range_last  : %4d\n", p->mo_val_range_last  ) ;
    fprintf( stdout, "c: p->mo_ext_range_first : %4d\n", p->mo_ext_range_first ) ;
    fprintf( stdout, "c: p->mo_ext_range_last  : %4d\n", p->mo_ext_range_last  ) ;
*/
///    for ( id = 0 ; id <  p->ncsf_symbol * ( p->ncsf_symbol + 1 ) / 2 ; id++ ) {
///        fprintf( stdout, "id: i, j = %10d: %6d %6d\n", id, p->tri[ id ][ 0 ], p->tri[ id ][ 1 ] ) ;
///    }
	return Qnil;
}

VALUE wrap_CSymbolic_set_csf_symbol( VALUE self, VALUE _csf_symbol )
{
    int i, j ;
	csymbolic_t *p;
///
	Data_Get_Struct( self, csymbolic_t, p ) ;
///
    for ( i = 0 ; i < p->ncsf_symbol ; i++ ) {
        VALUE tmp_ary = rb_ary_entry( _csf_symbol, i ) ;
        for ( j = 0 ; j < p->norb_symbol ; j++ ) {
            VALUE tmp = rb_ary_entry( tmp_ary, j ) ;
            char *str = StringValuePtr( tmp ) ;
            p->csf_symbolc[ i ][ j ] = str[ 0 ] ;
            if        ( str[ 0 ] == 'e' ) {
                p->csf_symbol[ i ][ j ] = 0 ;
            } else if ( str[ 0 ] == 'u' ) {
                p->csf_symbol[ i ][ j ] = 1 ;
            } else if ( str[ 0 ] == 'd' ) {
                p->csf_symbol[ i ][ j ] = 2 ;
            } else if ( str[ 0 ] == 'f' ) {
                p->csf_symbol[ i ][ j ] = 3 ;
            } else {
                fprintf( stderr, "Error: func: wrap_CSymbolic_set_csf_symbol: csf_list conversion : (i, j, char) = (%d, %d, %c)\n", i, j, str[ 0 ] ) ;
                exit( 1 ) ;
            }
        }
    }
/*
    for ( i = 0 ; i < p->ncsf_symbol ; i++ ) {
        fprintf( stdout, "csf_symbol: %10d:", i ) ;
        for ( j = 0 ; j < p->norb_symbol ; j++ ) {
            fprintf( stdout, "%2c", p->csf_symbolc[ i ][ j ] ) ;
        }
        fprintf( stdout, " -> " ) ;
        for ( j = 0 ; j < p->norb_symbol ; j++ ) {
            fprintf( stdout, "%2d", p->csf_symbol [ i ][ j ] ) ;
        }
        fprintf( stdout, "\n" ) ;
    }
*/
	return Qnil;
}

VALUE wrap_CSymbolic_set_csf_actual( VALUE self, VALUE _csf_actual )
{
    int i, j ;
	csymbolic_t *p;
///
	Data_Get_Struct( self, csymbolic_t, p ) ;
///
    for ( i = 0 ; i < p->ncsf_actual ; i++ ) {
        VALUE tmp_ary = rb_ary_entry( _csf_actual, i ) ;
        for ( j = 0 ; j < p->norb_actual ; j++ ) {
            VALUE tmp = rb_ary_entry( tmp_ary, j ) ;
            char *str = StringValuePtr( tmp ) ;
            p->csf_actualc[ i ][ j ] = str[ 0 ] ;
            if        ( str[ 0 ] == 'e' ) {
                p->csf_actual[ i ][ j ] = 0 ;
            } else if ( str[ 0 ] == 'u' ) {
                p->csf_actual[ i ][ j ] = 1 ;
            } else if ( str[ 0 ] == 'd' ) {
                p->csf_actual[ i ][ j ] = 2 ;
            } else if ( str[ 0 ] == 'f' ) {
                p->csf_actual[ i ][ j ] = 3 ;
            } else {
                fprintf( stderr, "Error: func: wrap_CSymbolic_csf_actual: csf_list conversion : (i, j, char) = (%d, %d, %c)\n", i, j, str[ 0 ] ) ;
                exit( 1 ) ;
            }
        }
    }
/*
    for ( i = 0 ; i < p->ncsf_actual ; i++ ) {
        fprintf( stdout, "csf_actual: %10d:", i ) ;
        for ( j = 0 ; j < p->norb_actual ; j++ ) {
            fprintf( stdout, "%2c", p->csf_actualc[ i ][ j ] ) ;
        }
        fprintf( stdout, " -> " ) ;
        for ( j = 0 ; j < p->norb_actual ; j++ ) {
            fprintf( stdout, "%2d", p->csf_actual [ i ][ j ] ) ;
        }
        fprintf( stdout, "\n" ) ;
    }
*/
	return Qnil;
}

VALUE wrap_CSymbolic_set_ij( VALUE self, VALUE _ij )
{
    int i ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    for ( i = 0 ; i < p->size_ij ; i++ ) {
        VALUE tmp  = rb_ary_entry( _ij, i ) ;
        p->ij[ i ] = FIX2INT( tmp ) ;
    }
///    for ( i = 0 ; i < p->size_ij ; i++ ) {
///        fprintf( stdout, "ij: %6d%18d\n", i, p->ij[ i ] ) ;
///    }
	return Qnil;
}

VALUE wrap_CSymbolic_set_addr( VALUE self, VALUE _addr )
{
    int i ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    for ( i = 0 ; i < p->size_addr ; i++ ) {
        VALUE tmp    = rb_ary_entry( _addr, i ) ;
        p->addr[ i ] = FIX2INT( tmp ) ;
    }
///    for ( i = 0 ; i < p->size_addr ; i++ ) {
///        fprintf( stdout, "addr: %6d%18d\n", i, p->addr[ i ] ) ;
///    }
	return Qnil;
}

VALUE wrap_CSymbolic_set_pqrs( VALUE self, VALUE _pqrs )
{
    int i ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    for ( i = 0 ; i < p->size_pqrs ; i++ ) {
        VALUE tmp    = rb_ary_entry( _pqrs, i ) ;
        p->pqrs[ i ] = FIX2INT( tmp ) ;
    }
///    for ( i = 0 ; i < p->size_pqrs ; i++ ) {
///        fprintf( stdout, "pqrs: %10d: %10d\n", i, p->pqrs[ i ] ) ;
///    }
	return Qnil;
}

VALUE wrap_CSymbolic_set_coef( VALUE self, VALUE _coef )
{
    int i ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    for ( i = 0 ; i < p->size_coef ; i++ ) {
        VALUE tmp    = rb_ary_entry( _coef, i ) ;
        p->coef[ i ] = rb_num2dbl( tmp ) ;
    }
///    for ( i = 0 ; i < p->size_coef ; i++ ) {
///        fprintf( stdout, "pqrs, coef: %6d%18d%23.8e\n", i, p->pqrs[ i ], p->coef[ i ] ) ;
///    }
	return Qnil;
}

VALUE wrap_CSymbolic_set_drt( VALUE self, VALUE _nelec, VALUE _spin2, VALUE _n_core, VALUE _n_val, VALUE _n_ext )
{
    int i, nelec, spin2, n_core, n_val, n_ext ;
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    nelec  = FIX2INT( _nelec  ) ;
    spin2  = FIX2INT( _spin2  ) ;
    n_core = FIX2INT( _n_core ) ;
    n_val  = FIX2INT( _n_val  ) ;
    n_ext  = FIX2INT( _n_ext  ) ;
    p->drt = edml_drt_create( nelec, spin2, n_core, n_val, n_ext ) ;
///    fprintf( stdout, "C-DRT:\n" ) ;
///    edml_drt_show( p->drt, stdout ) ;
	return Qnil;
}

VALUE wrap_CSymbolic_actual_expression( VALUE self )
{
	csymbolic_t *p;
	Data_Get_Struct( self, csymbolic_t, p ) ;
    main_actual_expression( p ) ;
	return Qnil;
}

static int dump_expr( csymbolic_t *p )
{
    int  i, j ;
    FILE *fp       = fopen( "Data_expr_c", "w" ) ;
    fprintf( fp, "%d %d %d\n", p->norb_actual, p->ncsf_actual, p->nexpr_actual ) ;
    for ( i = 0 ; i < p->ncsf_actual ; i++ ) {
        for ( j = 0 ; j < p->norb_actual ; j++ ) {
             fprintf( fp, "%2c", p->csf_actualc[ i ][ j ] ) ;
        }
         fprintf( fp, "\n" ) ;
    }
    for ( i = 0 ; i < p->nexpr_actual ; i++ ) {
        fprintf( fp, "%d %d %23.16e\n", (p->expr_actual[ i ]).ij, (p->expr_actual[ i ]).pqrs, (p->expr_actual[ i ]).coef ) ;
    }
    fclose( fp ) ;
///
    return 0 ;
}

VALUE wrap_CSymbolic_get_expression( VALUE self, VALUE _expr_actual_array )
{
    int i ;
	csymbolic_t *p;
    VALUE ary_ij   = rb_ary_new() ;
    VALUE ary_pqrs = rb_ary_new() ;
    VALUE ary_coef = rb_ary_new() ;
///
	Data_Get_Struct( self, csymbolic_t, p ) ;
///
///    dump_expr( p ) ;
///
    for ( i = 0 ; i < p->nexpr_actual ; i++ ) {
        rb_ary_push( ary_ij,   INT2NUM(      (p->expr_actual[ i ]).ij   ) ) ;
        rb_ary_push( ary_pqrs, INT2NUM(      (p->expr_actual[ i ]).pqrs ) ) ;
        rb_ary_push( ary_coef, rb_float_new( (p->expr_actual[ i ]).coef ) ) ;
    }
    rb_ary_push( _expr_actual_array, ary_ij   ) ;
    rb_ary_push( _expr_actual_array, ary_pqrs ) ;
    rb_ary_push( _expr_actual_array, ary_coef ) ;
///
 fprintf( stderr, "return: get_expression.\n" ) ;
	return Qnil;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
void Init_CSymbolic ( )
{
    VALUE c = rb_define_class( "CSymbolic",           rb_cObject ) ;
///
    rb_define_alloc_func    ( c,                      wrap_CSymbolic_allocate ) ;
    rb_define_private_method( c, "initialize",        wrap_CSymbolic_initialize,        14 ) ;
    rb_define_method(         c, "set_csf_symbol",    wrap_CSymbolic_set_csf_symbol,     1 ) ;
    rb_define_method(         c, "set_csf_actual",    wrap_CSymbolic_set_csf_actual,     1 ) ;
    rb_define_method(         c, "set_ij",            wrap_CSymbolic_set_ij,             1 ) ;
    rb_define_method(         c, "set_addr",          wrap_CSymbolic_set_addr,           1 ) ;
    rb_define_method(         c, "set_pqrs",          wrap_CSymbolic_set_pqrs,           1 ) ;
    rb_define_method(         c, "set_coef",          wrap_CSymbolic_set_coef,           1 ) ;
    rb_define_method(         c, "set_drt",           wrap_CSymbolic_set_drt,            5 ) ;
    rb_define_method(         c, "actual_expression", wrap_CSymbolic_actual_expression,  0 ) ;
    rb_define_method(         c, "get_expression",    wrap_CSymbolic_get_expression,     1 ) ;
}
