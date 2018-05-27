#include <stdio.h>
#include <stdlib.h>
#include "ruby.h"
#include "csymbolic.h"

FILE  *fp ;
const char *filename = "Dump_s2e_c" ;

static int get_pq( csymbolic_t *csymb, int pq, int *p_ptr, int *q_ptr )
{
///    fprintf ( stdout, "pq, p, q: %6d %6d %6d\n", pq, (csymb->tri[ pq ])[ 0 ], (csymb->tri[ pq ])[ 1 ] ) ;
    *p_ptr = csymb->tri[ pq ][ 0 ] ;
    *q_ptr = csymb->tri[ pq ][ 1 ] ;
    return 0 ;
}

static int mo_order( csymbolic_t *csymb, int mode, char *f, char *g,
              int *pdepth, int *f_order, int *g_order, char *f_state, char *g_state,
              int *len_fo, int *len_go, int *len_fs, int *len_gs )
{
    int i, k, depth, pos_fo, pos_go, pos_fs, pos_gs ;
    char cond ;
    if ( mode == 1 ) {
        cond = 'f' ;
    } else {
        cond = 'e' ;
    }
    pos_fo = pos_go = pos_fs = pos_gs = 0 ;
    depth  = 0 ;
    for ( k = csymb->mo_val_range_last + 1 ; k < csymb->mo_val_range_last + csymb->n_ext + 1 ; k++ ) {
        if ( f[ k - 1 ] != cond ) {
            f_order[ pos_fo++ ] = depth ;
            f_state[ pos_fs++ ] = f[ k - 1 ] ;
            if ( g[ k - 1 ] != cond ) {
                g_order[ pos_go++ ] = depth ;
                g_state[ pos_gs++ ] = g[ k - 1 ] ;
            }
            depth++ ;
        } else if ( g[ k - 1 ] != cond ) {
            g_order[ pos_go++ ] = depth ;
            g_state[ pos_gs++ ] = g[ k - 1 ] ;
            depth++ ;
        }
    }
    *pdepth = depth ;
    *len_fo = pos_fo ;
    *len_go = pos_go ;
    *len_fs = pos_fs ;
    *len_gs = pos_gs ;
/*
    fprintf( fp, "depth, f_order, g_order, f_state, g_state: %6d", *pdepth ) ;
    for ( k = 0 ; k < pos_fo ; k++ ) { fprintf( fp, "%3d", f_order[ k ] ) ; }
    fprintf( fp, ", " ) ;
    for ( k = 0 ; k < pos_go ; k++ ) { fprintf( fp, "%3d", g_order[ k ] ) ; }
    fprintf( fp, ", " ) ;
    for ( k = 0 ; k < pos_fs ; k++ ) { fprintf( fp, "%3c", f_state[ k ] ) ; }
    fprintf( fp, ", " ) ;
    for ( k = 0 ; k < pos_gs ; k++ ) { fprintf( fp, "%3c", g_state[ k ] ) ; }
    fprintf( fp, "\n" ) ;
*/
    return 0 ;
}

static int csf_upper_arc_weight( csymbolic_t *csymb, char *csf, int range_begin, int range_end, int _inode, int _id, int *inode_id )
{
    int i, j, walk ;
    int inode = _inode ;
    int id    = _id;
    for ( i = range_begin ; i < range_end ; i++ ) {
        if      ( csf[ i ] == 'e' ) { walk = 0 ; }
        else if ( csf[ i ] == 'u' ) { walk = 1 ; }
        else if ( csf[ i ] == 'd' ) { walk = 2 ; }
        else if ( csf[ i ] == 'f' ) { walk = 3 ; }
        else                        { fprintf( stderr, "Error: csf_upper_arc_weight(): csf[ %d ] = %c != eudf.\n", i, csf[ i ] ) ;
                                      exit ( 1 ) ; }
        id   += edml_drt_get_upper_arcweight( csymb->drt, inode, walk ) ;
        inode = edml_drt_get_upper_arc( csymb->drt, inode, walk ) ;
    }
    inode_id[ 0 ] = inode ;
    inode_id[ 1 ] = id ;
    return 0 ;
}

static int to_actual_mo( int *map, int x0, int x )
{
    int n ;
    if        ( x == ( x0 + 1 ) ) {
        n = map[ 0 ] ;
    } else if ( x == ( x0 + 2 ) ) {
        n = map[ 1 ] ;
    } else if ( x == ( x0 + 3 ) ) {
        n = map[ 2 ] ;
    } else if ( x == ( x0 + 4 ) ) {
        n = map[ 3 ] ;
    } else {
        n = x ;
    }
    return n ;
}

static int pq_tri( int p, int q )
{
    int n, max, min ;
    if ( p > q ) {
        max = p ;
        min = q ;
    } else {
        max = q ;
        min = p ;
    }
    return ( ( max * ( max - 1 ) ) / 2 + min ) ;
}

static int pqrs_tritri( int p, int q, int r, int s )
{
    int npq   = pq_tri ( p, q ) ;
    int nrs   = pq_tri ( r, s ) ;
    int npqrs = pq_tri ( npq, nrs ) ;
    return npqrs ;
}

int main_actual_expression( csymbolic_t *csymb )
{
    int  i, k, ij, ij_k, it, iloop, icsf, jcsf, depth ;
    int  *f_order, *g_order ;
    char *f, *g, *f_state, *g_state, *f_actual, *g_actual ;
    int  len_fo, len_go, len_fs, len_gs ;
    int  f_inode0, g_inode0, f_id0, g_id0, f_inode, g_inode, f_id, g_id, inode_id[ 2 ] ;
    int  ij_actual, p, q, r, s, pq, rs, p_actual, q_actual, r_actual, s_actual, pqrs_actual ;
    slib_pnode_t node ;

    multiloop_t **mos ;
///
/// fp = fopen( filename, "w" ) ;
///
    f_order  = ( int  * ) ruby_xmalloc( ( csymb->norb_symbol ) * sizeof( int  ) ) ;
    g_order  = ( int  * ) ruby_xmalloc( ( csymb->norb_symbol ) * sizeof( int  ) ) ;
    f_state  = ( char * ) ruby_xmalloc( ( csymb->norb_symbol ) * sizeof( char ) ) ;
    g_state  = ( char * ) ruby_xmalloc( ( csymb->norb_symbol ) * sizeof( char ) ) ;
    f_actual = ( char * ) ruby_xmalloc(     csymb->norb_actual * sizeof( char ) ) ;
    g_actual = ( char * ) ruby_xmalloc(     csymb->norb_actual * sizeof( char ) ) ;
///
    for ( ij_k = 0 ; ij_k < csymb->size_ij ; ij_k++ ) {
        ij   = csymb->ij[ ij_k ] ;
        icsf = csymb->tri[ ij - 1 ][ 0 ] - 1 ;
        jcsf = csymb->tri[ ij - 1 ][ 1 ] - 1 ;
        f    = csymb->csf_symbolc[ icsf ] ;
        g    = csymb->csf_symbolc[ jcsf ] ;
///
///     fprintf( fp, "ij, ij_k, %6d %6d:", ij, ij_k ) ;
///     for( i = 0 ; i < csymb->norb_symbol ; i++ ) {
///         fprintf( fp, "%2c", f[ i ] ) ;
///     }
///     fprintf( fp, ", " ) ;
///     for( i = 0 ; i < csymb->norb_symbol ; i++ ) {
///         fprintf( fp, "%2c", g[ i ] ) ;
///     }
///     fprintf( fp, "\n" ) ;
///
        mo_order( csymb, 2, f, g, &depth, f_order, g_order, f_state, g_state, &len_fo, &len_go, &len_fs, &len_gs ) ;

        for ( i = 0 ; i < csymb->mo_val_range_last ; i++ ) {
            f_actual[ i ] = f[ i ] ;
            g_actual[ i ] = g[ i ] ;
        }
        for ( i = csymb->mo_val_range_last ; i < csymb->mo_ext_range_last ; i++ ) {
            f_actual[ i ] = 'e' ;
            g_actual[ i ] = 'e' ;
        }
///     fprintf( fp, "f_actual: original: [" ) ;
///     for ( k = 0 ; k < csymb->norb_actual ; k++ ) {
///         fprintf( fp, "%2c", f_actual[ k ] ) ; 
///     }
///     fprintf( fp, "], [" ) ;
///     for ( k = 0 ; k < csymb->norb_actual ; k++ ) {
///         fprintf( fp, "%2c", g_actual[ k ] ) ; 
///     }
///     fprintf( fp, "]\n" ) ;

        csf_upper_arc_weight( csymb, f_actual, 0, csymb->mo_val_range_last, 0, 1, inode_id ) ;
        f_inode0 = inode_id[ 0 ] ;
        f_id0    = inode_id[ 1 ] ;
        csf_upper_arc_weight( csymb, g_actual, 0, csymb->mo_val_range_last, 0, 1, inode_id ) ;
        g_inode0 = inode_id[ 0 ] ;
        g_id0    = inode_id[ 1 ] ;
///        fprintf( fp, "f_inode0, f_id0, g_inode0, g_id0: %4d%4d, %4d%4d\n", f_inode0, f_id0, g_inode0, g_id0 ) ;
///        fflush ( fp ) ;

        for ( iloop = 0 ; iloop < (csymb->mos[ depth ]).nloop ; iloop++ ) {
            int *x = (csymb->mos[ depth ]).res + ( _MAX_DEPTH_ * iloop ) ;
///         fprintf( fp, "pnode: %4d:", iloop ) ; for ( i = 0 ; i < depth ; i++ ) { fprintf( fp, "%4d", x[ i ] ) ; } fprintf( fp, "\n" ) ;
///         fprintf( fp, "f_actual before: [" ) ; for ( k = 0 ; k < csymb->norb_actual ; k++ ) { fprintf( fp, "%2c", f_actual[ k ] ) ; } fprintf( fp, "]\n" ) ;
///         fprintf( fp, "g_actual before: [" ) ; for ( k = 0 ; k < csymb->norb_actual ; k++ ) { fprintf( fp, "%2c", g_actual[ k ] ) ; } fprintf( fp, "]\n" ) ;
///
            for ( k = 0 ; k < len_fo ; k++ ) {
                f_actual[ x[ f_order[ k ] ] - 1 ] = f_state[ k ] ;
            }
            for ( k = 0 ; k < len_go ; k++ ) {
                g_actual[ x[ g_order[ k ] ] - 1 ] = g_state[ k ] ;
            }
///
///         fprintf( fp, "fg_actual: [" ) ;
///         for ( k = 0 ; k < csymb->norb_actual ; k++ ) { fprintf( fp, "%2c", f_actual[ k ] ) ; } fprintf( fp, "], " ) ;
///         fprintf( fp, "[" ) ;
///         for ( k = 0 ; k < csymb->norb_actual ; k++ ) { fprintf( fp, "%2c", g_actual[ k ] ) ; } fprintf( fp, "]\n" ) ;
///         fprintf( fp, "range: %d, %d\n", csymb->mo_val_range_last, csymb->mo_ext_range_last ) ;
///
            csf_upper_arc_weight( csymb, f_actual, csymb->mo_val_range_last, csymb->mo_ext_range_last, f_inode0, f_id0, inode_id ) ;
            f_inode = inode_id[ 0 ] ;
            f_id    = inode_id[ 1 ] ;
            csf_upper_arc_weight( csymb, g_actual, csymb->mo_val_range_last, csymb->mo_ext_range_last, g_inode0, g_id0, inode_id ) ;
            g_inode = inode_id[ 0 ] ;
            g_id    = inode_id[ 1 ] ;
///         fprintf( fp, "f_inode, f_id, g_inode, g_id: %4d%4d, %4d%4d\n", f_inode, f_id, g_inode, g_id ) ;
///
            ij_actual = pq_tri( f_id, g_id ) ;
///            fprintf( fp, "ij_actual: %d\n", ij_actual ) ;
            for ( it = csymb->addr[ ij_k ] ; it < csymb->addr[ ij_k + 1 ] ; it++ ) {
///#printf( "%5d %5d  %5d  %24.14f\n", @ij[ it ], ij_actual, @pqrs[ it ], @coef[ it ] )
///$fout.printf( "it, @ij[it], ij_actual, @pqrs[ it ], @coef[ it ]: %s: %s %s %s %s\n", it, @ij, ij_actual, @pqrs, @coef )
                if ( csymb->pqrs[ it ] >= csymb->n_1el_symbol ) {
                        get_pq( csymb, csymb->pqrs[ it ] - csymb->n_1el_symbol, &pq, &rs ) ;
                        get_pq( csymb, pq - 1, &p, &q ) ;
                        get_pq( csymb, rs - 1, &r, &s ) ;
                        p_actual = to_actual_mo( x, csymb->mo_val_range_last, p ) ;
                        q_actual = to_actual_mo( x, csymb->mo_val_range_last, q ) ;
                        r_actual = to_actual_mo( x, csymb->mo_val_range_last, r ) ;
                        s_actual = to_actual_mo( x, csymb->mo_val_range_last, s ) ;
                        pqrs_actual = pqrs_tritri( p_actual, q_actual, r_actual, s_actual ) + csymb->n_1el - 1 ;
///                        fprintf( fp, "if1:  %4d: %6d %6d %6d %6d : %6d %6d : %6d %6d %6d %6d : %6d\n",
///                                         it, p, q, r, s, pq, rs, p_actual, q_actual, r_actual, s_actual, pqrs_actual ) ;
///#printf( "%5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
///$fout.printf( "pqrs, pqrs_actual:%5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
                } else {
                        get_pq( csymb, csymb->pqrs[ it ], &p, &q ) ;
                        p_actual = to_actual_mo( x, csymb->mo_val_range_last, p ) ;
                        q_actual = to_actual_mo( x, csymb->mo_val_range_last, q ) ;
                        pqrs_actual = pq_tri( p_actual, q_actual ) - 1 ;
///                        fprintf( fp, "else: %4d: %6d %6d %6d %6d : %6d\n", it, p, q, p_actual, q_actual, pqrs_actual ) ;
                }
///             fprintf( fp, "ij_actual, pqrs_actual, coef[ it ]: %10d %10d %23.16e\n", ij_actual, pqrs_actual, csymb->coef[ it ] ) ;
                if ( csymb->nexpr_actual < _MAX_NSIZE_EXPR_ ) {
                    (csymb->expr_actual[ csymb->nexpr_actual ]).ij   = ij_actual ;
                    (csymb->expr_actual[ csymb->nexpr_actual ]).pqrs = pqrs_actual ;
                    (csymb->expr_actual[ csymb->nexpr_actual ]).coef = csymb->coef[ it ] ;
                    (csymb->nexpr_actual)++ ;
                } else {
                    fprintf( stderr, "Error message by Developper: Number of expressions are too large.\n" ) ;
                    fprintf( stderr, "                             Please decrease correlating orbitals: ncore, nactive, nexternals.\n" ) ;
                    fprintf( stderr, "Error: main_symbolic: csymb->nexpr_actual >= _MAX_NSIZE_EXPR_: %d >= %d\n",
                             csymb->nexpr_actual, _MAX_NSIZE_EXPR_ ) ;
                    exit( 1 ) ;
                }
            }
            for ( i = csymb->mo_val_range_last ; i < csymb->mo_ext_range_last ; i++ ) {
                f_actual[ i ] = 'e' ;
                g_actual[ i ] = 'e' ;
            }
        }
    }
    ruby_xfree( f_order  ) ;
    ruby_xfree( g_order  ) ;
    ruby_xfree( f_state  ) ;
    ruby_xfree( g_state  ) ;
    ruby_xfree( f_actual ) ;
    ruby_xfree( g_actual ) ;

/// fclose( fp ) ;
    return 0 ;
}
