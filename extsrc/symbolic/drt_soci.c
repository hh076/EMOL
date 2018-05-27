#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "drt.h"

static char *str_eudf[] = { "e", "u", "d", "f" } ;

static int compare_nodes( const edml_drt_node_t node_a, const edml_drt_node_t node_b )
{
    int *codes_a = node_a->codes ;
    int *codes_b = node_b->codes ;
    if        ( codes_a[ 0 ] > codes_b[ 0 ] ) {
        return +1 ;
    } else if ( codes_a[ 0 ] < codes_b[ 0 ] ) {
        return -1 ;
    } else {
        if        ( codes_a[ 1 ] > codes_b[ 1 ] ) {
            return +1 ;
        } else if ( codes_a[ 1 ] < codes_b[ 1 ] ) {
            return -1 ;
        } else {
            if        ( codes_a[ 2 ] > codes_b[ 2 ] ) {
                return +1 ;
            } else if ( codes_a[ 2 ] < codes_b[ 2 ] ) {
                return -1 ;
            } else {
                return 0 ;
            }
        }
    }
}

static int imax( int i, int j )
{
    int n = i ;
    if ( j > i ) {
        n = j ;
    }
    return n ;
}

/*
int heap_sort( edml_drt_node_t *buff, int size )
{
    int i ;
    for ( i = size / 2 - 1 ; i >= 0 ; i-- ) {
        int c, n = i ;
        edml_drt_node_t x = buff[ n ] ;
        while ( ( c = 2 * n + 1 ) < size ) {
            fprintf( stdout, "a: c, c+1, buff: %4d%4d, %p, %p\n", c, c+1, buff + c, buff + c + 1 ) ;
            if ( ( c + 1 < size ) && ( compare_nodes( buff[ c ], buff[ c + 1 ] ) < 0 ) ) {
                c++ ;
            }
            if ( compare_nodes( x, buff[ c ] ) >= 0 ) {
                break;
            }
            buff[ n ] = buff[ c ] ;
            n = c ;
        }
        buff[ n ] = x ;
    }
    for ( i = size - 1 ; i >= 0 ; i-- ) {
        int c, n = 0;
        edml_drt_node_t x = buff[ i ] ;
        buff[ i ] = buff[ 0 ] ;
        while ( (c = 2 * n + 1) < i ) {
            fprintf( stdout, "b: c, c+1, buff: %4d%4d, %p, %p\n", c, c+1, buff + c, buff + c + 1 ) ;
            if ( ( c + 1 < i ) && ( compare_nodes( buff[ c ], buff[ c + 1 ] ) < 0 ) ) {
                c++ ;
            }
            if ( ( compare_nodes( x, buff[ c ]) >= 0 ) ) {
                break ;
            }
            buff[ n ] = buff[ c ] ;
            n = c ;
        }
        buff[ n ] = x ;
    }
    return 0 ;
}
*/

static int bubble_sort( int nnode, edml_drt_node_t *nodes )
{
    int i, j ;
    for ( i = 0 ; i < nnode - 1 ; i++ ) {
        for ( j = i + 1 ; j < nnode ; j++ ) {
            if ( compare_nodes( nodes[ i ], nodes[ j ] ) > 0 ) {
                edml_drt_node_t tmp = nodes[ i ] ;
                nodes[ i ] = nodes[ j ] ;
                nodes[ j ] = tmp ;
            }
        }
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static edml_drt_info_t edml_drt_info_create( int nelec, int spin2 )
{
    edml_drt_info_t info = ( edml_drt_info_t ) malloc( sizeof( edml_drt_info_item ) ) ;
    info->nelec = nelec ;
    info->spin2 = spin2 ;
    return info ;
}

static int edml_drt_info_free( edml_drt_info_t info )
{
    if ( info != ( edml_drt_info_t )NULL ) {
        free( info ) ;
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static edml_drt_mo_t edml_drt_mo_create( int ncor, int nval, int next )
{
    edml_drt_mo_t mo = ( edml_drt_mo_t ) malloc( sizeof( edml_drt_mo_item ) ) ;
    mo->ncor = ncor ;
    mo->nval = nval ;
    mo->next = next ;
    mo->norb = ncor + nval + next ;
    mo->cor_range_first = 1 ;
    mo->val_range_first = 1 + ncor ;
    mo->ext_range_first = 1 + ncor + nval ;
    mo->cor_range_last  = ncor ;
    mo->val_range_last  = ncor + nval ;
    mo->ext_range_last  = ncor + nval + next ;
    return mo ;
}

static int edml_drt_mo_free( edml_drt_mo_t mo )
{
    if ( mo != ( edml_drt_mo_t )NULL ) {
        free( mo ) ;
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int check_node( edml_drt_info_t info, edml_drt_mo_t mo, edml_drt_node_t node )
{
    int ret ;
    if ( ( ( node->spin2 < 0 ) || ( node->nelec > info->nelec ) || ( node->orb > mo->norb ) ) ||
         ( ( node->orb >= mo->val_range_last ) && ( node->nelec <= info->nelec - 3 ) ) ) {
        ret = 0 ;
    } else {
        ret = 1 ;
    }
    return ret ;
}

static int if_connect_node( edml_drt_node_t node, edml_drt_node_t other_node )
{
    int flg_e = ( ( node->orb + 1 == other_node->orb ) && ( node->nelec     == other_node->nelec ) && ( node->spin2     == other_node->spin2 ) ) ;
    int flg_u = ( ( node->orb + 1 == other_node->orb ) && ( node->nelec + 1 == other_node->nelec ) && ( node->spin2 + 1 == other_node->spin2 ) ) ;
    int flg_d = ( ( node->orb + 1 == other_node->orb ) && ( node->nelec + 1 == other_node->nelec ) && ( node->spin2 - 1 == other_node->spin2 ) ) ;
    int flg_f = ( ( node->orb + 1 == other_node->orb ) && ( node->nelec + 2 == other_node->nelec ) && ( node->spin2     == other_node->spin2 ) ) ;
    return flg_e || flg_u || flg_d || flg_f ;
}

static int if_connect_area( edml_drt_node_t node, int nnode, edml_drt_node_t *nodes )
{
    int i, id ;
    id = 0 ;
    for ( i = 0 ; i < nnode ; i++ ) {
        if ( if_connect_node( node, nodes[ i ] ) ) {
            goto found ;
        }
        id++ ;
    }
    return -1 ;
found:
    return id ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static edml_drt_node_t edml_drt_node_create( edml_drt_info_t info, edml_drt_mo_t mo,
                                      int orb, int nelec, int spin2 )
{
    edml_drt_node_t node = ( edml_drt_node_t ) malloc( sizeof( edml_drt_node_item ) ) ;
    node->arc_weights = ( int * ) malloc( 4 * sizeof( int ) ) ;
    node->upper_arcs  = ( int * ) malloc( 4 * sizeof( int ) ) ;
    node->lower_arcs  = ( int * ) malloc( 4 * sizeof( int ) ) ;
    node->upper_path_weights
                      = ( int * ) malloc( MAX_NWEIGHTS * sizeof ( int ) ) ;
    node->lower_path_weights
                      = ( int * ) malloc( MAX_NWEIGHTS * sizeof ( int ) ) ;
    node->codes       = ( int * ) malloc( 3 * sizeof ( int ) ) ;
    ///
    node->orb         = orb ;
    node->nelec       = nelec ;
    node->spin2       = spin2 ;
    node->codes[ 0 ]  = orb ;
    node->codes[ 1 ]  = nelec ;
    node->codes[ 2 ]  = spin2 ;
    node->number_upper_arcs = -1 ;
    node->number_lower_arcs = -1 ;
    node->n_upper_path_weights = 0 ;
    node->n_lower_path_weights = 0 ;
    // /
    node->status      = check_node( info, mo, node ) ;
    ///
    node->arc_weights[ 0 ] = -1 ; node->arc_weights[ 1 ] = -1 ; node->arc_weights[ 2 ] = -1 ; node->arc_weights[ 3 ] = -1 ;
    node->upper_arcs [ 0 ] = -1 ; node->upper_arcs [ 1 ] = -1 ; node->upper_arcs [ 2 ] = -1 ; node->upper_arcs [ 3 ] = -1 ;
    node->lower_arcs [ 0 ] = -1 ; node->lower_arcs [ 1 ] = -1 ; node->lower_arcs [ 2 ] = -1 ; node->lower_arcs [ 3 ] = -1 ;
    ///
    return node ;
}

static int edml_drt_node_free( edml_drt_node_t node )
{
    if ( node != ( edml_drt_node_t ) NULL ) {
        if ( node->arc_weights != ( int * ) NULL ) {
            free( node->arc_weights ) ;
        }
        if ( node->upper_arcs != ( int * ) NULL ) {
            free( node->upper_arcs ) ;
        }
        if ( node->lower_arcs != ( int * ) NULL ) {
            free( node->lower_arcs ) ;
        }
        if ( node->codes != ( int * ) NULL ) {
            free( node->codes ) ;
        }
        free( node ) ;
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int edml_drt_cor_nodes_create( edml_drt_info_t info, edml_drt_mo_t mo,
                               int *nnode,          edml_drt_node_t *nodes,
                               int *nnode_cor_last, edml_drt_node_t *nodes_cor_last )
{
    int k, ne ;
    ///
    nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, 0, 0,    0 ) ;
    ///
    nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, 1, 0,    0 ) ; /// bottom-empty
    nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, 1, 1,    1 ) ; /// bottom-up
    nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, 1, 2,    0 ) ; /// bottom-full
    ///
    for ( k = 2 ; k <= mo->cor_range_last ; k++ ) {
        ne = 2 * k ;
        nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, ne-2, 0 ) ;
        nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, ne-2, 2 ) ;
        nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, ne-1, 1 ) ;
        nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, ne,   0 ) ;
    }
    ///
    k = mo->cor_range_last ;
    if ( k > 1 ) {
        ne = 2 * k ;
        nodes_cor_last[ 0 ] = edml_drt_node_create( info, mo, k, ne-2, 0 ) ;
        nodes_cor_last[ 1 ] = edml_drt_node_create( info, mo, k, ne-2, 2 ) ;
        nodes_cor_last[ 2 ] = edml_drt_node_create( info, mo, k, ne-1, 1 ) ;
        nodes_cor_last[ 3 ] = edml_drt_node_create( info, mo, k, ne,   0 ) ;
        *nnode_cor_last     = 4 ;
    } else {
        nodes_cor_last[ 0 ] = edml_drt_node_create( info, mo, 1, 0,    0 ) ;
        nodes_cor_last[ 1 ] = edml_drt_node_create( info, mo, 1, 1,    1 ) ;
        nodes_cor_last[ 2 ] = edml_drt_node_create( info, mo, 1, 2,    0 ) ;
        *nnode_cor_last     = 3 ;
    }
    return 0 ;
}

/////////
/////////
/////////
static int edml_drt_ext_nodes_create( edml_drt_info_t info, edml_drt_mo_t mo,
                               int *nnode,           edml_drt_node_t *nodes,
                               int *nnode_ext_first, edml_drt_node_t *nodes_ext_first )
{
    int i, k ;
    int nelec  = info->nelec ;
    int spin22 = info->spin2 * 2 ;
    int nnode_first = *nnode ;
///
    for ( k = mo->ext_range_first ; k <= mo->ext_range_last ; k++ ) {
        if ( k < mo->ext_range_last ) {
            nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, k, nelec-2, spin22   ) ;
            if ( k < mo->ext_range_last-1 ) {
                nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, nelec-2, spin22-2 ) ;
                nodes[ (*nnode)++ ] = edml_drt_node_create( info, mo, k, nelec-2, spin22+2 ) ;
            }
            nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, k, nelec-1, spin22-1 ) ;
            nodes[ (*nnode)++ ]     = edml_drt_node_create( info, mo, k, nelec-1, spin22+1 ) ;
        }
        nodes[ (*nnode)++ ]         = edml_drt_node_create( info, mo, k, nelec,   spin22   ) ;
    }
    if ( mo->ext_range_last < mo->ext_range_first ) {
        nodes[ (*nnode)++ ]         = edml_drt_node_create( info, mo, mo->ext_range_last, nelec, spin22 ) ;
    }
///
    *nnode_ext_first = 0 ;
    for ( i = nnode_first ; i < *nnode ; i++ ) {
        if ( ( mo->next > 0 ) && ( nodes[ i ]->orb == mo->ext_range_first ) ) {
            nodes_ext_first[ (*nnode_ext_first)++ ] = edml_drt_node_create( info, mo, nodes[ i ]->orb, nodes[ i ]->nelec, nodes[ i ]->spin2 ) ;
        } else if ( ( mo->next == 0 ) && ( nodes[ i ]->orb == mo->ext_range_last ) ) {
            nodes_ext_first[ (*nnode_ext_first)++ ] = edml_drt_node_create( info, mo, nodes[ i ]->orb, nodes[ i ]->nelec, nodes[ i ]->spin2 ) ;
        }
    }
    return 0 ;
}

/////////
/////////
/////////
static int index_node( edml_drt_node_t node, int nnode, edml_drt_node_t *nodes )
{
    int i, id ;
    id = 0 ;
    for ( i = 0 ; i < nnode ; i++ ) {
        if ( compare_nodes( nodes[ i ], node ) == 0 ) {    // if equal.
            goto found ;
        }
        id++ ;
    }
    return -1 ;
found:
    return id ;
}

static int edml_drt_val_nodes_create( edml_drt_info_t info, edml_drt_mo_t mo, int orb, int nelec, int spin2,
                               int nnode_ext_first, edml_drt_node_t *nodes_ext_first, int *nnode, edml_drt_node_t *nodes )
{
    int inode, orb_last ;
    int empty, up, down, full ;
    edml_drt_node_t new_node ;

    new_node = edml_drt_node_create( info, mo, orb, nelec, spin2 ) ;
    orb_last = mo->val_range_last ;

    if ((inode = index_node( new_node, *nnode, nodes )) >= 0) {
        return nodes[ inode ]->status ;
    }

    if ( new_node->status == 0 ) {
        nodes[ (*nnode)++ ] = new_node ;
        return new_node->status ;
    }

    if ( orb == orb_last ) {
        if ( if_connect_area( new_node, nnode_ext_first, nodes_ext_first ) >= 0 ) {
            new_node->status = 1 ;
            nodes[ (*nnode)++ ] = new_node ;
            return new_node->status ;
        } else {
            new_node->status = 0 ;
            nodes[ (*nnode)++ ] = new_node ;
            return new_node->status ;
        }
    }

    empty = edml_drt_val_nodes_create( info, mo, orb+1, nelec,   spin2,   nnode_ext_first, nodes_ext_first, nnode, nodes ) ;
    up    = edml_drt_val_nodes_create( info, mo, orb+1, nelec+1, spin2+1, nnode_ext_first, nodes_ext_first, nnode, nodes ) ;
    down  = edml_drt_val_nodes_create( info, mo, orb+1, nelec+1, spin2-1, nnode_ext_first, nodes_ext_first, nnode, nodes ) ;
    full  = edml_drt_val_nodes_create( info, mo, orb+1, nelec+2, spin2,   nnode_ext_first, nodes_ext_first, nnode, nodes ) ;

    if (( empty == 0 ) && ( full == 0 ) && ( up == 0 ) && ( down == 0 )) {
        new_node->status = 0 ;
    } else {
        new_node->status = 1 ;
    }
    nodes[ (*nnode)++ ] = new_node ;
    return new_node->status ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int edml_drt_nodes_free( int nnode, edml_drt_node_t *nodes )
{
    if ( nodes != ( edml_drt_node_t * )NULL ) {
        int i ;
        for ( i = 0 ; i < nnode ; i++ ) {
            edml_drt_node_free( nodes[ i ] ) ;
        }
        free( nodes ) ;
    }
    return 0 ;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int set_arcs( edml_drt_info_t info, edml_drt_mo_t mo, int nnode, edml_drt_node_t *nodes )
{
    int i, j, orb, nelec, spin2 ;
    edml_drt_node_t node_eudf[ 4 ] ;

    for ( i = 0 ; i < nnode ; i++ ) {
        orb   = nodes[ i ]->codes[ 0 ] ;
        nelec = nodes[ i ]->codes[ 1 ] ;
        spin2 = nodes[ i ]->codes[ 2 ] ;
        nodes[ i ]->spin2 = spin2 ;
        node_eudf[ 0 ] = edml_drt_node_create( info, mo, orb+1, nelec,   spin2   ) ;
        node_eudf[ 1 ] = edml_drt_node_create( info, mo, orb+1, nelec+1, spin2+1 ) ;
        node_eudf[ 2 ] = edml_drt_node_create( info, mo, orb+1, nelec+1, spin2-1 ) ;
        node_eudf[ 3 ] = edml_drt_node_create( info, mo, orb+1, nelec+2, spin2   ) ;
        for ( j = 0 ; j < 4 ; j++ ) {
            nodes[ i ]->upper_arcs[ j ] = index_node( node_eudf[ j ], nnode, nodes ) ;
            ///fprintf( stdout, "i,j,code, %4d,%4d: %4d%4d%4d: %4d\n",
            ///         i, j, node_eudf[ j ]->codes[ 0 ], node_eudf[ j ]->codes[ 1 ], node_eudf[ j ]->codes[ 2 ], nodes[ i ]->upper_arcs[ j ] ) ;
        }
        for ( j = 0 ; j < 4 ; j++ ) {
            int upper_id = nodes[ i ]->upper_arcs[ j ] ;
            ///fprintf( stdout, "i,j, upper_id: %4d,%4d: %4d\n", i, j, upper_id ) ;
            if ( upper_id >= 0 ) {
                nodes[ upper_id ]->lower_arcs[ j ] = i ;
            }
        }
    }
    ///fprintf( stdout, "nodes: upper, lower_arcs: \n" ) ;
    ///for ( i = 0 ; i < nnode ; i++ ) {
    ///    fprintf( stdout, "%4d: %4d%4d%4d, %4d%4d%4d%4d, %4d%4d%4d%4d\n",
    ///             i, nodes[ i ]->codes[ 0 ], nodes[ i ]->codes[ 1 ], nodes[ i ]->codes[ 2 ],
    ///             nodes[ i ]->upper_arcs[ 0 ], nodes[ i ]->upper_arcs[ 1 ], nodes[ i ]->upper_arcs[ 2 ], nodes[ i ]->upper_arcs[ 3 ],
    ///             nodes[ i ]->lower_arcs[ 0 ], nodes[ i ]->lower_arcs[ 1 ], nodes[ i ]->lower_arcs[ 2 ], nodes[ i ]->lower_arcs[ 3 ] ) ;
    ///}
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

static int store_number_upper_arcs( int inode, int nnode, edml_drt_node_t *nodes )
{
    if ( nodes[ inode ]->number_upper_arcs == -1 ) {
        if ( inode == ( nnode - 1 ) ) {
            nodes[ inode ]->number_upper_arcs = 0 ;
        } else {
            int number_upper_arcs = 0 ;
            int arc_weights = 0 ;
            int j ;
            for ( j = 0 ; j < 4 ; j++ ) {
                int upper_id = nodes[ inode ]->upper_arcs[ j ] ; 
                if ( upper_id != -1 ) {
                    int ret_id = store_number_upper_arcs( upper_id, nnode, nodes ) ;
                    number_upper_arcs += imax( ret_id, 1 ) ;
                    nodes[ inode ]->arc_weights[ j ] = arc_weights ;
                    arc_weights += nodes[ upper_id ]->number_upper_arcs ;
                }
            }
            nodes[ inode ]->number_upper_arcs = number_upper_arcs ;
        }
    }
    return nodes[ inode ]->number_upper_arcs ;
}

static int store_number_lower_arcs( int inode, int nnode, edml_drt_node_t *nodes )
{
    if ( nodes[ inode ]->number_lower_arcs == -1 ) {
        if ( inode == 0 ) {
            nodes[ inode ]->number_lower_arcs = 0 ;
        } else {
            int number_lower_arcs = 0 ;
            int j ;
            for ( j = 0 ; j < 4 ; j++ ) {
                int lower_id = nodes[ inode ]->lower_arcs[ j ] ; 
                if ( lower_id != -1 ) {
                    int ret_id = store_number_lower_arcs( lower_id, nnode, nodes ) ;
                    number_lower_arcs += imax( ret_id, 1 ) ;
                }
            }
            nodes[ inode ]->number_lower_arcs = number_lower_arcs ;
        }
    }
    return nodes[ inode ]->number_lower_arcs ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int tree_search_lower_path( int nnode, edml_drt_node_t *nodes, int inode, int *weights )
{
    int j ;
    if ( inode == ( nnode - 1 ) ) {
        return 0 ;
    }
    for ( j = 0 ; j < 4 ; j++ ) {
        int next_node = nodes[ inode ]->upper_arcs[ j ] ;
        if ( next_node != -1 ) {
            *weights += nodes[ inode ]->arc_weights[ j ] ;
            nodes[ next_node ]->lower_path_weights[ (nodes[ next_node ]->n_lower_path_weights)++ ] = *weights ;
            ///fprintf( stdout, "%4d%4d%4d :", inode, j, next_node ) ;
            ///for ( i = 0 ; i < nodes[ next_node ]->n_lower_path_weights ; i++ ) {
            ///    fprintf( stdout, "%4d,", nodes[ next_node ]->lower_path_weights[ i ] ) ;
            ///}
            ///fprintf( stdout, "\n" ) ;
            tree_search_lower_path( nnode, nodes, next_node, weights ) ;
            *weights -= nodes[ inode ]->arc_weights[ j ] ;
        }
    }
    return 0 ;
}

static int set_lower_path_weights( int nnode, edml_drt_node_t *nodes, int *weights )
{
    *weights = 0 ;
    nodes[ 0 ]->n_lower_path_weights = 0 ;
    nodes[ 0 ]->lower_path_weights[ (nodes[ 0 ]->n_lower_path_weights)++ ] = 0 ;
    tree_search_lower_path( nnode, nodes, 0, weights ) ;
    return 0 ;
}

/////////////////////
/////////////////////
/////////////////////
static int tree_search_upper_path( int nnode, edml_drt_node_t *nodes, int inode, int *weights )
{
    int j ;
    if ( inode >= ( nnode - 1 ) ) {
        nodes[ inode ]->upper_path_weights[ (nodes[ inode ]->n_upper_path_weights)++ ] = *weights ;
        return 0 ;
    }
    for ( j = 0 ; j < 4 ; j++ ) {
        int next_node = nodes[ inode ]->upper_arcs[ j ] ;
        if ( next_node != -1 ) {
            *weights += nodes[ inode ]->arc_weights[ j ] ;
            tree_search_lower_path( nnode, nodes, inode, weights ) ;
            *weights -= nodes[ inode ]->arc_weights[ j ] ;
        }
    }
    return 0 ;
}

static int set_upper_path_weights( int nnode, edml_drt_node_t *nodes, int *weights )
{
    int i ;
    for ( i = 0 ; i < nnode ; i++ ) {
        tree_search_upper_path( nnode, nodes, i, weights ) ;
    }
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
edml_drt_t edml_drt_create( int nelec, int spin2, int ncor, int nval, int next )
{
    int             i, nnode_cor_last, nnode_ext_first, nnode2, weights ;
    edml_drt_t      drt ;
    edml_drt_node_t *nodes_cor_last, *nodes_ext_first, *nodes2 ;
    drt        = ( edml_drt_t ) malloc( sizeof( edml_drt_item ) ) ;
    ///
    ///
    drt->info  = edml_drt_info_create( nelec, spin2 ) ;
    drt->mo    = edml_drt_mo_create( ncor, nval, next ) ;
    drt->nodes = ( edml_drt_node_t * ) malloc( MAX_NNODE * sizeof( edml_drt_node_t ) ) ;
    drt->nnode = 0 ;
    nodes2     = ( edml_drt_node_t * ) malloc( MAX_NNODE * sizeof( edml_drt_node_t ) ) ;
    nnode2     = 0 ;
    ///
    nodes_cor_last  = ( edml_drt_node_t * ) malloc( 4   * sizeof( edml_drt_node_t ) ) ;
    nodes_ext_first = ( edml_drt_node_t * ) malloc( 100 * sizeof( edml_drt_node_t ) ) ;
    edml_drt_cor_nodes_create( drt->info, drt->mo, &drt->nnode, drt->nodes, &nnode_cor_last,  nodes_cor_last  ) ;
    edml_drt_ext_nodes_create( drt->info, drt->mo, &drt->nnode, drt->nodes, &nnode_ext_first, nodes_ext_first ) ;
    ///
    ///fprintf( stdout, "nodes: core, ext:\n" ) ;
    ///for ( i = 0 ; i < drt->nnode ; i++ ) {
    ///    fprintf( stdout, "%4d: %4d%4d%4d\n", i, drt->nodes[ i ]->codes[ 0 ], drt->nodes[ i ]->codes[ 1 ], drt->nodes[ i ]->codes[ 2 ] ) ;
    ///}
    for ( i = 0 ; i < nnode_cor_last ; i++ ) {
        int orb   = nodes_cor_last[ i ]->orb   ;
        int nelec = nodes_cor_last[ i ]->nelec ;
        int spin2 = nodes_cor_last[ i ]->spin2 ;
        edml_drt_val_nodes_create( drt->info, drt->mo, orb+1, nelec,   spin2,   nnode_ext_first, nodes_ext_first, &drt->nnode, drt->nodes ) ;
        edml_drt_val_nodes_create( drt->info, drt->mo, orb+1, nelec+1, spin2+1, nnode_ext_first, nodes_ext_first, &drt->nnode, drt->nodes ) ;
        edml_drt_val_nodes_create( drt->info, drt->mo, orb+1, nelec+1, spin2-1, nnode_ext_first, nodes_ext_first, &drt->nnode, drt->nodes ) ;
        edml_drt_val_nodes_create( drt->info, drt->mo, orb+1, nelec+2, spin2,   nnode_ext_first, nodes_ext_first, &drt->nnode, drt->nodes ) ;
    }
    ///
    ///fprintf( stdout, "nodes: core, ext, val:\n" ) ;
    ///for ( i = 0 ; i < drt->nnode ; i++ ) {
    ///    fprintf( stdout, "%4d: %4d%4d%4d\n", i, drt->nodes[ i ]->codes[ 0 ], drt->nodes[ i ]->codes[ 1 ], drt->nodes[ i ]->codes[ 2 ] ) ;
    ///}

    for ( i = 0 ; i < drt->nnode ; i++ ) {
        if ( drt->nodes[ i ]->status > 0 ) {
            nodes2[ nnode2++ ] = drt->nodes[ i ] ;
        }
    }
    ///fprintf( stdout, "nodes: deleted \n" ) ;
    ///for ( i = 0 ; i < nnode2 ; i++ ) {
    ///    fprintf( stdout, "%4d: %4d%4d%4d\n", i, nodes2[ i ]->codes[ 0 ], nodes2[ i ]->codes[ 1 ], nodes2[ i ]->codes[ 2 ] ) ;
    ///}

    bubble_sort( nnode2, nodes2 ) ;
    {
        ///edml_drt_nodes_free( drt->nnode, drt->nodes ) ;
        free( drt->nodes ) ;
        drt->nnode = nnode2 ;
        drt->nodes = nodes2 ;
    }
#ifdef DEBUG
///    fprintf( stdout, "nodes: sorted \n" ) ;
///    for ( i = 0 ; i < drt->nnode ; i++ ) {
///        fprintf( stdout, "%4d: %4d%4d%4d\n", i, drt->nodes[ i ]->codes[ 0 ], drt->nodes[ i ]->codes[ 1 ], drt->nodes[ i ]->codes[ 2 ] ) ;
///    }
#endif ///DEBUG
    
    set_arcs( drt->info, drt->mo, drt->nnode, drt->nodes ) ;
#ifdef DEBUG
///    fprintf( stdout, "nodes: upper, lower_arcs: \n" ) ;
///    for ( i = 0 ; i < drt->nnode ; i++ ) {
///        fprintf( stdout, "%4d: %4d%4d%4d, %5d%5d%5d%5d, %5d%5d%5d%5d\n",
///                 i, drt->nodes[ i ]->codes[ 0 ], drt->nodes[ i ]->codes[ 1 ], drt->nodes[ i ]->codes[ 2 ],
///                 drt->nodes[ i ]->upper_arcs[ 0 ], drt->nodes[ i ]->upper_arcs[ 1 ], drt->nodes[ i ]->upper_arcs[ 2 ], drt->nodes[ i ]->upper_arcs[ 3 ],
///                 drt->nodes[ i ]->lower_arcs[ 0 ], drt->nodes[ i ]->lower_arcs[ 1 ], drt->nodes[ i ]->lower_arcs[ 2 ], drt->nodes[ i ]->lower_arcs[ 3 ] ) ;
///    }
#endif ///DEBUG

    store_number_upper_arcs( 0,              drt->nnode, drt->nodes ) ;
    store_number_lower_arcs( drt->nnode - 1, drt->nnode, drt->nodes ) ;
    drt->ncsf = drt->nodes[ 0 ]->number_upper_arcs ;
#ifdef DEBUG
///    fprintf( stdout, "nodes: upper, lower_arcs: \n" ) ;
///    for ( i = 0 ; i < drt->nnode ; i++ ) {
///        fprintf( stdout, "%4d: %4d%4d%4d, %5d%5d%5d%5d, %5d%5d%5d%5d, %5d%5d%5d%5d, %6d%6d\n",
///                 i, drt->nodes[ i ]->codes[ 0 ], drt->nodes[ i ]->codes[ 1 ], drt->nodes[ i ]->codes[ 2 ],
///                 drt->nodes[ i ]->upper_arcs[ 0 ], drt->nodes[ i ]->upper_arcs[ 1 ], drt->nodes[ i ]->upper_arcs[ 2 ], drt->nodes[ i ]->upper_arcs[ 3 ],
///                 drt->nodes[ i ]->lower_arcs[ 0 ], drt->nodes[ i ]->lower_arcs[ 1 ], drt->nodes[ i ]->lower_arcs[ 2 ], drt->nodes[ i ]->lower_arcs[ 3 ],
///                 drt->nodes[ i ]->arc_weights[ 0 ], drt->nodes[ i ]->arc_weights[ 1 ], drt->nodes[ i ]->arc_weights[ 2 ], drt->nodes[ i ]->arc_weights[ 3 ],
///                 drt->nodes[ i ]->number_upper_arcs, drt->nodes[ i ]->number_lower_arcs ) ;
///    }
#endif ///DEBUG
    ///
    weights = 0 ;
    set_lower_path_weights( drt->nnode, drt->nodes, &weights ) ;
    ///set_upper_path_weights( drt->nnode, drt->nodes, &weights ) ;
#ifdef DEBUG
///    fprintf( stdout, "nodes: lower_path_weights: \n" ) ;
///    for ( i = 0 ; i < drt->nnode ; i++ ) {
///        int j ;
///        fprintf( stdout, "%4d: ", i ) ;
///        for ( j = 0 ; j < drt->nodes[ i ]->n_lower_path_weights ; j++ ) {
///            fprintf( stdout, "%4d,", drt->nodes[ i ]->lower_path_weights[ j ] ) ;
///        }
///        fprintf( stdout, "\n" ) ;
///    }
///    fprintf( stdout, "----\n" ) ;
///    for ( i = 0 ; i < drt->nnode ; i++ ) {
///        fprintf( stdout, "%4d: ", i ) ;
///        for ( j = 0 ; j < drt->nodes[ i ]->n_upper_path_weights ; j++ ) {
///            fprintf( stdout, "%4d,", drt->nodes[ i ]->upper_path_weights[ j ] ) ;
///        }
///        fprintf( stdout, "\n" ) ;
///    }
#endif ///DEBUG
    ///
    return drt ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
static int tree_search_csf( int inode, int *ncsf, int *acsf, edml_drt_t drt, int *max_ncsf )
{
    int i, j, orb, next_node ;
    if ( inode >= ( drt->nnode - 1 ) ) {
    ///if ( inode <= 0 ) {
        drt->csfs[ *ncsf ] = ( int * )malloc( drt->mo->norb * sizeof( int ) ) ;
        for ( i = 0 ; i < drt->mo->norb ; i++ ) {
            drt->csfs[ *ncsf ][ i ] = acsf[ i ] ;
        }
        (*ncsf)++ ;
        if ( *ncsf > *max_ncsf ) {
            fprintf( stderr, "Errror: edml_drt_gen_csfs: increase ncsf-max-value : %d\n", *max_ncsf ) ;
            exit( 1 ) ;
        }
        return 0 ;
    } else {
        for ( j = 0 ; j < 4 ; j++ ) {
            orb       = drt->nodes[ inode ]->orb ;
            next_node = drt->nodes[ inode ]->upper_arcs[ j ] ;
            ///next_node = drt->nodes[ inode ]->lower_arcs[ j ] ;
            if ( next_node != -1 ) {
                acsf[ orb ] = j ;
                tree_search_csf( next_node, ncsf, acsf, drt, max_ncsf ) ;
            }
        }
    }
    return 0 ;
}

int edml_drt_gen_csfs( edml_drt_t drt )
{
    int ncsf, *acsf ;
    drt->csfs = ( int ** )malloc( drt->ncsf     * sizeof( int * ) ) ;
    acsf      = ( int *  )malloc( drt->mo->norb * sizeof( int   ) ) ;
    ncsf      = 0 ;
    tree_search_csf( 0, &ncsf, acsf, drt, &(drt->ncsf) ) ;
    ///tree_search_csf( drt->nnode-1, &ncsf, acsf, drt, &(drt->ncsf) ) ;
    return ncsf ;
}

int edml_drt_get_ncsf( edml_drt_t drt )
{
    return drt->ncsf ;
}

int edml_drt_get_norb( edml_drt_t drt )
{
    return drt->mo->norb ;
}

int edml_drt_get_csf( edml_drt_t drt, int icsf, int *acsf )
{
    int i ;
    for ( i = 0 ; i < drt->mo->norb ; i++ ) {
        acsf[ i ] = drt->csfs[ icsf ][ i ] ;
    }
    return 0 ;
}

int edml_drt_get_nnode( edml_drt_t drt )
{
    return drt->nnode ;
}

int edml_drt_get_upper_level( edml_drt_t drt, int inode )
{
    return drt->nodes[ inode ]->orb ;
}

int edml_drt_get_upper_arc( edml_drt_t drt, int inode, int walk )
{
    return drt->nodes[ inode ]->upper_arcs[ walk ] ;
}

int edml_drt_get_upper_arcweight( edml_drt_t drt, int inode, int walk )
{
    return drt->nodes[ inode ]->arc_weights[ walk ] ;
}

int edml_drt_get_upper_spin( edml_drt_t drt, int inode )
{
    return drt->nodes[ inode ]->spin2 + 1 ;
}

int edml_drt_show_csfs( edml_drt_t drt, FILE *fp )
{
    int icsf, iorb, *acsf ;
    int ncsf = drt->ncsf ;
    int norb = drt->mo->norb ;
    fprintf( fp, "\n *** CSFs:\n     NCSF, NORB: %8d%8d\n", ncsf, norb ) ;
    acsf = ( int * )malloc( norb * sizeof ( int ) ) ;
    for ( icsf = 0 ; icsf < ncsf ; icsf++ ) {
        edml_drt_get_csf( drt, icsf, acsf ) ;
        fprintf( fp, "%8d: ", icsf+1 ) ;
        for ( iorb = 0 ; iorb < norb ; iorb++ ) {
            fprintf( fp, "%2d", acsf[ iorb ] ) ;
        }
        fprintf( fp, "  (" ) ;
        for ( iorb = 0 ; iorb < norb ; iorb++ ) {
            fprintf( fp, "%2s", str_eudf[ acsf[ iorb ] ] ) ;
        }
        fprintf( fp, " )\n" ) ;
    }
    fprintf( fp, "\n" ) ;
    free( acsf ) ;
    return 0 ;
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int edml_drt_show( edml_drt_t drt, FILE *fp )
{
    int inode ;
    fprintf( fp, "\n" ) ;
    fprintf( fp, " DISTINCT  ROW  TABLE --- SOCI\n\n" ) ;
    fprintf( fp, "   STATE information   NVE =  %4d,   S*2 = %4d\n", drt->info->nelec, drt->info->spin2 ) ;
    fprintf( fp, "   MO classification   NORB = %4d", drt->mo->norb ) ;
    fprintf( fp, "    ( NCORE=%4d, NVAL=%4d, NEXT=%4d )\n\n", drt->mo->ncor, drt->mo->nval, drt->mo->next ) ;
    fprintf( fp, "   Number of generated nodes = %8d\n", drt->nnode ) ;
    fprintf( fp, "   Number of generated csfs  = %8d\n", drt->ncsf ) ;
    fprintf( fp, "\n *** Nodes, Upper Arc Table, and Arc Weights\n" ) ;
    fprintf( fp, "                             upper arcs                #arcs         upper weight\n" ) ;
    fprintf( fp, "     n    mo  el  is         e      u      d      f                  e        u        d        f\n" ) ;
    for ( inode = drt->nnode - 1 ; inode >= 0 ; inode-- ) {
        int j ;
        int *codes            = drt->nodes[ inode ]->codes ;
        int number_upper_arcs = drt->nodes[ inode ]->number_upper_arcs ;
        int *arc_weights      = drt->nodes[ inode ]->arc_weights ;
        int *upper_arcs       = drt->nodes[ inode ]->upper_arcs ;
        ///
        fprintf( fp, "%6d ( %3d %3d %3d ) ", inode, codes[ 0 ], codes[ 1 ], codes[ 2 ] ) ;
        for ( j = 0 ; j < 4 ; j++ ) {
            int ua = edml_drt_get_upper_arc( drt, inode, j ) ;
           /// if ( upper_arcs[ j ] == -1 ) {
            if ( ua == -1 ) {
                fprintf( fp, "      -" ) ;
            } else {
                ///fprintf( fp, "  %5d", upper_arcs[ j ] ) ;
                fprintf( fp, "  %5d", ua ) ;
            }
        }
        fprintf( fp, "  %7d ", number_upper_arcs ) ;
        for ( j = 0 ; j < 4 ; j++ ) {
            int aw = edml_drt_get_upper_arcweight( drt, inode, j ) ;
            ///if ( arc_weights[ j ] == -1 ) {
            if ( aw == -1 ) {
                fprintf( fp, "        -" ) ;
            } else {
                fprintf( fp, "  %7d", aw ) ;
                ///fprintf( fp, "  %7d", arc_weights[ j ] ) ;
                ///fprintf( fp, "  %7d", edml_drt_get_arc_weight( drt, inode, j ) ) ;
            }
        }
        fprintf( fp, "\n" ) ;
    }
    ///
    fprintf( fp, "\n *** Lower Arc Table\n            lower arcs\n    n       e      u      d      f\n" ) ;
    for ( inode = drt->nnode - 1 ; inode >= 0 ; inode-- ) {
        int j ;
        int number_lower_arcs = drt->nodes[ inode ]->number_lower_arcs ;
        int *lower_arcs       = drt->nodes[ inode ]->lower_arcs ;
        ///
        fprintf( fp, "%5d ", inode ) ;
        for ( j = 0 ; j < 4 ; j++ ) {
            if ( lower_arcs[ j ] == -1 ) {
                fprintf( fp, "      -" ) ;
            } else {
                fprintf( fp, "  %5d", lower_arcs[ j ] ) ;
            }
        }
        fprintf( fp, "  %7d \n", number_lower_arcs ) ;
    }
    fprintf( fp, "End generate DRT.\n" ) ;
    return 0 ;
}


//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int edml_drt_free( edml_drt_t drt )
{
    if ( drt != NULL ) {
        edml_drt_info_free( drt->info ) ;
        edml_drt_mo_free( drt->mo ) ;
        ///edml_drt_nodes_free( drt->nnode, drt->nodes ) ;
        ///free( drt->nodes ) ;
        if ( drt->csfs != NULL ) {
            int i ;
            for ( i = 0 ; i < drt->ncsf ; i++ ) {
                if ( drt->csfs[ i ] != NULL ) {
                    free( drt->csfs[ i ] ) ;
                }
            }
            free( drt->csfs ) ;
        }
        free( drt ) ;
    }
    return 0 ;
}
