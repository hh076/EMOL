#ifndef __INCLUDE_EDUMOL_DRT__
#define __INCLUDE_EDUMOL_DRT__

#define MAX_NNODE    10000
#define MAX_NWEIGHTS 1000
#define MAX_NCSF     10000000

typedef struct {
    int nelec ;
    int spin2 ;
} edml_drt_info_item ;
typedef edml_drt_info_item* edml_drt_info_t ;

typedef struct {
    int norb ;
    int ncor ;
    int nval ;
    int next ;
    int cor_range_first ;
    int cor_range_last  ;
    int val_range_first ;
    int val_range_last  ;
    int ext_range_first ;
    int ext_range_last  ;
} edml_drt_mo_item ;
typedef edml_drt_mo_item* edml_drt_mo_t ;

typedef struct {
    int orb ;
    int nelec ;
    int spin2 ;
    int status ;
    int number_upper_arcs ;
    int number_lower_arcs ;
    int n_upper_path_weights  ;
    int n_lower_path_weights  ;
    int *codes ;
    int *arc_weights ;
    int *upper_arcs  ;
    int *lower_arcs  ;
    int *upper_path_weights  ;
    int *lower_path_weights  ;
    int check_node ;
} edml_drt_node_item ;
typedef edml_drt_node_item* edml_drt_node_t ;

typedef struct {
    int             nnode ;
    int             ncsf ;
    int             **csfs ;
    edml_drt_info_t info ;
    edml_drt_mo_t   mo ;
    edml_drt_node_t *nodes ;
} edml_drt_item ;
typedef edml_drt_item *edml_drt_t ;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

extern edml_drt_t edml_drt_create( int nelec, int spin2, int ncor, int nactive, int nexternal ) ;
extern int edml_drt_gen_csfs( edml_drt_t drt ) ;
extern int edml_drt_show( edml_drt_t drt, FILE *fp ) ;
extern int edml_drt_free( edml_drt_t drt ) ;
extern int edml_drt_get_lower_path_weights( edml_drt_t drt, int inode ) ;
extern int edml_drt_get_upper_path_weights( edml_drt_t drt, int inode ) ;
extern int edml_drt_get_nnode( edml_drt_t drt ) ;
extern int edml_drt_get_ncsf( edml_drt_t drt ) ;
extern int edml_drt_get_norb( edml_drt_t drt ) ;
extern int edml_drt_get_csf( edml_drt_t drt, int icsf, int *csf ) ;
extern int edml_drt_show_csfs( edml_drt_t drt, FILE *fp ) ;
extern int edml_drt_get_upper_level( edml_drt_t drt, int inode ) ;
extern int edml_drt_get_upper_arc( edml_drt_t drt, int inode, int walk ) ;
extern int edml_drt_get_upper_arcweight( edml_drt_t drt, int inode, int walk ) ;
extern int edml_drt_get_upper_spin( edml_drt_t drt, int inode ) ;

#endif
