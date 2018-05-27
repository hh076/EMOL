#ifndef __INCLUDE_SLIB_LIST__
#define __INCLUDE_SLIB_LIST__

typedef struct slib_object {
    size_t size ;
    void   *item ;
} slib_object_t ;
typedef slib_object_t *slib_pobj_t ;

typedef struct slib_node {
    struct slib_node *prev ;
    struct slib_node *next ;
    slib_pobj_t      obj ;
} slib_node_t ;
typedef slib_node_t *slib_pnode_t ;
typedef slib_pnode_t slib_list_it_t ;

typedef struct {
    slib_pnode_t head ;
    slib_pnode_t tail ;
    size_t       size ;
} slib_list_t ;
typedef slib_list_t *slib_plist_t ;

extern slib_plist_t   slib_create_list      ( ) ;
extern slib_pnode_t   slib_get_front_list   ( slib_plist_t list ) ;
extern slib_pnode_t   slib_get_back_list    ( slib_plist_t list ) ;
extern int            slib_push_front_list  ( slib_plist_t list, slib_pobj_t pobj ) ;
extern int            slib_push_back_list   ( slib_plist_t list, slib_pobj_t pobj ) ;
extern int            slib_insert_list      ( slib_plist_t list, slib_pnode_t it, slib_pobj_t pobj ) ;
extern int            slib_erase_item_list  ( slib_plist_t list, slib_pobj_t pobj, int (*compare) ( const void *a, const void *b ), void (*release) ( slib_pnode_t ) ) ;
extern int            slib_erase_this_list  ( slib_plist_t list, slib_pnode_t it, void (*release) ( slib_pnode_t ) ) ;
extern int            slib_map_list         ( slib_plist_t list, void (*action) ( slib_pnode_t ) ) ;
extern int            slib_revmap_list      ( slib_plist_t list, void (*action) ( slib_pnode_t ) ) ;
///
extern slib_pnode_t   slib_head_list        ( slib_plist_t list ) ;
extern slib_pnode_t   slib_tail_list        ( slib_plist_t list ) ;
extern slib_pnode_t   slib_incr_list        ( slib_pnode_t it ) ;
extern slib_pnode_t   slib_decr_list        ( slib_pnode_t it ) ;

#endif
