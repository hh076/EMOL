#include <stdio.h>
#include <stdlib.h>
#include "slib_list.h"

slib_plist_t slib_create_list()
{
    slib_plist_t plist = ( slib_plist_t ) malloc( sizeof( slib_list_t ) ) ;
///
    plist->head = NULL ;
    plist->tail = NULL ;
    plist->size = 0 ;
///
    return plist ;
}

slib_pnode_t slib_get_front_list ( slib_plist_t list )
{
    return list->head->next ;
}
 
slib_pnode_t slib_get_back_list ( slib_plist_t list )
{
    return list->tail->prev ;
}
 
int slib_push_front_list ( slib_plist_t list, slib_pobj_t obj )
{
    slib_pnode_t new_node ;
    ///fprintf( stderr, "push_front: list->head, list->tail: %16p, %16p\n", list->head, list->tail ) ;
    ///fprintf( stderr, "push_front: here 0, list, obj, obj->item: %16p, %16p, %16p\n", list, obj, obj->item ) ;
    if ( list == NULL || obj == NULL || obj->item == NULL ) {
        return 0 ;
    }
    ///fprintf( stderr, "push_front: here 1\n" ) ;
///
    new_node = ( slib_pnode_t ) malloc ( sizeof( slib_node_t ) ) ;
    if ( new_node == NULL ) {
        return 0 ;
    }
    if ( list->size == 0 ) {
      ///fprintf( stderr, "push_front: here 2\n" ) ;
        new_node->prev      = NULL ;
        new_node->next      = NULL ;
        new_node->obj->size = obj->size ;
        new_node->obj->item = obj->item ;
        ///
        list->head          = new_node ;
        list->tail          = new_node ;
      ///fprintf( stderr, "push_front: here 3: list->head, list->tail, list->head->prev, list->head->next: %16p, %16p, %16p, %16p\n", list->head, list->tail, list->head->prev, list->head->next ) ;
    } else {
      ///fprintf( stderr, "push_front: here 4\n" ) ;
        new_node->prev      = NULL ;
        new_node->next      = list->head ;
        new_node->obj->size = obj->size ;
        new_node->obj->item = obj->item ;
        ///
        list->head->prev    = new_node ;
        list->head          = new_node ;
      ///fprintf( stderr, "push_front: here 5: list->head, list->tail, list->head->prev, list->head->next: %16p, %16p, %16p, %16p\n", list->head, list->tail, list->head->prev, list->head->next ) ;
    }
    list->size++ ;
 
    return 1 ;
}

int slib_push_back_list ( slib_plist_t list, slib_pobj_t obj )
{
    slib_pnode_t new_node ;
    ///fprintf( stderr, "push_back: list->head, list->tail: %16p, %16p\n", list->head, list->tail ) ;
    ///fprintf( stderr, "push_back: here 0, list, obj, obj->item: %16p, %16p, %16p\n", list, obj, obj->item ) ;
    if ( list == NULL || obj == NULL || obj->item == NULL ) {
        return 0 ;
    }
    ///fprintf( stderr, "push_back: here 1\n" ) ;
    new_node = ( slib_pnode_t ) malloc ( sizeof( slib_node_t ) ) ;
    if ( new_node == NULL ) {
        return 0 ;
    }
    new_node->obj = ( slib_pobj_t ) malloc ( sizeof( slib_object_t ) ) ;
    if ( new_node->obj == NULL ) {
        return 0 ;
    }
    if ( list->size == 0 ) {
      ///fprintf( stderr, "push_back: here 2\n" ) ;
        new_node->prev      = NULL ;
        new_node->next      = NULL ;
        new_node->obj->size = obj->size ;
        new_node->obj->item = obj->item ;
        ///
        list->head          = new_node ;
        list->tail          = new_node ;
      ///fprintf( stderr, "push_back: here 3: list->head, list->tail, list->tail->prev, list->tail->next: %16p, %16p, %16p, %16p\n", list->head, list->tail, list->tail->prev, list->tail->next ) ;
    } else {
      ///fprintf( stderr, "push_back: here 4\n" ) ;
        new_node->prev      = list->head ;
        new_node->next      = NULL ;
        new_node->obj->size = obj->size ;
        new_node->obj->item = obj->item ;
        ///
        list->tail->next    = new_node ;
        list->tail          = new_node ;
      ///fprintf( stderr, "push_back: here 5: list->head, list->tail, list->tail->prev, list->tail->next: %16p, %16p, %16p, %16p\n", list->head, list->tail, list->tail->prev, list->tail->next ) ;
    }
    list->size++ ;
 
    return 1 ;
}

int slib_insert_list ( slib_plist_t list, slib_pnode_t it, slib_pobj_t obj )
{
    slib_pnode_t new_node ;
 
    if ( list == NULL || it == NULL || it->prev == NULL || it->prev->next == NULL) {
        return 0 ;
    }
    new_node = ( slib_pnode_t ) malloc ( sizeof( slib_node_t ) ) ;
    if ( new_node == NULL ) {
        return 0 ;
    }
    new_node->obj = ( slib_pobj_t ) malloc ( sizeof( slib_object_t ) ) ;
    if ( new_node->obj == NULL ) {
        return 0 ;
    }
    new_node->prev      = it->prev ;
    new_node->next      = it->prev->next ;
    new_node->obj->size = obj->size ;
    new_node->obj->item = obj->item ;
    it->prev->next      = new_node ;
    it->prev            = new_node ;
    list->size++ ;
 
    return 1 ;
}

slib_pnode_t slib_head_list ( slib_plist_t list )
{
    return list->head ;
}

slib_pnode_t slib_tail_list ( slib_plist_t list )
{
    return list->tail ;
}

slib_pnode_t slib_incr_list ( slib_pnode_t it )
{
    return it->next ;
}

slib_pnode_t slib_decr_list ( slib_pnode_t it )
{
    return it->prev ;
}
 
/*
int slib_erase_item_list ( slib_plist_t list, slib_pobj_t obj,
                           int (*compare) ( const void *a, const void *b ),
                           void (*release) ( slib_pnode_t ) )
{
    slib_pnode_t walk ;
 
    if ( list == NULL || obj == NULL || list->size == 0 ) {
        return 0 ;
    }
    for ( walk = list->head->next ; walk != list->tail ; walk = walk->next) {
        if ( (*compare) ( walk->item, old_item ) == 0 ) {
            walk->next->prev = walk->prev ;
            walk->prev->next = walk->next ;
            (*release) ( walk ) ;
            list->size-- ;
            return 1 ;
        }
    }
 
    return 0 ;
}
 
int slib_erase_this_list ( slib_plist_t list, slib_pnode_t it,
                           void (*release) ( slib_pnode_t ) )
{
    if ( it == NULL || list->size == 0 ) {
        return 0 ;
    }
    it->next->prev = it->prev ;
    it->prev->next = it->next ;
    (*release) ( it ) ;
    list->size-- ;
 
    return 1 ;
}
*/
 
int slib_map_list ( slib_plist_t list, void (*action) ( slib_pnode_t ) )
{
    slib_pnode_t walk ;
    if ( list == NULL ) {
        return 0 ;
    }
    for ( walk = slib_head_list( list ) ; walk ; walk = slib_incr_list( walk ) ) {
        (*action) ( walk ) ;
    }
 
    return 1 ;
}
 
int slib_revmap_list ( slib_plist_t list, void (*action) ( slib_pnode_t ) )
{
    slib_pnode_t walk ;
    if ( list == NULL ) {
        return 0 ;
    }
    for ( walk = slib_tail_list( list ) ; walk ; walk = slib_decr_list( walk ) ) {
        (*action) ( walk ) ;
    }
 
    return 1 ;
}
