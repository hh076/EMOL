#include<stdio.h>

int main ( )
{
    int i, j, id, nin, nout, np, mp, ndd ;
    double a[1000], b[1000], cp[10000], cpq[10000] ;
    nin  = 1 ;
    nout = 1 ;
    np   = 5 ;
    mp   = np ;
    ndd  = 1 ;

    id = 0 ;
    for ( i = 0 ; i < np ; i++ ) {
        for ( j = 0 ; j <= i ; j++ ) {
            if ( j == i ) {
                a [ id ] = i + 1 ;
            } else {
                a [ id ] = 0 ;
            }
            ++id ;
        }
    }
    for ( i = 0 ; i < np ; i++ ) {
        for ( j = 0 ; j < mp ; j++ ) {
            cp [ i * np + j ] = 0 ;
        }
        //cp [ i * np + i ] = i + 1 ;
        cp [ i * np + i ] = 1 ;
    }

    id = 0 ;
    for ( i = 0 ; i < np ; i++ ) {
        for ( j = 0 ; j <= i ; j++ ) {
            fprintf ( stdout, "  a: %5d%5d%5d: %8.2f\n", i, j, id, a [ id ] ) ;
            ++id ;
        }
    }
    for ( i = 0 ; i < np ; i++ ) {
        for ( j = 0 ; j < mp ; j++ ) {
            fprintf ( stdout, " cp: %5d%5d%5d: %8.2f\n", i, j, i * np + j, cp [ i * np + j ] ) ;
        }
    }

    trnpp_ ( &nin, &nout, &np, &mp, a, b, &ndd, cp, cpq ) ;

    id = 0 ;
    for ( i = 0 ; i < np ; i++ ) {
        for ( j = 0 ; j <= i ; j++ ) {
            fprintf ( stdout, "  b: %5d%5d%5d: %8.2f\n", i, j, id, b [ id ] ) ;
            ++id ;
        }
    }

    return 0 ;
}
