#include "ruby.h"
#include "gsl/gsl_vector.h"
#include "PSI_proto_trn.h"
#include "PSI_proto_wrapper_ruby_trn.h"

typedef struct NARRAY PSI_NA ;

//////////////////////////////////////////////////////////////////////////////////
/// Wigner coefficients INIT /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE inttrans_wrapper_init ( VALUE self )
{
    return Qnil ;
}

//////////////////////////////////////////////////////////////////////////////////
/// integral transformation library //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE inttrans_wrapper_trnpp  ( VALUE self, 
          VALUE _nin, VALUE _nout, VALUE _np,  VALUE _mp,
          VALUE _a,   VALUE _b,    VALUE _ndd, VALUE _cp, VALUE _cpq )
{
    PSI_LI nin, nout, np, mp, ndd ;
    gsl_vector *a = NULL, *b = NULL, *cp = NULL, *cpq = NULL ;

    nin  = NUM2INT ( _nin  ) ;
    nout = NUM2INT ( _nout ) ;
    np   = NUM2INT ( _np   ) ;
    mp   = NUM2INT ( _mp   ) ;
    ndd  = NUM2INT ( _ndd  ) ;
    ///CHECK_VECTOR( _a   ) ;
    ///CHECK_VECTOR( _b   ) ;
    ///CHECK_VECTOR( _cp  ) ;
    ///CHECK_VECTOR( _cpq ) ;
    Data_Get_Struct ( _a,   gsl_vector, a   ) ;
    Data_Get_Struct ( _b,   gsl_vector, b   ) ;
    Data_Get_Struct ( _cp,  gsl_vector, cp  ) ;
    Data_Get_Struct ( _cpq, gsl_vector, cpq ) ;

///    {
///        int i ;
///        for ( i = 0 ; i < ( np * ( np + 1 ) ) / 2 ; i++ ) {
///            fprintf ( stdout, " a-in[%d] = %f\n", i, a->data[ i ] ) ;
///        }
///        for ( i = 0 ; i < np * np ; i++ ) {
///            fprintf ( stdout, "cp-in[%d] = %f\n", i, cp->data[ i ] ) ;
///        }
///    }

    trnpp_ ( &nin, &nout, &np,  &mp,
             a->data, b->data, &ndd, cp->data, cpq->data ) ;

    return INT2FIX ( 0 ) ;
}

VALUE inttrans_wrapper_trnrr  ( VALUE self,
          VALUE _nin, VALUE _ntl, VALUE _np, VALUE _mp, 
          VALUE _a,   VALUE _b,   VALUE _cp, VALUE _cpq )
{
    PSI_LI nin, ntl, np, mp ;
    gsl_vector *a = NULL, *b = NULL, *cp = NULL, *cpq = NULL ;

    nin  = NUM2INT ( _nin  ) ;
    ntl  = NUM2INT ( _ntl ) ;
    np   = NUM2INT ( _np   ) ;
    mp   = NUM2INT ( _mp   ) ;
    ///CHECK_VECTOR( _a   ) ;
    ///CHECK_VECTOR( _b   ) ;
    ///CHECK_VECTOR( _cp  ) ;
    ///CHECK_VECTOR( _cpq ) ;
    Data_Get_Struct ( _a,   gsl_vector, a   ) ;
    Data_Get_Struct ( _b,   gsl_vector, b   ) ;
    Data_Get_Struct ( _cp,  gsl_vector, cp  ) ;
    Data_Get_Struct ( _cpq, gsl_vector, cpq ) ;

    trnrr_ ( &nin, &ntl, &np,  &mp,
             a->data, b->data, cp->data, cpq->data ) ;

    return INT2FIX ( 0 ) ;
}

VALUE inttrans_wrapper_trnsps ( VALUE self,
          VALUE _nr, VALUE _nc, VALUE _b )
{
    PSI_LI nr, nc ;
    gsl_vector *b = NULL ;

    nr = NUM2INT ( _nr ) ;
    nc = NUM2INT ( _nc ) ;
    ///CHECK_VECTOR( _b   ) ;
    Data_Get_Struct ( _b,   gsl_vector, b   ) ;

    trnsps_ ( &nr, &nc, b->data ) ;

    return INT2FIX ( 0 ) ;
}
