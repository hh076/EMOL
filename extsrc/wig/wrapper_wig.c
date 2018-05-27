#include "ruby.h"
#include "PSI_stan.h"
#include "PSI_proto_wig.h"
#include "PSI_proto_wrapper_ruby_wig.h"

//////////////////////////////////////////////////////////////////////////////////
/// Wigner coefficients INIT /////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE wigcoef_wrapper_init ( VALUE self )
{
    setwig_ ( ) ;
    return Qnil ;
}

//////////////////////////////////////////////////////////////////////////////////
/// wigcoef coefficients //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE wigcoef_wrapper_delta ( VALUE self,
              VALUE _ia, VALUE _ib, VALUE _ic )
{
    PSI_LI ia, ib, ic ;
    PSI_LD d ;
    ia = NUM2INT ( _ia ) ;
    ib = NUM2INT ( _ib ) ;
    ic = NUM2INT ( _ic ) ;

    d  = delta_ ( &ia, &ib, &ic ) ;
    return rb_float_new ( d ) ;
    ///return DBL2NUM ( d ) ;
}

VALUE wigcoef_wrapper_f3l ( VALUE self,
              VALUE _ia, VALUE _ib, VALUE _ic )
{
    PSI_LI ia, ib, ic ;
    PSI_LD d ;
    ia = NUM2INT ( _ia ) ;
    ib = NUM2INT ( _ib ) ;
    ic = NUM2INT ( _ic ) ;

    d  = f3l_ ( &ia, &ib, &ic ) ;
    return rb_float_new ( d ) ;
}

VALUE wigcoef_wrapper_wig6j ( VALUE self,
              VALUE _i1, VALUE _i2, VALUE _i3,
              VALUE _i4, VALUE _i5, VALUE _i6 )
{
    PSI_LI i1, i2, i3, i4, i5, i6 ;
    PSI_LD d ;
    i1 = NUM2INT ( _i1 ) ;
    i2 = NUM2INT ( _i2 ) ;
    i3 = NUM2INT ( _i3 ) ;
    i4 = NUM2INT ( _i4 ) ;
    i5 = NUM2INT ( _i5 ) ;
    i6 = NUM2INT ( _i6 ) ;

    d  = xwig6j_ ( &i1, &i2, &i3, &i4, &i5, &i6 ) ;
    return rb_float_new ( d ) ;
    ///return DBL2NUM ( d ) ;
}

VALUE wigcoef_wrapper_wig9x ( VALUE self,
              VALUE _i1, VALUE _i2, VALUE _i3,
              VALUE _i4, VALUE _i5, VALUE _i6,
              VALUE _i7, VALUE _i8, VALUE _i9 )
{
    PSI_LI i1, i2, i3, i4, i5, i6, i7, i8, i9 ;
    PSI_LD d ;
    i1 = NUM2INT ( _i1 ) ;
    i2 = NUM2INT ( _i2 ) ;
    i3 = NUM2INT ( _i3 ) ;
    i4 = NUM2INT ( _i4 ) ;
    i5 = NUM2INT ( _i5 ) ;
    i6 = NUM2INT ( _i6 ) ;
    i7 = NUM2INT ( _i7 ) ;
    i8 = NUM2INT ( _i8 ) ;
    i9 = NUM2INT ( _i9 ) ;

    d  = xwig9x_ ( &i1, &i2, &i3, &i4, &i5, &i6, &i7, &i8, &i9 ) ;
    ///return DBL2NUM ( d ) ;
    return rb_float_new ( d ) ;
}
