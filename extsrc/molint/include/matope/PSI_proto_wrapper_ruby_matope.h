/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      RHF Part of OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@c.csce.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_MATOPE_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_MATOPE_

#include "ruby.h"

VALUE rb_matope_init ( VALUE self, VALUE _nao ) ;
VALUE rb_matope_fin  ( VALUE self ) ;

VALUE rb_matope_scale_diag ( VALUE self,
                      VALUE _n,   VALUE _s,  VALUE _alpha ) ;
VALUE rb_matope_cholesky_dec ( VALUE self,
                        VALUE _n,   VALUE _s ) ;
VALUE rb_matope_eigen ( VALUE self,
                 VALUE _nao, VALUE _F, VALUE _V ) ;
VALUE rb_matope_diis_update ( VALUE self,
                 VALUE _nao, VALUE _S, VALUE _D, VALUE _F ) ;

#endif
