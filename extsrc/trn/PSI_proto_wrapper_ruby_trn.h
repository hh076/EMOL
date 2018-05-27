/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@soc.ait.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_TRN_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_TRN_

#include "ruby.h"

VALUE inttrans_wrapper_init   ( VALUE self ) ;

VALUE inttrans_wrapper_trnpp  ( VALUE self,
          VALUE _nin, VALUE _nout, VALUE _np,  VALUE _mp,
          VALUE _a,   VALUE _b,    VALUE _ndd, VALUE _cp, VALUE _cpp ) ;

VALUE inttrans_wrapper_trnrr  ( VALUE self,
          VALUE _nin, VALUE _ntl, VALUE _np, VALUE _mp,
          VALUE _a,   VALUE _b,   VALUE _cp, VALUE _cpq ) ;

VALUE inttrans_wrapper_trnsps ( VALUE self,
          VALUE _nr, VALUE _nc, VALUE _b ) ;

#endif
