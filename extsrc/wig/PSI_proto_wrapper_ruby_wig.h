/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@soc.ait.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_WIG_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_WIG_

#include "ruby.h"

VALUE wigcoef_wrapper_init   ( VALUE self ) ;

VALUE wigcoef_wrapper_delta  ( VALUE self,
          VALUE _ia, VALUE _ib, VALUE _ic ) ;

VALUE wigcoef_wrapper_f3l    ( VALUE self,
          VALUE _ia, VALUE _ib, VALUE _ic ) ;

VALUE wigcoef_wrapper_wig6j ( VALUE self,
          VALUE _i1, VALUE _i2, VALUE _i3,
          VALUE _i4, VALUE _i5, VALUE _i6 ) ;

VALUE wigcoef_wrapper_wig9x ( VALUE self,
          VALUE _ia, VALUE _ib, VALUE _ic,
          VALUE _id, VALUE _ie, VALUE _iv,
          VALUE _ig, VALUE _ih, VALUE _ii ) ;

#endif
