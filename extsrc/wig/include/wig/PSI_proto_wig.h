// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WIG_
#define _INCLUDE_PSI_PROTO_WIG_

#include "PSI_stan.h"

PSI_LI setwig_ ( ) ;
PSI_LD delta_  ( PSI_LI *ia, PSI_LI *ib, PSI_LI *ic ) ;
PSI_LD f3l_    ( PSI_LI *ia, PSI_LI *ib, PSI_LI *ic ) ;
PSI_LD xwig6j_ ( PSI_LI *i1, PSI_LI *i2, PSI_LI *i3,
                 PSI_LI *i4, PSI_LI *i5, PSI_LI *i6 ) ;
PSI_LD xwig9x_ ( PSI_LI *ia, PSI_LI *ib, PSI_LI *ic,
                 PSI_LI *id, PSI_LI *ie, PSI_LI *iv,
                 PSI_LI *ig, PSI_LI *ih, PSI_LI *ii ) ;
PSI_LD  wig6j_ ( PSI_LI *i1, PSI_LI *i2, PSI_LI *i3,
                 PSI_LI *i4, PSI_LI *i5, PSI_LI *i6 ) ;
PSI_LD  wig9x_ ( PSI_LI *ia, PSI_LI *ib, PSI_LI *ic,
                 PSI_LI *id, PSI_LI *ie, PSI_LI *iv,
                 PSI_LI *ig, PSI_LI *ih, PSI_LI *ii ) ;

#endif
