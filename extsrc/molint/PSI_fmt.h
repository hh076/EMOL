// ====================================================================
// 
//      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
// 
//      All rights reserved since Aug. 21th 2008.
//      RHF Part of OpenFMO project: http://www.openfmo.org
//      Correspondence: dahon@c.csce.kyushu-u.ac.jp
// ====================================================================

#ifndef _FMT_H_
#define _FMT_H_

typedef int fmt_int;
typedef double fmt_dbl;

extern void fmt_initialize(fmt_int);
extern void fmt(fmt_dbl*, fmt_int, fmt_dbl, fmt_dbl);
extern void fmt_finalize();

#endif
