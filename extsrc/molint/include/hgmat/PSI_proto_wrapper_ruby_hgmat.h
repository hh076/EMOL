/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      RHF Part of OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@c.csce.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_HGMAT_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_HGMAT_

#include "ruby.h"
#include "PSI_def_hgmat.h"

VALUE hgmat_wrapper_init_int ( VALUE self,
          VALUE _numb_atom,   VALUE _numb_shell,  VALUE _numb_prim,
          VALUE _shel_lqn,    VALUE _shel_tem,    VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg,  VALUE _atom_xyz,    VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _thr_ovch,    VALUE _level_print ) ;

VALUE hgmat_wrapper_init_hgmat ( VALUE self,
          VALUE _nshell,      VALUE _numb_ao,     VALUE _nsize_fock,
          VALUE _shel_lqn,    VALUE _shel_atm,    VALUE _shel_ini,
          VALUE _max_lij,     VALUE _ncsp,
          VALUE _csp_add,     VALUE _csp_tem,     VALUE _csp_ab,
          VALUE _thr_schwarz, VALUE _thr_eri,     VALUE _level_print ) ;

VALUE hgmat_wrapper_fin ( VALUE self ) ;

VALUE hgmat_wrapper_calc_hmat ( VALUE self ) ;
VALUE hgmat_wrapper_calc_gmat ( VALUE self, VALUE _dmat, VALUE _gmat ) ;

#endif
