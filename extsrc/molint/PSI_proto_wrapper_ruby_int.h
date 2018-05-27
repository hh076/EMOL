/// ====================================================================
/// 
///      Hiroaki HONDA, Yu-ichi INADOMI, and Jun MAKI
/// 
///      All rights reserved since Aug. 21th 2008.
///      RHF Part of OpenFMO project: http://www.openfmo.org
///      Correspondence: dahon@c.csce.kyushu-u.ac.jp
/// ====================================================================

#ifndef _INCLUDE_PSI_PROTO_WRAPPER_RUBY_INT_
#define _INCLUDE_PSI_PROTO_WRAPPER_RUBY_INT_

#include "ruby.h"
#include "PSI_def_int.h"

VALUE int1_s_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _TOL,        VALUE _LEVEL_DBG_PRINT ) ;
VALUE int1_s_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh ) ;

VALUE int1_t_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _TOL,        VALUE _LEVEL_DBG_PRINT ) ;
VALUE int1_t_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh ) ;

VALUE int1_v_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _TOL,        VALUE _LEVEL_DBG_PRINT ) ;
VALUE int1_v_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh ) ;

VALUE eri0_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _TOL,        VALUE _LEVEL_DBG_PRINT ) ;
VALUE eri0_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh, VALUE _ksh, VALUE _lsh ) ;
VALUE eri0_exchage_order ( VALUE self ) ;

#endif
