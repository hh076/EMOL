#include "ruby.h"
#include "PSI_def_int.h"
#include "PSI_proto_int1.h"
#include "PSI_proto_wrapper_ruby_int.h"

//////////////////////////////////////////////////////////////////////////////////
/// int1 s INIT //////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE int1_s_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _tol_eri,    VALUE _level_dbg_print )
{
    int i, numb_atom, numb_shell, numb_prim ;
    PSI_LI *shel_lqn,  *shel_tem, *shel_atm, *shel_add ;
    PSI_LD *atom_charg,*atom_xyz, *prim_exp, *prim_coe ;
    PSI_LD tol_eri ;
    PSI_LI level_dbg_print ;
    
    /// VALUE -> C-int ( using built-in ruby converter )
    numb_atom  = NUM2INT ( _numb_atom  ) ;
    numb_shell = NUM2INT ( _numb_shell ) ;
    numb_prim  = NUM2INT ( _numb_prim  ) ;
    tol_eri    = NUM2DBL ( _tol_eri    ) ;
    level_dbg_print
               = NUM2INT ( _level_dbg_print ) ;

    shel_lqn   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_tem   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_atm   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_add   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    atom_charg = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) ) ;
    atom_xyz   = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) * 3 ) ;
    prim_exp   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    prim_coe   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    
    /// VALUE -> C-array int or double
    for ( i = 0 ; i < numb_shell ; i++ ) {
        VALUE val_lqn, val_tem, val_atm, val_add ;
#ifdef RUBY_18
        val_lqn = RARRAY ( _shel_lqn )->ptr [ i ] ;
        val_tem = RARRAY ( _shel_tem )->ptr [ i ] ;
        val_atm = RARRAY ( _shel_atm )->ptr [ i ] ;
        val_add = RARRAY ( _shel_add )->ptr [ i ] ;
#else /// RUBY_19
        val_lqn = *( RARRAY_PTR ( _shel_lqn ) + i ) ;
        val_tem = *( RARRAY_PTR ( _shel_tem ) + i ) ;
        val_atm = *( RARRAY_PTR ( _shel_atm ) + i ) ;
        val_add = *( RARRAY_PTR ( _shel_add ) + i ) ;
#endif
        shel_lqn [ i ] = NUM2INT ( val_lqn ) ;
        shel_tem [ i ] = NUM2INT ( val_tem ) ;
        shel_atm [ i ] = NUM2INT ( val_atm ) ;
        shel_add [ i ] = NUM2INT ( val_add ) ;
    }
    for ( i = 0 ; i < numb_atom ; i++ ) {
        int j ;
        VALUE val ;
#ifdef RUBY_18
        val = RARRAY ( _atom_charg )->ptr [ i ] ;
#else /// RUBY_19
        val = *( RARRAY_PTR ( _atom_charg ) + i ) ;
#endif
        atom_charg [ i ] = NUM2DBL ( val ) ;
        for ( j = 0 ; j < 3 ; j++ ) {
#ifdef RUBY_18
            val = RARRAY ( _atom_xyz )->ptr [ i * 3 + j ] ;
#else /// RUBY_19
            val = *( RARRAY_PTR ( _atom_xyz ) + i * 3 + j ) ;
#endif
            atom_xyz [ i * 3 + j ] = NUM2DBL ( val ) ;
        }
    }
    for ( i = 0 ; i < numb_prim ; i++ ) {
        VALUE val_exp, val_coe ;
#ifdef RUBY_18
        val_exp = RARRAY ( _prim_exp )->ptr [ i ] ;
        val_coe = RARRAY ( _prim_coe )->ptr [ i ] ;
#else /// RUBY_19
        val_exp = *( RARRAY_PTR ( _prim_exp ) + i ) ;
        val_coe = *( RARRAY_PTR ( _prim_coe ) + i ) ;
#endif
        prim_exp [ i ] = NUM2DBL ( val_exp ) ;
        prim_coe [ i ] = NUM2DBL ( val_coe ) ;
    }

    int1_s_init ( numb_atom, numb_shell, numb_prim,
                shel_lqn, shel_tem, shel_atm, shel_add,
                atom_charg, atom_xyz,
                prim_exp, prim_coe,
                tol_eri,  level_dbg_print ) ;
/*
    {
        int i ;
        fprintf ( stdout, "c: %5d%5d%5d\n", numb_atom, numb_shell, numb_prim ) ;
        for ( i = 0 ; i < numb_shell ; i++ ) {
            fprintf ( stdout, "c: %5d%5d%5d%5d\n", 
                      shel_lqn [ i ], shel_tem [ i ], shel_atm [ i ], shel_add [ i ] ) ;
        }
        for ( i = 0 ; i < numb_atom ; i++ ) {
            fprintf ( stdout, "c: %5d %7.2f %23.15e %23.15e %23.15e\n",
                      i, atom_charg [ i ],       atom_xyz [ i * 3 + 0 ],
                         atom_xyz [ i * 3 + 1 ], atom_xyz [ i * 3 + 2 ] ) ;
        }
        for ( i = 0 ; i < numb_prim ; i++ ) {
            fprintf ( stdout, "c: %5d %23.15e %23.15e\n",
                      i, prim_exp [ i ], prim_coe [ i ] ) ; 
        }
    }
*/

    free ( shel_lqn ) ;
    free ( shel_tem ) ;
    free ( shel_atm ) ;
    free ( shel_add ) ;
    free ( atom_charg ) ;
    free ( atom_xyz ) ;
    free ( prim_exp ) ;
    free ( prim_coe ) ;

    return Qnil ;
}


VALUE int1_s_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh )
{
    VALUE _eris, _results ;
    PSI_LI i, ish, jsh, inttype, nsize_int ;
    PSI_LD dint [ _N_STV_ ] ;

    ish = NUM2INT ( _ish ) ;
    jsh = NUM2INT ( _jsh ) ;
    int1_s_calc ( ish, jsh, &inttype, &nsize_int, dint ) ;

    _eris = rb_ary_new2 ( nsize_int ) ;
    for ( i = 0 ; i < nsize_int ; i++ ) {
        rb_ary_store ( _eris, i, rb_float_new ( dint [ i ] ) ) ;
    }

    _results = rb_ary_new2 ( 3 ) ;
    rb_ary_store ( _results, 0, INT2NUM ( inttype   ) ) ;
    rb_ary_store ( _results, 1, INT2NUM ( nsize_int ) ) ;
    rb_ary_store ( _results, 2, _eris ) ;

/////
///
///    fprintf ( stdout, "===============================================\n" ) ;
/*
    fprintf ( stdout, "c: Orb: ( %d %d ): type = %d nsize = %d\n", ish, jsh, inttype, nsize_int ) ;
    for ( i = 0 ; i < nsize_int ; i++ ) {
        fprintf ( stdout, "c: s  %6d%26.16e\n", i, dint [ i ] ) ;
    }
*/
    return _results ;
}

//////////////////////////////////////////////////////////////////////////////////
/// int1 t INIT //////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE int1_t_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _tol_eri,    VALUE _level_dbg_print )
{
    int i, numb_atom, numb_shell, numb_prim ;
    PSI_LI *shel_lqn,  *shel_tem, *shel_atm, *shel_add ;
    PSI_LD *atom_charg,*atom_xyz, *prim_exp, *prim_coe ;
    PSI_LD tol_eri ;
    PSI_LI level_dbg_print ;
    
    /// VALUE -> C-int ( using built-in ruby converter )
    numb_atom  = NUM2INT ( _numb_atom  ) ;
    numb_shell = NUM2INT ( _numb_shell ) ;
    numb_prim  = NUM2INT ( _numb_prim  ) ;
    tol_eri    = NUM2DBL ( _tol_eri    ) ;
    level_dbg_print
               = NUM2INT ( _level_dbg_print ) ;

    shel_lqn   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_tem   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_atm   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_add   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    atom_charg = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) ) ;
    atom_xyz   = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) * 3 ) ;
    prim_exp   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    prim_coe   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    
    /// VALUE -> C-array int or double
    for ( i = 0 ; i < numb_shell ; i++ ) {
        VALUE val_lqn, val_tem, val_atm, val_add ;
#ifdef RUBY_18
        val_lqn = RARRAY ( _shel_lqn )->ptr [ i ] ;
        val_tem = RARRAY ( _shel_tem )->ptr [ i ] ;
        val_atm = RARRAY ( _shel_atm )->ptr [ i ] ;
        val_add = RARRAY ( _shel_add )->ptr [ i ] ;
#else ///RUBY_18
        val_lqn = *( RARRAY_PTR ( _shel_lqn ) + i ) ;
        val_tem = *( RARRAY_PTR ( _shel_tem ) + i ) ;
        val_atm = *( RARRAY_PTR ( _shel_atm ) + i ) ;
        val_add = *( RARRAY_PTR ( _shel_add ) + i ) ;
#endif
        shel_lqn [ i ] = NUM2INT ( val_lqn ) ;
        shel_tem [ i ] = NUM2INT ( val_tem ) ;
        shel_atm [ i ] = NUM2INT ( val_atm ) ;
        shel_add [ i ] = NUM2INT ( val_add ) ;
    }
    for ( i = 0 ; i < numb_atom ; i++ ) {
        int j ;
        VALUE val ;
#ifdef RUBY_18
        val = RARRAY ( _atom_charg )->ptr [ i ] ;
#else /// RUBY_19
        val = *( RARRAY_PTR ( _atom_charg ) + i ) ;
#endif
        atom_charg [ i ] = NUM2DBL ( val ) ;
        for ( j = 0 ; j < 3 ; j++ ) {
#ifdef RUBY_18
            val = RARRAY ( _atom_xyz )->ptr [ i * 3 + j ] ;
#else /// RUBY_19
            val = *( RARRAY_PTR ( _atom_xyz ) + i * 3 + j ) ;
#endif
            atom_xyz [ i * 3 + j ] = NUM2DBL ( val ) ;
        }
    }
    for ( i = 0 ; i < numb_prim ; i++ ) {
        VALUE val_exp, val_coe ;
#ifdef RUBY_18
        val_exp = RARRAY ( _prim_exp )->ptr [ i ] ;
        val_coe = RARRAY ( _prim_coe )->ptr [ i ] ;
#else /// RUBY_19
        val_exp = *( RARRAY_PTR ( _prim_exp ) + i ) ;
        val_coe = *( RARRAY_PTR ( _prim_coe ) + i ) ;
#endif
        prim_exp [ i ] = NUM2DBL ( val_exp ) ;
        prim_coe [ i ] = NUM2DBL ( val_coe ) ;
    }

    int1_t_init ( numb_atom, numb_shell, numb_prim,
                shel_lqn, shel_tem, shel_atm, shel_add,
                atom_charg, atom_xyz,
                prim_exp, prim_coe,
                tol_eri,  level_dbg_print ) ;
/*
    {
        int i ;
        fprintf ( stdout, "c: %5d%5d%5d\n", numb_atom, numb_shell, numb_prim ) ;
        for ( i = 0 ; i < numb_shell ; i++ ) {
            fprintf ( stdout, "c: %5d%5d%5d%5d\n", 
                      shel_lqn [ i ], shel_tem [ i ], shel_atm [ i ], shel_add [ i ] ) ;
        }
        for ( i = 0 ; i < numb_atom ; i++ ) {
            fprintf ( stdout, "c: %5d %7.2f %23.15e %23.15e %23.15e\n",
                      i, atom_charg [ i ],       atom_xyz [ i * 3 + 0 ],
                         atom_xyz [ i * 3 + 1 ], atom_xyz [ i * 3 + 2 ] ) ;
        }
        for ( i = 0 ; i < numb_prim ; i++ ) {
            fprintf ( stdout, "c: %5d %23.15e %23.15e\n",
                      i, prim_exp [ i ], prim_coe [ i ] ) ; 
        }
    }
*/

    free ( shel_lqn ) ;
    free ( shel_tem ) ;
    free ( shel_atm ) ;
    free ( shel_add ) ;
    free ( atom_charg ) ;
    free ( atom_xyz ) ;
    free ( prim_exp ) ;
    free ( prim_coe ) ;

    return Qnil ;
}


VALUE int1_t_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh )
{
    VALUE _eris, _results ;
    PSI_LI i, ish, jsh, inttype, nsize_int ;
    PSI_LD dint [ _N_STV_ ] ;

    ish = NUM2INT ( _ish ) ;
    jsh = NUM2INT ( _jsh ) ;
    int1_t_calc ( ish, jsh, &inttype, &nsize_int, dint ) ;

    _eris = rb_ary_new2 ( nsize_int ) ;
    for ( i = 0 ; i < nsize_int ; i++ ) {
        rb_ary_store ( _eris, i, rb_float_new ( dint [ i ] ) ) ;
    }

    _results = rb_ary_new2 ( 3 ) ;
    rb_ary_store ( _results, 0, INT2NUM ( inttype   ) ) ;
    rb_ary_store ( _results, 1, INT2NUM ( nsize_int ) ) ;
    rb_ary_store ( _results, 2, _eris ) ;

/////
///
///    fprintf ( stdout, "===============================================\n" ) ;
///    fprintf ( stdout, "c: (%3d,%3d ):%5d%5d\n", ish, jsh, inttype, nsize_int ) ;
///    for ( i = 0 ; i < nsize_int ; i++ ) {
///        fprintf ( stdout, "c: %5d: %23.15e\n", i, dint [ i ] ) ;
///    }
    return _results ;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/// int1 v INIT //////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
VALUE int1_v_wrapper_init ( VALUE self,
          VALUE _numb_atom,  VALUE _numb_shell, VALUE _numb_prim,
          VALUE _shel_lqn,   VALUE _shel_tem,   VALUE _shel_atm,  VALUE _shel_add,
          VALUE _atom_charg, VALUE _atom_xyz,   VALUE _prim_exp,  VALUE _prim_coe,
          VALUE _tol_eri,    VALUE _level_dbg_print )
{
    int i, numb_atom, numb_shell, numb_prim ;
    PSI_LI *shel_lqn,  *shel_tem, *shel_atm, *shel_add ;
    PSI_LD *atom_charg,*atom_xyz, *prim_exp, *prim_coe ;
    PSI_LD tol_eri ;
    PSI_LI level_dbg_print ;
    
    /// VALUE -> C-int ( using built-in ruby converter )
    numb_atom  = NUM2INT ( _numb_atom  ) ;
    numb_shell = NUM2INT ( _numb_shell ) ;
    numb_prim  = NUM2INT ( _numb_prim  ) ;
    tol_eri    = NUM2DBL ( _tol_eri    ) ;
    level_dbg_print
               = NUM2INT ( _level_dbg_print ) ;

    shel_lqn   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_tem   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_atm   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    shel_add   = ( PSI_LI * ) malloc ( numb_shell * sizeof ( PSI_LI ) ) ;
    atom_charg = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) ) ;
    atom_xyz   = ( PSI_LD * ) malloc ( numb_atom * sizeof ( PSI_LD ) * 3 ) ;
    prim_exp   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    prim_coe   = ( PSI_LD * ) malloc ( numb_prim * sizeof ( PSI_LD ) ) ;
    
    /// VALUE -> C-array int or double
    for ( i = 0 ; i < numb_shell ; i++ ) {
        VALUE val_lqn, val_tem, val_atm, val_add ;
#ifdef RUBY_18
        val_lqn = RARRAY ( _shel_lqn )->ptr [ i ] ;
        val_tem = RARRAY ( _shel_tem )->ptr [ i ] ;
        val_atm = RARRAY ( _shel_atm )->ptr [ i ] ;
        val_add = RARRAY ( _shel_add )->ptr [ i ] ;
#else /// RUBY_19
        val_lqn = *( RARRAY_PTR ( _shel_lqn ) + i ) ;
        val_tem = *( RARRAY_PTR ( _shel_tem ) + i ) ;
        val_atm = *( RARRAY_PTR ( _shel_atm ) + i ) ;
        val_add = *( RARRAY_PTR ( _shel_add ) + i ) ;
#endif
        shel_lqn [ i ] = NUM2INT ( val_lqn ) ;
        shel_tem [ i ] = NUM2INT ( val_tem ) ;
        shel_atm [ i ] = NUM2INT ( val_atm ) ;
        shel_add [ i ] = NUM2INT ( val_add ) ;
    }

    for ( i = 0 ; i < numb_atom ; i++ ) {
        int j ;
        VALUE val ;
#ifdef RUBY_18
        val = RARRAY ( _atom_charg )->ptr [ i ] ;
#else /// RUBY_19
        val = *( RARRAY_PTR ( _atom_charg ) + i ) ;
#endif
        atom_charg [ i ] = NUM2DBL ( val ) ;
        for ( j = 0 ; j < 3 ; j++ ) {
#ifdef RUBY_18
            val = RARRAY ( _atom_xyz )->ptr [ i * 3 + j ] ;
#else /// RUBY_19
            val = *( RARRAY_PTR ( _atom_xyz ) + i * 3 + j ) ;
#endif
            atom_xyz [ i * 3 + j ] = NUM2DBL ( val ) ;
        }
    }
    for ( i = 0 ; i < numb_prim ; i++ ) {
        VALUE val_exp, val_coe ;
#ifdef RUBY_18
        val_exp = RARRAY ( _prim_exp )->ptr [ i ] ;
        val_coe = RARRAY ( _prim_coe )->ptr [ i ] ;
#else /// RUBY_19
        val_exp = *( RARRAY_PTR ( _prim_exp ) + i ) ;
        val_coe = *( RARRAY_PTR ( _prim_coe ) + i ) ;
#endif
        prim_exp [ i ] = NUM2DBL ( val_exp ) ;
        prim_coe [ i ] = NUM2DBL ( val_coe ) ;
    }

    int1_v_init ( numb_atom, numb_shell, numb_prim,
                shel_lqn, shel_tem, shel_atm, shel_add,
                atom_charg, atom_xyz,
                prim_exp, prim_coe,
                tol_eri,  level_dbg_print ) ;
/*
    {
        int i ;
        fprintf ( stdout, "cv: %5d%5d%5d\n", numb_atom, numb_shell, numb_prim ) ;
        for ( i = 0 ; i < numb_shell ; i++ ) {
            fprintf ( stdout, "cv: %5d%5d%5d%5d\n", 
                      shel_lqn [ i ], shel_tem [ i ], shel_atm [ i ], shel_add [ i ] ) ;
        }
        for ( i = 0 ; i < numb_atom ; i++ ) {
            fprintf ( stdout, "cv: %5d %7.2f %23.15e %23.15e %23.15e\n",
                      i, atom_charg [ i ],       atom_xyz [ i * 3 + 0 ],
                         atom_xyz [ i * 3 + 1 ], atom_xyz [ i * 3 + 2 ] ) ;
        }
        for ( i = 0 ; i < numb_prim ; i++ ) {
            fprintf ( stdout, "cv: %5d %23.15e %23.15e\n",
                      i, prim_exp [ i ], prim_coe [ i ] ) ; 
        }
    }
*/

    free ( shel_lqn ) ;
    free ( shel_tem ) ;
    free ( shel_atm ) ;
    free ( shel_add ) ;
    free ( atom_charg ) ;
    free ( atom_xyz ) ;
    free ( prim_exp ) ;
    free ( prim_coe ) ;

    return Qnil ;
}


VALUE int1_v_wrapper_calc ( VALUE self,
          VALUE _ish, VALUE _jsh )
{
    VALUE _eris, _results ;
    PSI_LI i, ish, jsh, inttype, nsize_int ;
    PSI_LD dint [ _N_STV_ ] ;

    ish = NUM2INT ( _ish ) ;
    jsh = NUM2INT ( _jsh ) ;
    int1_v_calc ( ish, jsh, &inttype, &nsize_int, dint ) ;

    _eris = rb_ary_new2 ( nsize_int ) ;
    for ( i = 0 ; i < nsize_int ; i++ ) {
        rb_ary_store ( _eris, i, rb_float_new ( dint [ i ] ) ) ;
    }

    _results = rb_ary_new2 ( 3 ) ;
    rb_ary_store ( _results, 0, INT2NUM ( inttype   ) ) ;
    rb_ary_store ( _results, 1, INT2NUM ( nsize_int ) ) ;
    rb_ary_store ( _results, 2, _eris ) ;

/////
///
///    fprintf ( stdout, "===============================================\n" ) ;
///    fprintf ( stdout, "cv: (%3d,%3d ):%5d%5d\n", ish, jsh, inttype, nsize_int ) ;
///    for ( i = 0 ; i < nsize_int ; i++ ) {
///        fprintf ( stdout, "cv: %5d%26.16e\n", i, dint [ i ] ) ;
///    }
    return _results ;
}
