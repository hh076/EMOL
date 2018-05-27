#include "ruby.h"
#include "PSI_def_int.h"
#include "PSI_proto_int1.h"
#include "PSI_proto_int2.h"
#include "PSI_proto_wrapper_ruby_int.h"

//////////////////////////////////////////////////////////////////////////////////
/// CLASS molint /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void Init_Molint ( )
{
    VALUE cl_sint, cl_tint, cl_vint, cl_eri0 ;
    cl_sint  = rb_define_class ( "Sint0", rb_cObject ) ;
    cl_tint  = rb_define_class ( "Tint0", rb_cObject ) ;
    cl_vint  = rb_define_class ( "Vint0", rb_cObject ) ;
    cl_eri0  = rb_define_class ( "Eri0" , rb_cObject ) ;
///
    rb_define_method ( cl_sint, "init", int1_s_wrapper_init, 13 ) ;
    rb_define_method ( cl_sint, "calc", int1_s_wrapper_calc,  2 ) ;
    rb_define_method ( cl_tint, "init", int1_t_wrapper_init, 13 ) ;
    rb_define_method ( cl_tint, "calc", int1_t_wrapper_calc,  2 ) ;
    rb_define_method ( cl_vint, "init", int1_v_wrapper_init, 13 ) ;
    rb_define_method ( cl_vint, "calc", int1_v_wrapper_calc,  2 ) ;
///
    rb_define_method ( cl_eri0, "init", eri0_wrapper_init,   13 ) ;
    rb_define_method ( cl_eri0, "calc", eri0_wrapper_calc,    4 ) ;
    rb_define_method ( cl_eri0, "exchange_order",
                                        eri0_exchage_order,   0 ) ;
}
