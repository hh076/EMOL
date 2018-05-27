#include "ruby.h"
#include "PSI_proto_trn.h"
#include "PSI_proto_wrapper_ruby_trn.h"

//////////////////////////////////////////////////////////////////////////////////
/// CLASS molint /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void Init_Inttrans ( )
{
    VALUE module ;
    module = rb_define_module ( "Inttrans" ) ;
///
    rb_define_module_function ( module, "init",   inttrans_wrapper_init,   0 ) ;
    rb_define_module_function ( module, "trnpp",  inttrans_wrapper_trnpp,  9 ) ;
    rb_define_module_function ( module, "trnrr",  inttrans_wrapper_trnrr,  8 ) ;
    rb_define_module_function ( module, "trnsps", inttrans_wrapper_trnsps, 3 ) ;
}
