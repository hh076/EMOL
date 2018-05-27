#include "ruby.h"
#include "PSI_proto_wig.h"
#include "PSI_proto_wrapper_ruby_wig.h"

//////////////////////////////////////////////////////////////////////////////////
/// CLASS molint /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void Init_Wigcoef ( )
{
    VALUE module ;
    module = rb_define_module ( "Wigcoef" ) ;
///
    rb_define_module_function ( module, "init",  wigcoef_wrapper_init,  0 ) ;
    rb_define_module_function ( module, "delta", wigcoef_wrapper_delta, 3 ) ;
    rb_define_module_function ( module, "f3l",   wigcoef_wrapper_f3l,   3 ) ;
    rb_define_module_function ( module, "wig6j", wigcoef_wrapper_wig6j, 6 ) ;
    rb_define_module_function ( module, "wig9x", wigcoef_wrapper_wig9x, 9 ) ;
}
