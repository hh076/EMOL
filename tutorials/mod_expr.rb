
#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require 'emol_drt_soci'
require 'emol_molinfo'
require 'emol_expression'
require 'emol_futil'

load "Data_molinfo"
molinfo = Data_molinfo::molinfo_data

nel            = molinfo.nel
spin           = molinfo.spin
n_frozen       = molinfo.n_frozen
n_core         = molinfo.n_core
n_active       = molinfo.n_active
n_external     = molinfo.n_external
nve            = nel - 2 * n_frozen
flg_use_symbol = molinfo.flg_use_symbol

printf( "nve, ncore, nactive, nexternal: %4d, %4d, %4d, %4d\n", nve, n_core, n_active, n_external )

drt         = Second_order_CI.new( nve, spin, n_core, n_active, n_external )
br          = BrooksCases.new( drt )
expr_actual = br.expr
expr_one    = br.expr_one

#####################################################################
modulename_expr = "Data_expression"
filename_expr   = "Data_expression"
fwrite_multiobject( [ expr_actual ], [ "expr" ], 
                    modulename_expr, filename_expr )
######################################################################
######################################################################
modulename_expr = "Data_expression_one"
filename_expr   = modulename_expr
fwrite_multiobject( [ expr_one ], [ "expr" ], 
                    modulename_expr, filename_expr )
######################################################################
#
printf( "\nEXPRESSION : completed.\n" )
