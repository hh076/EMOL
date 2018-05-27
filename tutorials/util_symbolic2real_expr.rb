#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

#require_relative 'emol_trans'
require 'emol_drt_soci'
require 'emol_molinfo'
require 'emol_symbolic_expression'
require 'emol_expression'
require 'emol_symbol2realx'
require 'emol_futil'
require 'emol_consts'

load "Data_molinfo"
#load "Data_symbolic_expr"
load "Data_symbolic_expr_1el"

molinfo = Data_molinfo::molinfo_data
#symbol = Data_symbolic_expr::symbolic_expr_data
symbol = Data_symbolic_expr_1el::symbol_1el_data
expr = Symbolic_expression.new( [ symbol.norb, symbol.csf_list, symbol.ij, symbol.addr, symbol.pqrs, symbol.coef ] )

nel = molinfo.nel
spin = molinfo.spin
n_core = molinfo.n_core
n_active = molinfo.n_active
n_external = molinfo.n_external
nve = nel - 2 * molinfo.n_frozen
drt = Second_order_CI.new( nve, spin, n_core, n_active, n_external )
#actual_expr = EmolTools::symbolic2expression( expr, drt )
#actual_expr.show

p1 = Proc.new{ |actual_ary, ij_actual, pqrs_actual, coef, idummy1, dummy2, idummy3, dummy4, dummy5|
                    null = 0
                    actual_ary[ 0 ].push( ij_actual ) 
                    actual_ary[ 1 ].push( pqrs_actual ) 
                    actual_ary[ 2 ].push( coef ) 
                    null
}

p2 = Proc.new{ | dummy1, idummy2, idummy3, dummy4 | 
}

load "Data_trnint"
require 'emol_trans'
trnint = Data_trnint::trnint_data
core_energy = trnint.core_energy
h_mo = trnint.h_mo
eris_mo = trnint.eris_mo
a = EmolTools::symbolic2x( expr, drt, h_mo, eris_mo, p1, p2 )
actual_expr = Expression.new( [ drt.get_mo.norb, drt.mk_csf_list, a[ 0 ], a[ 1 ], a[ 2 ] ] )
actual_expr.show

#####################################################################
#modulename_data_expr  = "Data_expression"
#filename_expr  = modulename_data_expr
#fwrite_multiobject( [ actual_expr ], [ "expr" ], 
#                    modulename_data_expr, filename_expr )
#####################################################################

#actual_expr_one = actual_expr.extract_1el_expression
actual_expr_one = actual_expr
actual_expr_one.show
#####################################################################
modulename_data_expr  = "Data_expression_one"
filename_expr  = modulename_data_expr
fwrite_multiobject( [ actual_expr_one ], [ "expr" ], 
                    modulename_data_expr, filename_expr )
#####################################################################

print "\n EXPR     : "
