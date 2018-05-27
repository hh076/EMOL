#!/usr/bin/ruby

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

#require_relative 'emol_drt_soci'
#require_relative 'emol_molinfo'
#require_relative 'emol_symbolic_expression'
#require_relative 'emol_symbol2realx'
#require_relative 'emol_hmatrix'
#require_relative 'emol_consts'
#require_relative 'emol_futil'
require 'emol_hmatrix_symbolic'

hmatrix = H_matrix_symbolic.new( ["Data_molinfo", "Data_trnint", "Data_symbolic_expr"] )

#####################################################################
modulename_data_hmatrix  = "Data_hmatrix"
filename_hmatrix  = modulename_data_hmatrix
fwrite_multiobject( [ hmatrix ], [ "hmatrix" ], 
                    modulename_data_hmatrix, filename_hmatrix )
#####################################################################

print "\n HMATRIX_SYMBOLIC  : "
