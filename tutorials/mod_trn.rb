#!/usr/bin/ruby

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require 'emol_trans'
require 'emol_futil'

trnint = Transformation.new( [ "Data_molinfo", "Data_int", "Data_rhf" ] )

modulename_trnint = "Data_trnint"
filename_trnint   = "Data_trnint"
fwrite_multiobject( [ trnint ], [ "trnint" ], modulename_trnint, filename_trnint )

print "\n TRAN     : completed.\n"
