
require 'emol_no'
require 'emol_futil'

no = NO_Analysis.new( ["Data_molinfo", "Data_rhf", "Data_expression_one", "Data_hmatrix"] )

#####################################################################
modulename_no = "Data_no"
filename_no   = "Data_no"
fwrite_multiobject( [ no ], [ "no" ], modulename_no, filename_no )
#####################################################################

print "\n NO       : completed.\n"
