
require 'emol_hmatrix'
require 'emol_futil'

hmatrix = H_matrix.new( ["Data_molinfo", "Data_trnint", "Data_expression"] )

#####################################################################
modulename_hmatrix = "Data_hmatrix"
filename_hmatrix   = "Data_hmatrix"
fwrite_multiobject( [ hmatrix ], [ "hmatrix" ], 
                    modulename_hmatrix, filename_hmatrix )
#####################################################################

print "\n HMATRIX  : completed.\n"
