# coding: utf-8
require 'emol_nl'
require 'emol_superci'
require_relative 'superci_driver'

istate, eval, cao, sqcdf =  driver( "RESTART" )
#####################################################################
modulename = "Data_mcscf"
filename   = "Data_mcscf"
fwrite_multiobject( [ eval[ istate - 1 ], cao, sqcdf ],
                    [ "mcscf_energy", "cao", "sqcdf" ],
                    modulename, filename )
#####################################################################

load "Data_molinfo"
molinfo = Data_molinfo::molinfo_data
if sqcdf > molinfo.get_thresh_mcscf then
	printf( "\n Causion :: Number of iterations exceeds max_iteration \n\n" )
else
	printf( "\n MCSCF     : completed. \n")
end
