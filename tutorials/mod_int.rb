
require 'emol_molinfo'
require 'emol_int'
require 'emol_futil'

load 'Data_molinfo'
molinfo = Data_molinfo::molinfo_data

s, t, v, h, nz_g = emol_int( molinfo )

#####################################################################
modulename_int = "Data_int"
filename_int   = "Data_int"
fwrite_multiobject( [ molinfo, s, t, v, h, nz_g ], [ "molinfo", "s", "t", "v", "h", "nz_g" ], modulename_int, filename_int )

#####################################################################
printf "\n MOLINT      : completed.\n"
