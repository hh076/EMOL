
require 'emol_molinfo'
require 'emol_int'
require 'emol_rhf'
require 'emol_futil'

load 'Data_molinfo'
load 'Data_int'

molinfo  = Data_molinfo::molinfo_data
s        = Data_int::s_data
h        = Data_int::h_data
nz_g     = Data_int::nz_g_data

elec_energy, hf_energy, mo_energies, f, g, cao, dao, tm = rhf( molinfo, s, h, nz_g )

#####################################################################
modulename_rhf = "Data_rhf"
filename_rhf   = "Data_rhf"
fwrite_multiobject( [ elec_energy, hf_energy, mo_energies, f, g, cao, dao, tm ],
                    [ "elec_energy", "hf_energy", "mo_energies", "f", "g", "cao", "dao", "tm" ],
                    modulename_rhf, filename_rhf )

#####################################################################
printf "\n RHF      : completed.\n"
