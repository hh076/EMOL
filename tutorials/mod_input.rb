#!/usr/bin/ruby

require 'emol_molinfo'
require 'emol_futil'

if ( ARGV.length <= 0 ) then
    $stderr.printf( "error: spcify inputfile.\n" )
    $stderr.printf( "ruby %s filename_input\n", $0 )
    exit
end
inputfilename = ARGV[ 0 ]
printf "input: %s\n", inputfilename

molinfo = Molinfo.new( inputfilename )

#####################################################################
modulename_molinfo = "Data_molinfo"
filename_molinfo   = "Data_molinfo"
fwrite_multiobject( [ molinfo ], [ "molinfo" ], modulename_molinfo, filename_molinfo )

#####################################################################
printf "\n MOLINFO      : completed.\n"
