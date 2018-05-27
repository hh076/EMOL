#!/usr/bin/ruby

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative 'emol_drt_soci'
require_relative 'emol_molinfo'
require_relative 'emol_symbolic_expression'
require_relative 'emol_symbol2realx'
require_relative 'emol_hmatrix'
require_relative 'emol_consts'
require_relative 'emol_futil'

load "Data_molinfo"
load "Data_symbolic_expr"

class H_matrix_symbolic < H_matrix
    def initialize( obj = nil )
        if ( obj.instance_of?( Array ) ) then
            load obj[ 0 ]
            molinfo = Data_molinfo::molinfo_data

            load obj[ 1 ]
            molint = Data_trnint::trnint_data
            core_energy = molint.core_energy
            h_mo = molint.h_mo
            eris_mo = molint.eris_mo

            load obj[ 2 ]
            expr_symbol = Data_symbolic_expr::symbolic_expr_data

            nel = molinfo.nel
            spin = molinfo.spin
            n_core = molinfo.n_core
            n_active = molinfo.n_active
            n_external = molinfo.n_external
            nve = nel - 2 * molinfo.n_frozen
            actual_drt = Second_order_CI.new( nve, spin, n_core, n_active, n_external )
            @n = actual_drt.get_ncsfs

            mkhmat( expr_symbol, actual_drt, h_mo, eris_mo )
            nstate = molinfo.get_nstate
#           thresh = EmolConsts::THR_CONV_WF
	    thresh = molinfo.get_thresh_ci
            eigen_val, @eigen_vector = Liu.new( self, nstate, thresh ).solve

            @eigen_value = GSL::Vector.alloc( nstate )
            (0...nstate).each do | k |
                @eigen_value[ k ] = eigen_val[ k ] + core_energy
            end
            printf( "\n\n" )
            printf( "Eigenvalues\n" )
            (0...nstate).each do | k |
                printf( "  %5d           ", k + 1 )
            end
            printf( "\n" )

            (0...nstate).each do | k |
                printf( "  %14.10f ", @eigen_value[ k ] )
            end
            printf( "\n\n" )
        end
    end

    def mkhmat( expr_symbol, actual_drt, h_mo, eris_mo )
        p1 = Proc.new{ |dummy, idummy, pqrs_actual, coef, pqrs_symbol, n_1el_symbol, n_1el_actual, h_mo, eris_mo|
               #printf( ">>> sum  %25.14f \n", sum )
            if pqrs_actual >= n_1el_actual
                  #printf( " 2el %24.14f     %24.14f\n", coef,   eris_mo[ pqrs_actual - n_1el_actual ])
                  # printf( " 2el %5d\n", pqrs_actual - n_1el_actual )
		    	cint = coef * eris_mo[ pqrs_actual - n_1el_actual ]
            else
                  #printf( " 1el %5d\n", pqrs_actual )
                  # printf( " 1el %24.14f     %24.14f\n", coef,   h_mo[ pqrs_actual ])
                cint = coef * h_mo[ pqrs_actual ]
            end
            cint
        }
        p2 = Proc.new{ | actual_ary, f_id, g_id, value |
#               printf( " p2 %5d   %5d   %24.14f\n", f_id,  g_id, value)
           if f_id == g_id
                    actual_ary[ 3 ][ f_id - 1 ] = value
                elsif value.abs > EmolConsts::THR_ZERO
                    actual_ary[ 0 ].push( f_id )
                    actual_ary[ 1 ].push( g_id )
                    actual_ary[ 2 ].push( value )
                end
        }

        actual_ary = EmolTools::symbolic2x( expr_symbol, actual_drt, h_mo, eris_mo, p1, p2 )
        @h_diag = actual_ary[ 3 ].clone
        @index_i = actual_ary[ 0 ].clone
        @index_j = actual_ary[ 1 ].clone
        @value = actual_ary[ 2 ].clone
#       printf( "Diagonal elements of Hamiltonian \n" )
#       @h_diag.each_with_index do | x, i | 
#           printf( " %5d     %24.14f\n", i + 1, x )
#       end
#       printf( "\n" )
#       printf( "Off-diagonal elements of Hamiltonian : i, j, value \n" )
#       @index_i.each_with_index do | x, i |
#           printf( " %5d   %5d    %24.14f\n", x, @index_j[ i ], @value[ i ] )
#       end
    end
end
