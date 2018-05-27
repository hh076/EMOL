#!/usr/bin/ruby
# -*- coding: utf-8 -*-

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative './emol_expression'
require_relative './emol_hmatrix'
require_relative './emol_hmatrix_symbolic'
require_relative './emol_futil'
require_relative './emol_consts'
require_relative './emol_molinfo'

require 'gsl'

class NO_Analysis
    attr_accessor :nbf, :nob, :occ, :no

    def initialize( obj = nil )
        if obj.nil? then
            @nbf = 0
            @nob = 0
            @nstate = 0
            @occ = [] 
            @no = []
        end

        if ( obj.instance_of?( Array ) ) then   #
            load obj[ 0 ]
            molinfo = Data_molinfo::molinfo_data

            load obj[ 1 ]
            cao = Data_rhf::cao_data

            load obj[ 2 ]
            expr_one = Data_expression_one::expr_data

            load obj[ 3 ]
            hmat = Data_hmatrix::hmatrix_data

            @nbf = molinfo.nbf
            bs_name = molinfo.get_bs_name
            n_frozen = molinfo.n_frozen
            n_core = molinfo.n_core
            n_active = molinfo.n_active
            n_external = molinfo.n_external
            @nob = molinfo.nob
            evec = hmat.eigen_vector
            eval = hmat.eigen_value
            csf_list = expr_one.csf_list

            @nstate = hmat.eigen_value.size
            @occ = [] 
            @no = []
            (0...@nstate).each do | istate |
                no, occ = get_no( nbf, n_frozen, nob, expr_one, evec, cao, istate )
                print_ci_result( istate + 1, eval, evec, csf_list, n_core, n_active, n_external )
                EmolUtil::print_2dim_ary_with_value_moname( no.to_a, occ.to_a, bs_name, "Natural Orbitals and Occupation Numbers ", 8, "%12.6f" )
                @no.push( no )
                @occ.push( occ )
            end
    
        elsif ( obj.instance_of?( NO_Analysis ) ) then
            @nbf = obj.nbf
            @nob = obj.nob
            @occ = obj.occ.clone
            @no  = obj.no.clone
        else
        end
    end

    def get_no( nbf, n_frozen, nob, expr_one, evec, cao, istate ) 
        pair = Array.new( nob * ( nob + 1 ) / 2 )
        pq = -1
        for p in 0..nbf-1
            for q in 0..p
                pq += 1
                pair[ pq ] = [p, q]
            end
        end
#                     make density matrix in mo representation
        density = GSL::Matrix.alloc(nob, nob)
        i_prev = 1; j_prev = 1; ij_prev = 1; ij_end = 0
        expr_one.ij.each_with_index do |ij, k|
            num = expr_one.pqrs[ k ]
            coef = expr_one.coef[ k ]
            if ij == ij_prev
                i = i_prev; j = j_prev 
            elsif ij <= ij_end + i_prev
                i = i_prev
                j = ij - ij_end
            else
                ij_end += i_prev
                i = i_prev + 1
                j = ij - ij_end
            end
            p = pair[ num ][ 0 ]
            q = pair[ num ][ 1 ]
            density[ p, q ] += evec[ i-1, istate ] * evec[ j-1, istate ] * coef
            ij_prev = ij
            i_prev = i
            j_prev = j
        end
#  EmolUtil::print_2dim_ary_tri( density.to_a, "Density-matrix", 8, "%12.6f" )

#                     no in mo representation
        eval, evec = GSL::Eigen::symmv( density )
        GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_DESC )
#  EmolUtil::print_2dim_ary_with_value( evec.to_a, eval.to_a, "Natural Orbitals and Occupation Numbers", 8, "%12.6f" )

#                     no in ao representation
        cao_active = GSL::Matrix.alloc(nbf, nob)
        for i in 0..nob - 1
            for p in 0..nbf - 1
                cao_active[ p, i ] = cao[ p, n_frozen + i ]
            end
        end
#  EmolUtil::print_2dim_ary( cao_active.to_a, "Active MOs", 8, "%12.6f" )
        no_ao = cao_active * evec
#  EmolUtil::print_2dim_ary( no_ao.to_a, "NOs-AOs", 8, "%12.6f" )

#                     add frozen core 

        if n_frozen > 0
            no = cao.submatrix(nil, 0, n_frozen).horzcat( no_ao )
            fc = GSL::Vector.alloc( n_frozen )
            fc.set_all( 2.0 )
            occ = fc.connect( eval )
        else
            no = no_ao
            occ = eval
        end

#       no = GSL::Matrix.alloc(nbf, n_frozen + nob)
#       for i in 0..n_frozen - 1
#           for p in 0..nbf - 1
#               no[ p, i ] = cao[ p, i ]
#           end
#       end
#       for i in 0..nob - 1
#           for p in 0..nbf - 1
#               no[ p, n_frozen + i ] = no_ao[ p, i ]
#           end
#       end
#       occ = GSL::Vector.alloc( n_frozen + nob )
#       for i in 0..n_frozen-1
#           occ[ i ] = 2.0
#       end
#       for i in 0..nob - 1
#           occ[ n_frozen + i ] = eval[ i ]
#       end

        return no, occ
    end

    def print_ci_result( istate, eval, evec, csf_list, n_core, n_active, n_external )
        if istate == 1
            printf("\n %2dst CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
        elsif istate == 2
            printf("\n %2dnd CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
        elsif istate == 3
            printf("\n %2drd CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
        else
            printf("\n %2dth CI EIGENSTATE   TOTAL ENERGY =  %16.10f \n\n", istate, eval[ istate - 1 ])
        end

        printf("    CSF       COEF     ( CORE / ACTIVE / EXTERNAL )\n")
        printf("    ---    ----------   ") 
        n_core.times do | i |
            printf("--")
        end
        printf("  ")
        n_active.times do | i |
            printf("--")
        end
        printf("  ")
        n_external.times do | i |
            printf("--")
        end
        printf("  \n")

        vec_print = {}
        vec = evec.column( istate - 1 )
        vec.to_a.each_with_index do |x, i|
            if x.abs > 0.01
                vec_print.store( i,  x )
            end
        end

        vec_print.sort{ |(k1, v1), (k2, v2)| v2.abs <=> v1.abs }.each do | x |
            printf("  %5d  ", x[ 0 ] + 1)
            printf(" %10.6f   ", x[ 1 ])
            a = csf_list[ x[ 0 ] ]
 
            k = 0
            n_core.times do | i |
                printf("%2s", a[ k ] )
                k += 1
            end
            printf("  ")
            n_active.times do | i |
                printf("%2s", a[ k ] )
                k += 1
            end
            printf("  ")
            n_external.times do | i |
                printf("%2s", a[ k ] )
                k += 1
            end
            printf("  \n")
        #   csf_list[ x[ 0 ] ].each do |y|
        #       printf("%2s", y)
        #   end
        end
    end
end

#NO_Analysis.new( [ "Data_rhf", "Data_trnint", "Data_expression_one", "Data_hmatrix" ] )
