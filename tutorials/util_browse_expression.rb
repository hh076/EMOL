#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require 'emol_symbol2realx'
require 'emol_molinfo'

FILE_MOLINFO  = "Data_molinfo"
FILE_EXPR     = "Data_expression"
FILE_EXPR_ONE = "Data_expression_one"

class BrowseExpression
    def initialize( mode ) # one or all
        load FILE_MOLINFO
        @molinfo = Data_molinfo::molinfo_data
        if mode == "one"
            load FILE_EXPR_ONE
            expr = Data_expression_one::expr_data
        else
            load FILE_EXPR
            expr = Data_expression::expr_data
        end

        @norb = expr.norb
        coef = expr.coef
        ij = expr.ij
        pqrs = expr.pqrs
        @csf_list = expr.csf_list

        ncsf = @csf_list.size
        @npq = @norb * ( @norb + 1 ) / 2
        npqrs = @npq * ( @npq + 1 ) / 2
        nmax = ncsf
        if @npq > nmax
            nmax = @npq
        end

        @tri = Pair.new( nmax )

        @expr_sort = []
        for k in 0...ij.size do
            @expr_sort.push( [ ij[ k ], pqrs[ k ], coef[ k ] ] )
        end
        @expr_sort.sort!
    end

    def show_all
        prev_ij = 1
        buf = []
        printf( "Energy Expressions ( norb = %3d, core:%2d,act:%2d,ext:%2d )\n", 
			@norb, @molinfo.n_core, @molinfo.n_active, @molinfo.n_external )
        printf( "*** IJ,  ICSF, JCSF \n" )
        printf( "       ( p, q ) or ( p, q, r, s ), Coef\n" )

        for k in 0...@expr_sort.size do
            ij = @expr_sort[ k ][ 0 ]
            if ij != prev_ij
                show_formula( buf, prev_ij )
                prev_ij = ij
                buf = []
            end
            buf.push( [ @expr_sort[ k ][ 1 ], @expr_sort[ k ][ 2 ] ] ) 
        end
        show_formula( buf, prev_ij )
    end

    def show_individual_csf_pair( icsf, jcsf )
        if icsf > jcsf
            target = icsf * ( icsf - 1 ) / 2 + jcsf
        else
            target = jcsf * ( jcsf - 1 ) / 2 + icsf
        end

        buf = []
        printf( "Energy Expressions ( norb = %3d, core:%2d,act:%2d,ext:%2d )\n", 
			@norb, @molinfo.n_core, @molinfo.n_active, @molinfo.n_external )
        printf( "ICSF : %s\n", @csf_list[ icsf - 1] )
        printf( "JCSF : %s\n", @csf_list[ jcsf - 1] )

        k = 0
        until @expr_sort[ k ][ 0 ] == target
            k += 1
        end

        while @expr_sort[ k ][ 0 ] == target
            buf.push( [ @expr_sort[ k ][ 1 ], @expr_sort[ k ][ 2 ] ] ) 
            k += 1
        end
        show_formula( buf, target )
    end

    def show_formula( buf, ij )
        buf.sort!
        printf( "***%8d, %5d, %5d\n", ij, @tri.get_pq( ij - 1 )[ 0 ], @tri.get_pq( ij - 1 )[ 1 ] )
        buf.each do | b |
            if b[ 0 ] < @npq
                p, q = @tri.get_pq( b[ 0 ] )
                printf( "     ( %3d, %3d )            %20.14f \n", p, q, b[ 1 ] )
            else
                pq, rs = @tri.get_pq( b[ 0 ] - @npq )
                p, q = @tri.get_pq( pq - 1 )
                r, s = @tri.get_pq( rs - 1 )
                printf( "     ( %3d, %3d, %3d, %3d )  %20.14f \n", p, q, r, s, b[ 1 ] )
            end
        end
    end
end 

#a = BrowseExpression.new( "all" )
a = BrowseExpression.new( "one" )
a.show_all
#a.show_individual_csf_pair( 2, 1 )
