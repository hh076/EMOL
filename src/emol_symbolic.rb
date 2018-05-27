#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative 'emol_drt_soci'
require_relative 'emol_molinfo'
require_relative 'emol_expression'

class PairSymbolic
    def initialize( nmax )
        @max = nmax
        @pair = Array( nmax * (nmax + 1) / 2 )   
        pq = -1
        for p in 1..nmax do
           for q in 1..p do
              pq += 1
              @pair[ pq ] = [ p, q ]
           end
        end
    end

    def get_pq( pq )
        @pair[ pq ]
    end

    def get_pqrs( pqrs )
        pq = get_pq( @pair[ pqrs ][ 0 ] )
        rs = get_pq( @pair[ pqrs ][ 1 ] )
        [ pq, rs ].flatten
    end
end

class Symbolic
    attr_accessor :norb, :ij, :addr, :pqrs, :coef, :csf_list
    def initialize( nve, spin, n_core, n_active, n_external ) 
        to_initialize( nve, spin, n_core, n_active, n_external )
    end
#    def initialize( obj = nil )
#        if ( obj.instance_of?( Array ) && obj.length == 5 ) then   #
#            nve        = obj[ 0 ]
#            spin       = obj[ 1 ]
#            n_core     = obj[ 2 ]
#            n_active   = obj[ 3 ]
#            n_external = obj[ 4 ]
#            to_initialize( nve, spin, n_core, n_active, n_external )
#        end
##       if obj.nil? then
#           @norb = 0
#           @csf_list = []
#           @ij = []
#           @addr = []
#           @pqrs = []
#           @coef = []
#       end
#
#       if ( obj.instance_of?( Molinfo ) ) then   #
#           molinfo = obj
#           nel = molinfo.nel
#           spin = molinfo.spin
#           n_frozen = molinfo.n_frozen
#           n_core = molinfo.n_core
#           n_active = molinfo.n_active
#           n_external = molinfo.n_external
#           nve = nel - 2 * n_frozen
#           to_initialize( nve, spin, n_core, n_active, n_external )
#       elsif ( obj.instance_of?( Array ) && obj.length == 5 ) then   #
#           nve        = obj[ 0 ]
#           spin       = obj[ 1 ]
#           n_core     = obj[ 2 ]
#           n_active   = obj[ 3 ]
#           n_external = obj[ 4 ]
#           to_initialize( nve, spin, n_core, n_active, n_external )
#       elsif ( obj.instance_of?( Array ) && obj.length = 6 ) then   #
#           @norb = obj[ 0 ]
#           @csf_list = obj[ 1 ].clone
#           @ij = obj[ 2 ].clone
#           @addr = obj[ 3 ].clone
#           @pqrs = obj[ 4 ].clone
#           @coef = obj[ 5 ].clone
#           #self.show
#       elsif ( obj.instance_of?( Symbolic ) ) then
#           @norb = obj.norb
#           @csf_list = obj.csf_list.clone
#           @ij = obj.ij.clone
#           @addr = obj.nterm.clone
#           @pqrs = obj.pqrs.clone
#           @coef = obj.coef.clone
#           #self.show
#       else
#       end
#	end

    def to_initialize( nve, spin, n_core, n_active, n_external )
        if n_external >= 4
            n_external_symbol = 4
        else
            n_external_symbol = n_external
        end

        #drt = Second_order_CI.new( nve, spin, n_core, n_active, n_external )
        #drt.mk_csf_list
        drt_symbol = Second_order_CI.new( nve, spin, n_core, n_active, n_external_symbol )
        drt_symbol.mk_csf_list
        br_symbol = BrooksCases.new( drt_symbol )
        @norb, @csf_list, @ij, @addr, @pqrs, @coef = to_nonredundant_symbolic_expression( br_symbol )
        #br.expr.show
        #br_symbol.expr.show
        #self.show
    end

    def both_fg_packed?( f, g, cond, range)  
        previous = "hp"
        for k in range do
            if f[ k - 1 ] == cond && g[ k - 1 ] == cond
                previous = "ef"       
            elsif previous == "ef"
                return false
            end 
        end
        return true
    end

    def to_nonredundant_symbolic_expression( br_symbol )
        drt = br_symbol.drt
        expr = br_symbol.expr
        n = expr.pqrs.size

        norb = drt.get_norb
        ncsf = drt.get_ncsf
        mo = drt.get_mo
        tri = PairSymbolic.new( ncsf )
 
        ij_ary = []
        addr_ary = []
        pqrs_ary = []
    	coef_ary = []

        ij_prev = 0 
        nterm = 0
        addr_ary.push( 0 )
        for i in 0...n do
    #   printf( "%10d%10d   %24.14e \n", expr.ij[i], expr.pqrs[i], expr.coef[i] )
            ij = expr.ij[ i ] 
                if ij > ij_prev
                    if nterm > 0
                    addr_ary.push( addr_ary.last + nterm ) 
                    nterm = 0
                end
                csf_pair = tri.get_pq( ij - 1 )
                f = drt.get_csf( csf_pair[ 0 ] - 1 )
                g = drt.get_csf( csf_pair[ 1 ] - 1 )
                representative = both_fg_packed?( f, g, "e", mo.ext_range ) 
                if representative
                    ij_ary.push( ij )
                    ij_prev = ij
                end
            end
            if representative
                nterm += 1
                pqrs_ary.push( expr.pqrs[ i ] )
                coef_ary.push( expr.coef[ i ] )
            end
        end
        addr_ary.push( addr_ary.last + nterm ) 

        return norb, drt.mk_csf_list, ij_ary, addr_ary, pqrs_ary, coef_ary
    end

    def extract_1el_expression
        ij_1el = []
        addr_1el = []
        pqrs_1el = []
        coef_1el = []
        nint_1el = @norb * ( @norb + 1 ) / 2
        i1 = -1
        i2 = -1 
        addr_1el.push( 0 )
        @ij.each_with_index do | ij, k |
            nt = 0
            for it in @addr[ k ]...@addr[ k + 1 ] do
                if @pqrs[ it ] < nint_1el
			 	    i2 += 1
                    nt += 1
                    pqrs_1el[ i2 ] = @pqrs[ it ]
                    coef_1el[ i2 ] = @coef[ it ]
                end
            end
            if nt > 0
                i1 += 1
                ij_1el[ i1 ] = ij
                addr_1el.push( addr_1el.last + nt )
            end
        end
        return Symbolic.new( [ @norb, @csf_list, ij_1el, addr_1el, pqrs_1el, coef_1el ] )
    end

    def show
        printf( " CSF List \n" )
        @csf_list.each do | x |
            p x
        end
        printf( "\nNon-Redundant Symbolic Energy Expressions \n" )
        printf( "\n        IJ    NTERM    PQRS      COEFFICIENT\n" )
        @ij.each_with_index do | ij, k |
            printf( "%10d %5d ", ij, @addr[ k + 1 ] - @addr[ k ] )
            for it  in @addr[ k ]...@addr[ k + 1 ] do 
                if it == @addr[ k ]
                    printf( "%10d  %24.14e \n", @pqrs[ it ], @coef[ it ] )
                else
                    printf( "                 %10d  %24.14e \n", @pqrs[ it ], @coef[ it ] )
                end
            end
        end
    end
end
#
#load "Data_molinfo"
#molinfo = Data_molinfo::molinfo_data
#symbol = Symbolic.new( molinfo )
#symbol_1el = symbol.extract_1el_expression
#printf( "start output\n" )
#symbol.show
#
#require_relative '../emol_futil'
######################################################################
#modulename_data_expr  = "Data_symbol"
#filename_expr  = modulename_data_expr
#fwrite_multiobject( [ symbol ], [ "symbol" ], 
#                    modulename_data_expr, filename_expr )
######################################################################
#
######################################################################
#modulename_data_expr  = "Data_symbol_1el"
#filename_expr  = modulename_data_expr
#fwrite_multiobject( [ symbol_1el ], [ "symbol_1el" ], 
#                    modulename_data_expr, filename_expr )
######################################################################
