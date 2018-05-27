#printf "ruby version: %s\n", RUBY_VERSION

$fout = nil ;

if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
#   str2 = "./" + str
    str2 = "../" + str
    require str2
  end
end
require_relative 'emol_expression'
require_relative 'CSymbolic'

class MultiLoops
   def initialize( depth, nmin, nmax )
      @depth = depth
      @nmin = nmin
      @nmax = nmax
      @val = Array[ depth ]
      @res = []
      self.loopfun( 0 )
   end

   def loopfun( k )
      if k == @depth
          tmp = []
          @val.each do |x|
             tmp.push( x )
          end
          @res.push( tmp )
      else
         if k == 0
            nin = @nmin
         else
            nin = @val[k - 1] + 1
         end
         for kk in nin..@nmax do
            @val[ k ] = kk
            loopfun( k + 1 )
         end
      end
   end

   def get_res
      @res
   end

   def show( os = STDOUT )
      @res.each do | x |
         x.each do | y |
            os.printf( "%2d", y )
         end
         os.printf( " ¥n" )
      end
   end
end

class Pair
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
#$fout.printf( "pq, p, q: %6d %6d %6d\n", pq, @pair[ pq ][ 0 ], @pair[ pq ][ 1 ] )
        @pair[ pq ]
    end

    def get_pqrs( pqrs )
        pq = get_pq( @pair[ pqrs ][ 0 ] )
        rs = get_pq( @pair[ pqrs ][ 1 ] )
        [ pq, rs ].flatten
    end

    def show()
        pq = 0
        for p in 1..nmax do
           for q in 1..p do
             printf( "p,q,pair: %6d%6d: %10d\n", p, q, pair[ pq ] )
             pq += 1
           end
        end

    end
end

module EmolTools
    def mo_order( mode, range, f, g )
        wk = Array[ 4 ]
        if mode == 1
            cond = "f"
        else
            cond = "e"
        end
        f_order = []; f_state = []
        g_order = []; g_state = []

        i = -1
#$fout.printf( "%s\n", range ) ;
        for k in range
            if f[ k - 1 ] != cond
                i += 1
                f_order.push( i ); f_state.push( f[ k - 1 ])
                if g[ k - 1 ] != cond
                    g_order.push( i ); g_state.push( g[ k - 1 ])
                end
            elsif
                if g[ k - 1 ] != cond
                    i += 1
                    g_order.push( i ); g_state.push( g[ k - 1 ] )
                end
		    end
        end
        return i+1, f_order, g_order, f_state, g_state
    end

    def c_symbolic2expression( symbol, drt_actual, info_drt, _expr_array )
        norb = symbol.norb
        csf_list = symbol.csf_list
        ij = symbol.ij
        addr = symbol.addr
        pqrs = symbol.pqrs
        coef = symbol.coef
        
        n_1el_symbol = norb * ( norb + 1 ) / 2
        norb_actual = drt_actual.get_norb
        n_1el = norb_actual * ( norb_actual + 1 ) / 2
        csf_actual = drt_actual.mk_csf_list
        mo = drt_actual.get_mo
        
        ncsf = csf_list.size
        n_ext = norb - mo.val_range.last
#        tri = Pair.new( ncsf )
        ext_range = mo.val_range.last+1...mo.val_range.last+n_ext+1

#$fout.printf( "norb: %4d\n", norb )
#$fout.printf( "size: csf_list: %4d\n", ncsf )
#$fout.printf( "size: ij: %4d\n", ij.length )
#$fout.printf( "size: addr, pqrt, coef: %d %d %d\n", addr.length, pqrs.length, coef.length )
#$fout.printf( "n_1el_symbol: %4d\n", n_1el_symbol )
#$fout.printf( "norb_actual: %4d\n", norb_actual )
#$fout.printf( "n_1el: %4d\n", n_1el )
#$fout.printf( "size: csf_actual: %4d\n", csf_actual.length )
#$fout.printf( "n_ext: %4d\n", n_ext )
#$fout.printf( "mo.val_range_last:  %4d\n", mo.val_range.last )
#$fout.printf( "mo.ext_range_first: %4d\n", mo.ext_range.first )
#$fout.printf( "mo.ext_range_last:  %4d\n", mo.ext_range.last )
#for pq in 0...(ncsf * ( ncsf + 1 ) / 2) do
#    $fout.printf( "%10d %10s\n", pq, tri.get_pq( pq ) )
#end

#$fout.printf( "csf_list:\n%s\n", csf_list )
#$fout.printf( "ij:\n%s\n", ij )
#$fout.printf( "addr:\n%s\n", addr )
#$fout.printf( "pqrs:\n%s\n", pqrs )
#$fout.printf( "coef:\n%s\n", coef )
#$fout.printf( "csf_actual:\n%s\n", csf_actual )
#for i in 0...pqrs.length do
#    $fout.printf( "pqrs: %10d: %10d\n", i, pqrs[ i ] )
#end


        csymb = CSymbolic.new( norb, norb_actual, ncsf, csf_actual.length, 
                                ij.length, addr.length, pqrs.length, coef.length,
                                n_1el_symbol, n_1el, n_ext, 
                                mo.val_range.last, mo.ext_range.first, mo.ext_range.last )
        csymb.set_ij( ij )
        csymb.set_addr( addr )
        csymb.set_pqrs( pqrs )
        csymb.set_coef( coef )
        csymb.set_csf_symbol( csf_list )
        csymb.set_csf_actual( csf_actual )
        csymb.set_drt( info_drt[ 0 ], info_drt[ 1 ], info_drt[ 2 ], info_drt[ 3 ], info_drt[ 4 ] )
        csymb.actual_expression( )
        #ij_actual_array   = []
        #pqrs_actual_array = []
        #coef_actual_array = []
        expr_actual_array   = []
#STDERR.printf( "enter get_expression\n" )
        csymb.get_expression( expr_actual_array )
#STDERR.printf( "out   get_expression\n" )
        ij_actual_array   = expr_actual_array[ 0 ]
        pqrs_actual_array = expr_actual_array[ 1 ]
        coef_actual_array = expr_actual_array[ 2 ]
        #STDERR.printf( "after get_expression\n" ) ;
        _expr_array.push( norb_actual )
        #STDERR.printf( "here 1\n" ) ;
        _expr_array.push( csf_actual )
        #STDERR.printf( "here 2\n" ) ;
        _expr_array.push( ij_actual_array )
        #STDERR.printf( "here 3\n" ) ;
        _expr_array.push( pqrs_actual_array )
        #STDERR.printf( "here 4\n" ) ;
        _expr_array.push( coef_actual_array )
        #STDERR.printf( "here 5\n" ) ;

#$fout = open( "Dump_s2e_crr", "w" )
#$fout.printf( "Actual Expressions:\n" )
#for i in 0...ij_actual_array.length do
#    #$fout.printf( "i, ij, pqrs, coef: %10d: %10d %10d %23.16e\n", i, ij_actual_array[ i ], pqrs_actual_array[ i ], coef_actual_array[ i ] )
#    $fout.printf( "i, ij, pqrs, coef: %s: %s %s %s\n", i, ij_actual_array[ i ], pqrs_actual_array[ i ], coef_actual_array[ i ] )
#    $fout.flush
#end
#$fout.close
#
#        return Expression.new( [ norb_actual, csf_actual, ij_actual_array, pqrs_actual_array, coef_actual_array ] )
#       return [ norb_actual, csf_actual, ij_actual_array, pqrs_actual_array, coef_actual_array ]
        return nil
    end

    def rb_symbolic2expression( symbol, drt_actual, info_drt, _expr_array )
#$fout = open( "Dump_s2e_r", "w" )
        norb = symbol.norb
        csf_list = symbol.csf_list
        ij = symbol.ij
        addr = symbol.addr
        pqrs = symbol.pqrs
        coef = symbol.coef
        
        n_1el_symbol = norb * ( norb + 1 ) / 2
        norb_actual = drt_actual.get_mo.norb
        n_1el = norb_actual * ( norb_actual + 1 ) / 2
        csf_actual = drt_actual.mk_csf_list
        mo = drt_actual.get_mo
        
        ncsf = csf_list.size
        n_ext = norb - mo.val_range.last
        tri = Pair.new( ncsf )
        ext_range = mo.val_range.last+1...mo.val_range.last+n_ext+1

#$fout.printf( "Actual Expressions:\n" )
        ij_actual_array   = []
        pqrs_actual_array = []
        coef_actual_array = []
        f_actual = []; g_actual = []
        mos = []
        for depth in 0..4 do
            mos.push( MultiLoops.new( depth, mo.ext_range.first, mo.ext_range.last ).get_res )
        end
#for i in 0..4 do
#    #$fout.printf( "mos[ %d ].length = %4d\n", i, mos[ i ].length )
#    STDERR.printf( "mos[ %d ].length = %4d\n", i, mos[ i ].length )
#    mos[ i ].each do |x|
#        #$fout.printf( "r: %s\n", x )
#        STDERR.printf( "r: %s\n", x )
#    end
#end

        ij.each_with_index do | ij, ij_k |
            csf_pair = tri.get_pq( ij - 1 )
            f = csf_list[ csf_pair[ 0 ] - 1 ]
            g = csf_list[ csf_pair[ 1 ] - 1 ]
#$fout.printf( "ij, ij_k, f, g: %s, %s :%s, %s \n", ij, ij_k, f, g )
            depth, f_order, g_order, f_state, g_state = EmolTools::mo_order( 2, ext_range, f, g )
#$fout.printf( "depth; f_order, g_order; f_state, g_state: %6d; %s, %s; %s, %s \n", depth, f_order, g_order, f_state, g_state )

	        for it in addr[ ij_k ]...addr[ ij_k + 1 ] do
        #       printf( "******** %5d  %5d  %24.14f\n", ij, @pqrs[ it ], @coef[ it ] )
            end
            for i in 0...mo.val_range.last do
                f_actual[ i ] = f[ i ]
                g_actual[ i ] = g[ i ]
            end
            for i in mo.val_range.last...mo.ext_range.last do
                f_actual[ i ] = "e"
                g_actual[ i ] = "e"
            end
#$fout.printf( "fg_actual: original: %s, %s\n", f_actual, g_actual )

            f_inode0, f_id0 = csf_upper_arc_weight( drt_actual, f_actual, 0...mo.val_range.last, 0, 1 )
            g_inode0, g_id0 = csf_upper_arc_weight( drt_actual, g_actual, 0...mo.val_range.last, 0, 1 )
#$fout.printf( "f_inode0, f_id0, g_inode0, g_id0: %4d%4d, %4d%4d\n", f_inode0, f_id0, g_inode0, g_id0 )

            mos[ depth ].each do | x |
$fout.printf( "pnode: %s\n", x ) ;
#$fout.printf( "f_actual, before : %s\n", f_actual ) ;
                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = f_state[ k ]
#$fout.printf( "f_actual : %4d: %s, %s, %s, %s\n", k, i, x[ i ], f_state[ k ], f_actual[ x[ i ] - 1 ] ) ;
                end
#$fout.printf( "f_actual, after  : %s\n", f_actual ) ;
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = g_state[ k ]
                end
#$fout.printf( "fg_actual : %s, %s\n", f_actual, g_actual ) ;
$fout.printf( "x, fg_actual : %s => %s, %s\n", x, f_actual, g_actual ) ;
#$fout.printf( "%d...%d\n", mo.val_range.last, mo.ext_range.last ) ;
                f_inode, f_id = csf_upper_arc_weight( drt_actual, f_actual, mo.val_range.last...mo.ext_range.last, f_inode0, f_id0 )
                g_inode, g_id = csf_upper_arc_weight( drt_actual, g_actual, mo.val_range.last...mo.ext_range.last, g_inode0, g_id0 )
#$fout.printf( "f_inode, f_id, g_inode, g_id: %4d%4d, %4d%4d\n", f_inode, f_id, g_inode, g_id )
                ij_actual = pq( f_id, g_id )
#$fout.printf( "ij_actual: %s\n", ij_actual )
	            for it in addr[ ij_k ]...addr[ ij_k + 1 ] do
#printf( "%5d %5d  %5d  %24.14f\n", @ij[ it ], ij_actual, @pqrs[ it ], @coef[ it ] )
#$fout.printf( "it, @ij[it], ij_actual, @pqrs[ it ], @coef[ it ]: %s: %s %s %s %s\n", it, @ij, ij_actual, @pqrs, @coef )
                    p = q = r = s = -1
                    if pqrs[ it ] >= n_1el_symbol
                        pq, rs = tri.get_pq( pqrs[ it ] - n_1el_symbol )
                        p, q = tri.get_pq( pq - 1 )
                        r, s = tri.get_pq( rs - 1 )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        r_actual = to_actual_mo( x, mo.val_range.last, r )
                        s_actual = to_actual_mo( x, mo.val_range.last, s )
                        pqrs_actual = pqrs( p_actual, q_actual, r_actual, s_actual ) + n_1el - 1
#printf( "%5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
#$fout.printf( "pqrs, pqrs_actual:%5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
                    else
                        p, q = tri.get_pq( pqrs[ it ] )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        pqrs_actual = pq( p_actual, q_actual ) - 1
                    end
#printf( "  %6d %10d  %13.8f\n", ij_actual, pqrs_actual, coef[ it ] )
#$fout.printf( "ij_actual, pqrs_actual, coef[ it ]: %10d %10d %23.16e\n", ij_actual, pqrs_actual, coef[ it ] )
$fout.printf( "ij_actual, pqrs_actual, coef[ it ]: %10d %10d(%4d%4d%4d%4d) %23.16e\n", ij_actual, pqrs_actual, p, q, r, s,  coef[ it ] )

                    ij_actual_array.push( ij_actual ) 
                    pqrs_actual_array.push( pqrs_actual ) 
                    coef_actual_array.push( coef[ it ] ) 
                end
                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = "e"
                end
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = "e"
                end
            end
        end
#for i in 0...ij_actual_array.length do
#    $fout.printf( "i, ij, pqrs, coef: %10d: %10d %10d %23.16e\n", i, ij_actual_array[ i ], pqrs_actual_array[ i ], coef_actual_array[ i ] )
#end
$fout.close
        _expr_array.push( norb )
        _expr_array.push( csf_actual )
        _expr_array.push( ij_actual_array )
        _expr_array.push( pqrs_actual_array )
        _expr_array.push( coef_actual_array )
#       return [ norb, csf_actual, ij_actual_array, pqrs_actual_array, coef_actual_array ]
#       return Expression.new( [ norb, csf_actual, ij_actual_array, pqrs_actual_array, coef_actual_array ] )
        return nil
    end 

    def to_actual_mo( map, x0, x )
        case x
        when x0 + 1
            map[ 0 ]
        when x0 + 2
            map[ 1 ]
        when x0 + 3
            map[ 2 ]
        when x0 + 4
            map[ 3 ]
        else
            x
        end
    end

    def pq( p, q )
        if p > q
            max = p
            min = q
        else
            max = q
            min = p
        end
        max * ( max - 1 ) / 2 + min
    end

    def pqrs( p, q, r, s )
        pq( pq( p, q ), pq( r, s ) )
    end

    def csf_upper_arc_weight( drt, csf, range, _inode, _id )
        walk_num = { "e" => 0, "u" => 1, "d" => 2, "f" => 3 }
        inode = _inode
        id = _id
        for i in range do
            walk = walk_num[ csf[ i ] ]
            id += drt.get_arc_weight( inode, walk )
            inode = drt.get_upper_arc( inode, walk )
        end
        return inode, id
    end

    module_function :c_symbolic2expression
    module_function :rb_symbolic2expression
    module_function :mo_order
    module_function :to_actual_mo
    module_function :pq
    module_function :pqrs
    module_function :csf_upper_arc_weight
end

##printf "ruby version: %s\n", RUBY_VERSION
#if ( RUBY_VERSION.to_f <= 1.8 ) then
#  def require_relative( str )
##   str2 = "./" + str
#    str2 = "../" + str
#    require str2
#  end
#end
#
#require_relative './emol_drt_soci'
#require_relative './emol_molinfo'
#require_relative './emol_symbolic'
#require_relative './emol_expression'
#require_relative './emol_futil'
#
#load "Data_molinfo"
#load "Data_symbol"
#
#
#molinfo = Data_molinfo::molinfo_data
#symbol = Data_symbol::symbol_data
#
#expr = Symbolic.new( [ symbol.norb, symbol.csf_list, symbol.ij, symbol.addr, symbol.pqrs, symbol.coef ] )
#
#nel = molinfo.nel
#spin = molinfo.spin
#n_core = molinfo.n_core
#n_active = molinfo.n_active
#n_external = molinfo.n_external
#nve = nel - 2 * molinfo.n_frozen
#info_drt = [ nve, spin, n_core, n_active, n_external ]
#drt = Second_order_CI.new( nve, spin, n_core, n_active, n_external )
#actual_expr = EmolTools::symbolic2expression( expr, drt, info_drt )
##ij_array    = actual_expr[ 2 ]
##pqrs_array  = actual_expr[ 3 ]
##coef_array  = actual_expr[ 4 ]
##$fout.printf( "Actual Expressions:\n" )
##for i in 0...ij_array.length do
##    $fout.printf( "i, ij, pqrs, coef: %10d: %10d %10d %23.16e\n", i, ij_array[ i ], pqrs_array[ i ], coef_array[ i ] )
##end
#
