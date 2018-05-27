#!/usr/bin/ruby
#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end
require_relative 'emol_expression'

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

   def show
      @res.each do | x |
         x.each do | y |
            printf( "%2d", y )
         end
         printf( " ¥n" )
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
        @pair[ pq ]
    end

    def get_pqrs( pqrs )
        pq = get_pq( @pair[ pqrs ][ 0 ] )
        rs = get_pq( @pair[ pqrs ][ 1 ] )
        [ pq, rs ].flatten
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

    def symbolic2expression( symbol, drt_actual )
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

        printf( "Actual expression \n" )
        ij_actual_ary = []; pqrs_actual_ary = []; coef_actual_ary = []
        f_actual = []; g_actual = []
        mos = []
        for depth in 0..4 do
            mos.push( MultiLoops.new( depth, mo.ext_range.first, mo.ext_range.last ).get_res )
        end

        ij.each_with_index do | ij, ij_k |
            csf_pair = tri.get_pq( ij - 1 )
            f = csf_list[ csf_pair[ 0 ] - 1 ]
            g = csf_list[ csf_pair[ 1 ] - 1 ]
            depth, f_order, g_order, f_state, g_state = EmolTools::mo_order( 2, ext_range, f, g )

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

            f_inode0, f_id0 = csf_upper_arc_weight( drt_actual, f_actual, 0...mo.val_range.last, 0, 1 )
            g_inode0, g_id0 = csf_upper_arc_weight( drt_actual, g_actual, 0...mo.val_range.last, 0, 1 )
            mos[ depth ].each do | x |
                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = f_state[ k ]
                end
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = g_state[ k ]
                end
                f_inode, f_id = csf_upper_arc_weight( drt_actual, f_actual, mo.val_range.last...mo.ext_range.last, f_inode0, f_id0 )
                g_inode, g_id = csf_upper_arc_weight( drt_actual, g_actual, mo.val_range.last...mo.ext_range.last, g_inode0, g_id0 )
                ij_actual = pq( f_id, g_id )
	            for it in addr[ ij_k ]...addr[ ij_k + 1 ] do
                #   printf( "%5d %5d  %5d  %24.14f\n", @ij[ it ], ij_actual, @pqrs[ it ], @coef[ it ] )
                    if pqrs[ it ] >= n_1el_symbol
                        pq, rs = tri.get_pq( pqrs[ it ] - n_1el_symbol )
                        p, q = tri.get_pq( pq - 1 )
                        r, s = tri.get_pq( rs - 1 )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        r_actual = to_actual_mo( x, mo.val_range.last, r )
                        s_actual = to_actual_mo( x, mo.val_range.last, s )
                #       printf( "*** %5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
                        pqrs_actual = pqrs( p_actual, q_actual, r_actual, s_actual ) + n_1el - 1
                    else
                        p, q = tri.get_pq( pqrs[ it ] )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        pqrs_actual = pq( p_actual, q_actual ) - 1
                    end
                #   printf( "  %6d %10d  %13.8f\n", ij_actual, pqrs_actual, coef[ it ] )
                    ij_actual_ary.push( ij_actual ) 
                    pqrs_actual_ary.push( pqrs_actual ) 
                    coef_actual_ary.push( coef[ it ] ) 
                end
                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = "e"
                end
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = "e"
                end
            end
        end
        return Expression.new( [ norb, drt_actual.mk_csf_list, ij_actual_ary, pqrs_actual_ary, coef_actual_ary ] )
    end 

    def symbolic2x( symbol, drt_actual, h_mo, eris_mo, p1, p2 )
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

        printf( "Actual expression ; Actual Hamiltonian \n" )
        actual_ary = [ [], [], [], Array.new( csf_actual.size ) ]
        f_actual = []; g_actual = []
        mos = []
        for depth in 0..4 do
            mos.push( MultiLoops.new( depth, mo.ext_range.first, mo.ext_range.last ).get_res )
        end

        ij.each_with_index do | ij, ij_k |
            csf_pair = tri.get_pq( ij - 1 )
            f = csf_list[ csf_pair[ 0 ] - 1 ]
            g = csf_list[ csf_pair[ 1 ] - 1 ]
            depth, f_order, g_order, f_state, g_state = EmolTools::mo_order( 2, ext_range, f, g )

            for i in 0...mo.val_range.last do
                f_actual[ i ] = f[ i ]
                g_actual[ i ] = g[ i ]
            end
            for i in mo.val_range.last...mo.ext_range.last do
                f_actual[ i ] = "e"
                g_actual[ i ] = "e"
            end

            f_inode0, f_id0 = csf_upper_arc_weight( drt_actual, f_actual, 0...mo.val_range.last, 0, 1 )
            g_inode0, g_id0 = csf_upper_arc_weight( drt_actual, g_actual, 0...mo.val_range.last, 0, 1 )
            mos[ depth ].each do | x |
                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = f_state[ k ]
                end
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = g_state[ k ]
                end
                f_inode, f_id = csf_upper_arc_weight( drt_actual, f_actual, mo.val_range.last...mo.ext_range.last, f_inode0, f_id0 )
                g_inode, g_id = csf_upper_arc_weight( drt_actual, g_actual, mo.val_range.last...mo.ext_range.last, g_inode0, g_id0 )
                ij_actual = pq( f_id, g_id )
                sum = 0.0
	            for it in addr[ ij_k ]...addr[ ij_k + 1 ] do
                    if pqrs[ it ] >= n_1el_symbol
                        pq, rs = tri.get_pq( pqrs[ it ] - n_1el_symbol )
                        p, q = tri.get_pq( pq - 1 )
                        r, s = tri.get_pq( rs - 1 )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        r_actual = to_actual_mo( x, mo.val_range.last, r )
                        s_actual = to_actual_mo( x, mo.val_range.last, s )
                #       printf( "*** %5d  %5d  %5d  %5d : %5d  %5d  %5d  %5d \n", p, q, r, s, p_actual, q_actual, r_actual, s_actual )
                        pqrs_actual = pqrs( p_actual, q_actual, r_actual, s_actual ) + n_1el - 1
                    else
                        p, q = tri.get_pq( pqrs[ it ] )
                        p_actual = to_actual_mo( x, mo.val_range.last, p )
                        q_actual = to_actual_mo( x, mo.val_range.last, q )
                        pqrs_actual = pq( p_actual, q_actual ) - 1
                    end
                    sum += p1.call( actual_ary, ij_actual, pqrs_actual, coef[ it ], pqrs[ it ], n_1el_symbol, n_1el, h_mo, eris_mo )
                    #   -> expression
                    #   p1( actual_ary, ij_actual, pqrs_actual, coef, idummy, dummy, idummy, dummy )
                    #   printf( "  %6d %10d  %13.8f\n", ij_actual, pqrs_actual, coef )
                    #   actual_ary[ 0 ].push( ij_actual ) 
                    #   actual_ary[ 1 ].push( pqrs_actual ) 
                    #   actual_ary[ 2 ].push( coef ) 
                    #   
                    #   -> hmatrix
                    #   p1( dummy, idummy, dummy, pqrs_actual, coef, pqrs, n_1el_symbol, n_1el, sum )
                    #   if pqrs > n_1el_symbol
                    #      sum += coef * 2el[ pqrs_actual - n_1el ]
                    #   else
                    #      sum += coef * 1el[ pqrs_actual ]
                    #   end
                end
                p2.call( actual_ary, f_id, g_id, sum )
                #   -> expression
                #   p2( dummy, idummy, idummy, dummy )
                #   
                #   -> hmatrix
                #   p2( actual_ary, f_id, g_id, value )
                #   actual_ary[ 0 ].push( f_id )
                #   actual_ary[ 1 ].push( g_id )
                #   actual_ary[ 2 ].push( sum )

                f_order.each_with_index do | i, k |
                    f_actual[ x[ i ] - 1 ] = "e"
                end
                g_order.each_with_index do | i, k |
                    g_actual[ x[ i ] - 1 ] = "e"
                end
            end
        end

        work = []
        for k in 0...actual_ary[ 0 ].size
            work.push( [ actual_ary[ 0 ][ k ], actual_ary[ 1 ][ k ], actual_ary[ 2 ][ k ] ] )
        end
        work.sort! { |a, b| a[ 0 ] <=> b [ 0 ] }
        work.each_with_index do | w, k |
            actual_ary[ 0 ][ k ] = w[ 0 ]
            actual_ary[ 1 ][ k ] = w[ 1 ]
            actual_ary[ 2 ][ k ] = w[ 2 ]
        end

        return actual_ary
        # return [ ij_actual_ary, pqrs_actual_ary, coef_actual_ary ] ]
        # return [ index_i, index_j, value, diag_h ] ]
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

    module_function :symbolic2expression
    module_function :symbolic2x
    module_function :mo_order
    module_function :to_actual_mo
    module_function :pq
    module_function :pqrs
    module_function :csf_upper_arc_weight
end
