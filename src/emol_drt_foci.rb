# -*- coding: utf-8 -*-

require 'emol_drt_soci'

class First_order_CI < Second_order_CI

  def show_drt_title
    print " \n  D I S T I N C T  R O W  T A B L E --- F O C I                                      coded by T. Noro \n \n"
  end
#
#  Total spin = S, 電子数 = Nve
#  core 領域の占有数 : Nve, Nve-1 ( 0, 1 holes )
#
#  core 領域の nodes ( orb#, nve, 2S )
#  orb# 0 : ( 0, 0, 0 ) ... bottom
#           up           full
#       1 : ( 1, 1, 1 ), ( 1, 2, 0 )
#       k >= 2
#         : ( k, 2k - 1, 1 ),
#         : ( k, 2k, 0 )
#  説明  k ( >= 2)-th 軌道までの電子数は、
#       既に1電子励起した場合、励起がない場合があり、それぞれの場合の
#       電子数は 2k-1、2k である。スピン量子数は、1電子励起では
#       二重項(S = 1/2 ) そして励起のない場合は一重項( S = 0 ) の
#       node ができる。
#
  def generate_core_node( h_nodes )
    h_nodes.store( [0, 0, 0], Node.new( 0, 0, 0 ) )  # bottom node
    for k in $mo.core_range do
      if k == 1
        h_nodes.store( [1, 1, 1], Node.new( 1, 1, 1 ) )   # bottom-up
        h_nodes.store( [1, 2, 0], Node.new( 1, 2, 0 ) )   # bottom-full
      else
        nve = 2 * k
        h_nodes.store( [k, nve-1, 1], Node.new( k, nve-1, 1 ) )
        h_nodes.store( [k, nve, 0], Node.new( k, nve, 0 ) )
      end
    end
  end
#
#  external 領域の nodes
#       0, 1 particles
#  orb# k <= N-1 : ( k, ne-1, 2S+1 ), ( k, ne-1, 2S-1 ), ( N-1, ne, 2S )
#       N : ( N, ne, 2S )  ... top
#  説明  external の node の電子数は全電子数 ne, ne-1 が 
#       許される。ne の場合のスピンは S、ne-1 の場合はさらに一電子
#       が占有する( up, down )のでスピンは S+1/2, S-1/2 となる。
#       ne-2 の場合はさらに 2 電子が占有されるので (up,up), (down,
#       down), (up, down), (down, up), (full)、スピンとしては 
#       S+1, S-1, S となる。ただし、軌道が N-1 の場合は 2 個の電子
#       が N 番目の軌道にしか占有することができないので、(full) に
#       対応してスピンは S しか許されない。
#
  def generate_external_node( h_nodes )
    nve = $Info.nve
    sp2 = $Info.is * 2
    for k in $mo.ext_range do
      if k < $mo.ext_range.last
        h_nodes.store( [k, nve-1, sp2-1], Node.new( k, nve - 1, sp2 - 1 ) )
        h_nodes.store( [k, nve-1, sp2+1], Node.new( k, nve - 1, sp2 + 1 ) )
      end
      h_nodes.store( [k, nve, sp2], Node.new( k, nve, sp2 ) )
    end
    if $mo.ext_range.last < $mo.ext_range.first # no external の場合 top
      # node を格納
      h_nodes.store( [$mo.ext_range.last, nve, sp2], Node.new( $mo.ext_range.last, nve, sp2 ) )
    end
  end

end

# a = First_order_CI.new( _nve, _is, _ncore, _nval, _next )
# foci = First_order_CI.new( 10, 0, 3, 4, 16 )

#print "\n\n Lower Path Weights \n"
#(0...foci.get_nnodes).each do | inode |
#  printf( "%3d [", inode )
#  foci.get_lower_path_weights( inode ).each do | x |
#    printf( "%3d,", x )
#  end
#  printf( " ]\n" )
#end

#print "\n\n Upper Path Weights \n"
#(0...foci.get_nnodes).each do | inode |
#  printf( "%3d [", inode )
#  foci.get_upper_path_weights( inode ).each do | x |
#    printf( "%3d,", x )
#  end
#  printf( " ]\n" )
#end

#foci.mk_csf_list
#print "\n\n CSF \n"
#(0...foci.get_ncsf).each do | i |
#  printf( "%5d ", i )
#  csf = foci.get_csf( i )
#  csf.each do | p |
#    printf( " %s", p )
#  end
#  printf( "\n" )  
#end
