# -*- coding: utf-8 -*-
def max(i, j)
  if i >= j
    return i
  else
    return j
  end
end

class StateInfo
  attr_reader :nve, :is
  def initialize(_nve, _is)
#           _nve : 凍結殼以外の電子数 ( number of valence electrons )
#           _is  : S(スピン量子数) * 2 
    @nve = _nve
    @is = _is
  end
end

class MO_classification
   attr_reader :norb, :ncore, :nval, :next, :core_range, :val_range, :ext_range
   def initialize( _ncore, _nval, _nexternal )
#           _ncore     : 0-2 電子励起可能な軌道数
#           _nval      : 自由な個数の電子占有可能な軌道数
#           _nexternal : 0-2 電子占有可能な軌道数
      @ncore = _ncore
      @nval  = _nval
      @next  = _nexternal
      @norb  = @ncore + @nval + @next
      @core_range = 1..@ncore
      @val_range = @ncore+1..@ncore+@nval
      @ext_range = @ncore+@nval+1..@ncore+@nval+@next
   end
end

class Node
  attr_reader :nve, :is, :orb, :code, :status, :upper_arc, :number_upper_arcs, :arc_weight, :lower_arc, :number_lower_arcs

  def initialize( _orb, _nve, _is )
    @orb = _orb
    @nve = _nve
    @is = _is
    @code = [ @orb, @nve, @is ]
    @status = check_node
    @number_upper_arcs = nil
    @arc_weight = [nil, nil, nil, nil]
    @lower_arc = [nil, nil, nil, nil]
  end

  def check_node
    if ( @is < 0 || @nve > $Info.nve || @orb > $mo.norb ) || (@orb >= $mo.val_range.last && @nve <= $Info.nve - 3)
#                            external の電子数は 2 個以下
      false
    else
      true
    end
  end

  def set_status(st)
    @status = st
  end

  def set_upper_arc(upper_arc)
    @upper_arc = upper_arc.clone
  end

  def set_arc_weight(i, weight)
    @arc_weight[i] = weight
  end

  def set_number_upper_arcs(num)
    @number_upper_arcs = num
  end

  def set_lower_arc(i, from_node)
    @lower_arc[ i ] = from_node
  end

  def set_number_lower_arcs(num)
    @number_lower_arcs = num
  end

  def if_connect_node( other_node )    #  他のノードとの連結の有無
     empty = [ @orb + 1, @nve, @is ] == other_node.code
     up    = [ @orb + 1, @nve+1, @is+1 ] == other_node.code
     down  = [ @orb + 1, @nve+1, @is-1 ] == other_node.code
     full  = [ @orb + 1, @nve+2, @is ] == other_node.code
     return empty || up || down || full
  end

  def if_connect_area( nodes )    #  ノードの集合(ハッシュ)との連結の有無  
     connect = false
     nodes.each do | k, node |
        connect = connect || self.if_connect_node( node )
     end
     return connect
  end
  
  def show
    printf(" %5d  %5d  %5d     %s\n", @orb, @nve, @is, @status)
  end

  def show_upper_arc
    @upper_arc.each do | x |
      if x == nil
        printf("      -")
      else
        printf("  %5d", x)
      end
    end
    printf("  %7d ", @number_upper_arcs)
    @arc_weight.each do | x |
      if x == nil
        printf("        -")
      else
        printf("  %7d", x)
      end
    end
    printf("\n")
  end

  def show_lower_arc
    @lower_arc.each do | x |
      if x == nil
        printf("      -")
      else
        printf("  %5d", x)
      end
    end
    printf("  %7d ", @number_lower_arcs)
    printf("\n")
  end
end

class Second_order_CI
  def initialize( _nve, _is, _ncore, _nval, _next )
    $Info = StateInfo.new( _nve, _is)
    $mo = MO_classification.new( _ncore, _nval, _next )

    h_nodes = {}
    generate_core_node( h_nodes )
#              core_nodes の k = ncore を抜き出す。ncore = 0 では、bottom
#              node を格納
    nodes_core_last = {}
    h_nodes.each do | k, node |
#    h_nodes.sort.each do | k, node |
#      printf(" %s ", k) 
#      node.show
      if k[0] == $mo.core_range.last
        nodes_core_last.store( k, node )
      end
    end
#    p "p h_nodes"
#    p h_nodes
 
    generate_external_node( h_nodes )
#              ext_nodes の k = ncore + nval + 1 を抜き出す。next = 0 では、top
#              node を格納
    nodes_ext_first = {}
    h_nodes.each do | k, node |
      if $mo.next > 0 && k[0] == $mo.ext_range.first
        nodes_ext_first.store( k, node )
      elsif $mo.next == 0 && k[0] == $mo.ext_range.last
        nodes_ext_first.store( k, node )
      end
    end

#              core_nodes の最終断面から出発して ext_nodes の第一断面に到達する valence 領域の
#              nodes の生成
    nodes_core_last.each do | k, node |
      generate_valence_node(node.orb+1, node.nve, node.is, nodes_ext_first, h_nodes)
      generate_valence_node(node.orb+1, node.nve+1, node.is+1, nodes_ext_first, h_nodes)
      generate_valence_node(node.orb+1, node.nve+1, node.is-1, nodes_ext_first, h_nodes)
      generate_valence_node(node.orb+1, node.nve+2, node.is, nodes_ext_first, h_nodes)
    end

#              h_nodes ( unordered hash ) => nodes ( ordered array )
#              生成した全ノードから不正ノードを消去したうえで、整列化して配列 @nodes を生成
    h_nodes.delete_if {|x, v| v.status == false}
    tmp = h_nodes.values
    @nodes = tmp.sort{|a,b| a.code <=> b.code} 

#   各ノードの upper_arc の設定
#              { [ orb, nel, is ] => number } ハッシュ作成
    code2num = {}
    @nodes.each_with_index do | node, i |
      code2num.store( node.code, i )
    end
    @nodes.each_with_index do | node, inode |
      empty = code2num[ [node.orb + 1, node.nve, node.is] ]
      up    = code2num[ [node.orb + 1, node.nve + 1, node.is + 1] ]
      down  = code2num[ [node.orb + 1, node.nve + 1, node.is - 1] ]
      full  = code2num[ [node.orb + 1, node.nve + 2, node.is] ]
      node.set_upper_arc( [empty, up, down, full] )

      [empty, up, down, full].each_with_index do | jnode, j |
        if jnode != nil
          @nodes[ jnode ].set_lower_arc( j, inode ) 
        end 
      end
    end
#            各ノードの number_upper_arcs の設定, 再帰的処理
    num_upper_arcs = store_number_upper_arcs( 0 )
    num_lower_arcs = store_number_lower_arcs( @nodes.size - 1 )
    set_lower_path_weights
    set_upper_path_weights
#
    show_drt
  end

  def show_drt_title
    print " \n  D I S T I N C T  R O W  T A B L E --- S O C I                                      coded by T. Noro \n \n"
  end

  def show_drt
    show_drt_title
    printf( "   STATE information   NVE = %3d,  S*2 = %3d \n", $Info.nve, $Info.is )
    printf( "   MO classification   NORB = %3d", $mo.norb )
    printf( "   ( NCORE=%3d, NVAL=%3d, NEXT=%3d )\n\n", $mo.ncore, $mo.nval, $mo.next )

    printf( "   Number of generated nodes = %7d \n", get_nnodes )
    printf( "   Number of generated csfs  = %7d \n", get_ncsfs )

    print "\n *** Nodes, Upper Arc Table, and Arc Weights\n"
    print "                             upper arcs                #arcs         upper weight\n"
    print "     n    mo  el  is         e      u      d      f                  e        u        d        f\n"
    (0..@nodes.size-1).reverse_each do | node |
      c = @nodes[ node ].code
      printf("%6d ( %3d %3d %3d ) ", node, c[0], c[1], c[2], c[3] )
      @nodes[ node ].show_upper_arc
    end
#
    print "\n *** Lower Arc Table\n"
    print "            lower arcs\n"
    print "    n       e      u      d      f\n"
     (0..@nodes.size-1).reverse_each do | node |
      printf("%5d ", node)
      @nodes[ node ].show_lower_arc
    end
  end

  def get_nnodes
    @nodes.size
  end

  def nnode
    get_nnodes
  end

  def get_ncsfs
    @nodes[ 0 ].number_upper_arcs
  end

  def r_ncsf
    get_ncsfs
  end

  def get_mo
    $mo
  end

  def get_norb
    $mo.norb  
  end

  def norb
    get_norb
  end

  def get_orb( inode )
    @nodes[ inode ].orb
  end
  
#  r_upperlevel( inode )  ... orb# Honda
  def r_upperlevel( inode )
    get_orb( inode )
  end

  def get_nve( inode )
    @nodes[ inode ].nve
  end

  def get_spin( inode )
    @nodes[ inode ].is + 1 # 2S + 1
  end
#  r_upperspin( inode ) ... spin ( 2S + 1 ) Honda
  def r_upperspin( inode )
    get_spin( inode )
  end

  def get_upper_arc( inode, walk )
    @nodes[ inode ].upper_arc[ walk ]
  end
#  r_uppercon( inode, walk ) ... next node# Honda
  def r_uppercon( inode, walk )
    get_upper_arc( inode, walk )
  end

  def get_arc_weight( inode, walk )
    @nodes[ inode ].arc_weight[ walk ]
  end
#  r_upperarcw( inode, walk ) ... weight of next node#
  def r_upperarcw( inode, walk )
    get_arc_weight( inode, walk )
  end
#
#  Total spin = S, 電子数 = Nve
#  core 領域の占有数 : Nve, Nve-1, nve-2 ( 0, 1, 2 holes )
#
#  core 領域の nodes ( orb#, nve, 2S )
#  orb# 0 : ( 0, 0, 0 ) ... bottom
#           empty        up           down  
#       1 : ( 1, 0, 0 ), ( 1, 1, 1 ), ( 1, 2, 0 )
#       k >= 2
#         : ( k, 2k - 2, 0 ), ( k, 2k - 2, 0 )  
#         : ( k, 2k - 1, 1 ),
#         : ( k, 2k, 0 )
#  説明  k ( >= 2)-th 軌道までの電子数は、既に2電子励起した場合、
#       既に1電子励起した場合、励起がない場合があり、それぞれの場合の
#       電子数は 2k-2、2k-1、2k である。スピン量子数は、2電子励起では
#       一重項(閉殼、開殼 S = 0)と三重項(開殼 S = 1)、1電子励起では
#       二重項(S = 1/2 ) そして励起のない場合は一重項( S = 0 ) の
#       node ができる。
#
  def generate_core_node( h_nodes )
    h_nodes.store( [0, 0, 0], Node.new( 0, 0, 0 ) )  # bottom node
    for k in $mo.core_range do
      if k == 1
        h_nodes.store( [1, 0, 0], Node.new( 1, 0, 0 ) )   # bottom-empty
        h_nodes.store( [1, 1, 1], Node.new( 1, 1, 1 ) )   # bottom-up
        h_nodes.store( [1, 2, 0], Node.new( 1, 2, 0 ) )   # bottom-full
      else
        nve = 2 * k
        h_nodes.store( [k, nve-2, 0], Node.new( k, nve-2, 0 ) )
        h_nodes.store( [k, nve-2, 2], Node.new( k, nve-2, 2 ) )
        h_nodes.store( [k, nve-1, 1], Node.new( k, nve-1, 1 ) )
        h_nodes.store( [k, nve, 0], Node.new( k, nve, 0 ) )
      end
    end
  end
#
#  external 領域の nodes
#       0, 1, 2 particles
#  orb# k < N-1 : ( k, ne-2, 2S+2 ), ( k, ne-2, 2S ), ( k, ne-2, 2S-2 )
#         : ( k, ne-1, 2S+1 ), ( k, ne-1, 2S-1 )
#         : ( N-1, ne, 2S )
#  orb# N-1 : ( N-1, ne-2, 2S )
#         : ( N-1, ne-1, 2S+1 ), ( N-1, ne-1, 2S-1 )
#         : ( N-1, ne, 2S )
#       N : ( N, ne, 2S )  ... top
#  説明  external の node の電子数は全電子数 ne, ne-1, ne-2 が 
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
        h_nodes.store( [k, nve-2, sp2], Node.new( k, nve - 2, sp2 ) )
        if k < $mo.ext_range.last-1
          h_nodes.store( [k, nve-2, sp2-2], Node.new( k, nve - 2, sp2 - 2 ) ) 
          h_nodes.store( [k, nve-2, sp2+2], Node.new( k, nve - 2, sp2 + 2 ) )
        end 
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
#
#  valence 領域の nodes 作成
#    external の nodes に接続する nodes を再帰的に作成
#
  def generate_valence_node(orb, nve, is, dest, h_nodes)
    new_node = Node.new(orb, nve, is)
    orb_last = $mo.val_range.last

    if h_nodes.key?(new_node.code)
      return h_nodes[new_node.code].status
    end

    if !new_node.status
      h_nodes.store(new_node.code, new_node)
      return new_node.status
    end

    if orb == orb_last
      if new_node.if_connect_area( dest )
        new_node.set_status(true)
        h_nodes.store(new_node.code, new_node)
        return new_node.status
      else
        new_node.set_status(false)
        h_nodes.store(new_node.code, new_node)
        return new_node.status
      end
    end

    empty = generate_valence_node(orb+1, nve,   is,   dest, h_nodes)
    up    = generate_valence_node(orb+1, nve+1, is+1, dest, h_nodes)
    down  = generate_valence_node(orb+1, nve+1, is-1, dest, h_nodes)
    full  = generate_valence_node(orb+1, nve+2, is,   dest, h_nodes)
    if !empty && !full && !up && !down
      new_node.set_status(false)
    else
      new_node.set_status(true)
    end
    h_nodes.store(new_node.code, new_node)
    return new_node.status
  end
#
#     各 node の number_upper_arcs, arc_weight の設定
#     再帰的処理
#     引数：node 番号、戻り値： number_upper_arcs
#     node = 0 からの戻り値 number_upper_arcs が全 csf 数
#
  def store_number_upper_arcs( node )
    if @nodes[node].number_upper_arcs == nil
      if node == (@nodes.size - 1) 
        @nodes[node].set_number_upper_arcs(0)
      else  
        num_upper_arcs = 0
        weight = 0
        @nodes[node].upper_arc.each_with_index do |v, i|
          if v != nil
            num_upper_arcs += max(store_number_upper_arcs( v ), 1)
            @nodes[node].set_arc_weight(i, weight)
            weight += @nodes[v].number_upper_arcs
          end
        end
        @nodes[node].set_number_upper_arcs(num_upper_arcs)
      end
    end
    return @nodes[node].number_upper_arcs
  end
#
#     各 node の number_lower_arcs の設定
#     再帰的処理
#     引数：node 番号、戻り値： number_lower_arcs
#     最終 node からの戻り値 number_upper_arcs が全 csf 数
#
  def store_number_lower_arcs( node )
    if @nodes[node].number_lower_arcs == nil
      if node == 0
        @nodes[node].set_number_lower_arcs(0)
      else  
        num_lower_arcs = 0
        @nodes[node].lower_arc.each_with_index do |v, i|
          if v != nil
            num_lower_arcs += max(store_number_lower_arcs( v ), 1)
          end
        end
        @nodes[node].set_number_lower_arcs(num_lower_arcs)
      end
    end
    return @nodes[node].number_lower_arcs
  end
#
#     bottom から各 node にいたる全ての paths の weights(arc weights の和) の配列
#     再帰的処理
#     引数：none、戻り値：none
#
  def set_lower_path_weights
    @weight = 0
    @lower_path_weights = []
    @nodes.size.times do | inode |
      @lower_path_weights.push( [] )
    end
    @lower_path_weights[0].push( 0 )
    tree_search_lower_path( 0 )
  end

  def tree_search_lower_path( inode )
    if inode >= @nodes.size - 1
      return
    else
      [0, 1, 2, 3].each do |step|
        next_node = @nodes[ inode ].upper_arc[ step ]
        if next_node != nil
          @weight += @nodes[ inode ].arc_weight[ step ]
          @lower_path_weights[ next_node ].push( @weight )
          tree_search_lower_path( next_node )
          @weight -= @nodes[ inode ].arc_weight[ step ]
        end
      end
    end
  end 

  def get_lower_path_weights( inode )
    return @lower_path_weights[ inode ]
  end

#
#     各 node から top にいたる全ての paths の weights(arc weights の和) の配列
#     再帰的処理
#     引数：none、戻り値：none
#
  def set_upper_path_weights
    @upper_path_weights = []
    @nodes.size.times do | inode |
       arr = []
       tree_search_upper_path( inode, arr )
       @upper_path_weights.push( arr )
    end    
  end

  def tree_search_upper_path( inode, arr )
    if inode >= @nodes.size - 1
      arr.push( @weight )
      return
    else
      [0, 1, 2, 3].each do |step|
        next_node = @nodes[ inode ].upper_arc[ step ]
        if next_node != nil
          @weight += @nodes[ inode ].arc_weight[ step ]
          tree_search_upper_path( next_node, arr )
          @weight -= @nodes[ inode ].arc_weight[ step ]
        end
      end
    end
  end

  def get_upper_path_weights( inode )
    return @upper_path_weights[ inode ]
  end

#
#     全 CSF リストの作成
#     再帰的処理
#     引数：none、戻り値：none
#
  def mk_csf_list
    @walk = []
    @csf_list = []
    @csf = []
    tree_search_csf( 0 )
    @csf_list
  end

  def tree_search_csf( inode )
    wk = [ "e", "u", "d", "f" ]
    if inode >= @nodes.size - 1
      @csf_list.push( @csf.clone )
    else
      [0, 1, 2, 3].each do |step|
        next_node = @nodes[ inode ].upper_arc[ step ]
        if next_node != nil
          @walk.push( next_node )
          @csf.push( wk[ step ] )
          tree_search_csf( next_node )
          @walk.pop
          @csf.pop
        end
      end
    end
  end

  def get_ncsf
    return @csf_list.size
  end

  def get_csf( i )
    return @csf_list[ i ]
  end
end

# a = Second_order_CI.new( _nve, _is, _ncore, _nval, _next )
# soci = Second_order_CI.new(  8, 0, 2, 4, 6 )
# soci = Second_order_CI.new(  8, 0, 1, 4, 6 )
# soci = Second_order_CI.new(  8, 0, 2, 4, 6 )
# soci = Second_order_CI.new(  4, 0, 0, 4, 0 )
# soci = Second_order_CI.new( 12, 0, 4, 4, 4 )
# soci = Second_order_CI.new( 14, 0, 4, 6, 10 )

#print "\n\n Lower Path Weights \n"
#(0...soci.get_nnodes).each do | inode |
#  printf( "%3d [", inode )
#  soci.get_lower_path_weights( inode ).each do | x |
#    printf( "%3d,", x )
#  end
#  printf( " ]\n" )
#end

#print "\n\n Upper Path Weights \n"
#(0...soci.get_nnodes).each do | inode |
#  printf( "%3d [", inode )
#  soci.get_upper_path_weights( inode ).each do | x |
#    printf( "%3d,", x )
#  end
#  printf( " ]\n" )
#end

#soci.mk_csf_list
#print "\n\n CSF \n"
#(0...soci.get_ncsf).each do | i |
#  printf( "%5d ", i )
#  csf = soci.get_csf( i )
#  csf.each do | p |
#    printf( " %s", p )
#  end
#  printf( "\n" )  
#end
