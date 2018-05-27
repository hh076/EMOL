#!/usr/bin/ruby
# -*- coding: utf-8 -*-

#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
#   str2 = "./" + str
    str2 = "../" + str
    require str2
  end
end

require_relative 'emol_drt_soci'
require_relative 'emol_const_ops'
require_relative 'Wigcoef'

#require 'profile'
require 'gsl'

#$opt_debug_print = true    # or false
$opt_debug_print = false    # or false

Wigcoef.init

class InvalidCombination      < RuntimeError ; end
class InvalidLkValue          < RuntimeError ; end
class InvalidArgumentRange    < RuntimeError ; end
class InvalidLkValueAttention < RuntimeError ; end
class InvalidBySystemError    < RuntimeError ; end

##############################################
##############################################
##############################################
#  HMLKVA                                                         
#     Function:     Calculate Lk value                            
#                                                                 
#     Parameters:          Explanation                   (I/U/O/W)
#           GF             segment walk of psi - Q part       (I)
#           FF             segment walk of phi - R part       (I)
#           HF             operator h (k)                     (I)
#           PSI            spin multiplicity of psi (k)       (I)
#           PHI            spin multiplicity of phi (k)       (I)
#           THETA          spin multiplicity of theta (k)     (I)
#           THETA_PREV     spin multiplicity of theta (k-1)   (I)
#           W              return value Lk                    (O)
#           RTNCODE        0 ---> normal return               (O)
#                          1 ---> error ( ERROR condition )
#                          2 ---> error ( ATTENTION )
#                          3 ---> *****
#                         10 ---> miscellaneous

def hmlkva( gf, ff, hf, psi, phi, theta, theta_prev )
    begin
        w = 0.0e0
    
    #   #RTNCODE = 0 ;
        rtn_code = 0
    
        #H = REPR_OP ( HF );
        #F = REPR_STATE ( FF );
        #G = REPR_STATE ( GF );
        h = REPR_OP[    hf ]
        f = REPR_STATE[ ff ]
        g = REPR_STATE[ gf ]
    
        #IF        GF = EMPTY ! GF = FULL   THEN  PSI_PREV = PSI ;
        #ELSE IF   GF = UP                  THEN  PSI_PREV = PSI - 1 ;
        #ELSE IF   GF = DOWN                THEN  PSI_PREV = PSI + 1 ;
        #
        #IF        FF = EMPTY ! FF = FULL   THEN  PHI_PREV = PHI ;
        #ELSE IF   FF = UP                  THEN  PHI_PREV = PHI - 1 ;
        #ELSE IF   FF = DOWN                THEN  PHI_PREV = PHI + 1 ;
    
        if ( gf == EMPTY || gf == FULL ) then
            psi_prev = psi
        elsif ( gf == UP ) then
            psi_prev = psi - 1
        elsif ( gf == DOWN ) then
            psi_prev = psi + 1
        end
        if ( ff == EMPTY || ff == FULL ) then
            phi_prev = phi
        elsif ( ff == UP ) then
            phi_prev = phi - 1
        elsif ( ff == DOWN ) then
            phi_prev = phi + 1
        end
    
        #IF     PSI   <= 0 ! PSI_PREV   <= 0
        #     ! PHI   <= 0 ! PHI_PREV   <= 0
        #     ! THETA <= 0 ! THETA_PREV <= 0
        #THEN DO ;
        #   W = 0 ;
    
        if ( psi   <= 0 || psi_prev   <= 0 ||
             phi   <= 0 || phi_prev   <= 0 ||
             theta <= 0 || theta_prev <= 0 ) then
            w = 0
            raise InvalidArgumentRange
        else
            #/* multiply    sqrt (theta ( k-1 ) / theta ( k ) )            */
            #/* when theta ( k-1 ) > theta ( k )                           */
            #/* i. e.   twice   for  (x(x(x x)s)d)s type                   */
            #/*         sqrt(3) for  (x(x(x x)t)d)s type                   */
            #/* where x is a or c.                                         */
            #/* modified 84-08-06  SAS                                     */
            #/* re-modified 84-08-06  SAS                                  */
            #/* THETA_GT = MAX ( THETA, THETA_PREV )                       */
            #W =  (-1)**( ( F - 1 ) * ( THETA_PREV - 1 ) )
            # * MATRIX_EL ( GF, HF, FF )
            # * SQRT  ( FLOAT ( PSI * PHI_PREV * THETA, 53 ) )
            # * WIG9X ( H,   THETA_PREV,  THETA,
            #           F,   PHI_PREV,    PHI,
            #           G,   PSI_PREV,    PSI   ) ; /* rev. 84-09-04  */
    
            phase   = (-1)**( ( f - 1 ) * ( theta_prev - 1 ) )
            matelem = matrix_el( gf, hf, ff )
            mult    = Math.sqrt( ( psi * phi_prev * theta ).to_f )
            ninej   = Wigcoef.wig9x( h,  theta_prev, theta,
                                     f,  phi_prev,   phi,
                                     g,  psi_prev,   psi   )
    
            w = (-1)**( ( f - 1 ) * ( theta_prev - 1 ) ) *
                matrix_el( gf, hf, ff ) *
                Math.sqrt( ( psi * phi_prev * theta ).to_f ) *
                Wigcoef.wig9x( h,  theta_prev, theta,
                               f,  phi_prev,   phi,
                               g,  psi_prev,   psi   )
    
            if $opt_debug_print then
                printf( "  HF THETA THETA_PREV :  %6s %4d %4d\n", TBL_OPSTR[hf], theta, theta_prev )
                printf( "  GF PSI     PSI_PREV :  %6s %4d %4d\n", TBL_WFSTR[gf], psi,   psi_prev   )
                printf( "  FF PHI     PHI_PREV :  %6s %4d %4d\n", TBL_WFSTR[ff], phi,   phi_prev   )
                printf( "  phase = %14.6e, mult = %14.6e, 9j = %14.6e, mat = %14.6e, w = %14.6e\n",
                        phase, mult, ninej, matelem, w )
            end
        end
 
    rescue InvalidArgumentRange
        rtn_code = 1
        $stderr.print "HMLKVA: INVALID ARGUMENTS ?\n"
        $stderr.printf( "  HF THETA THETA_PREV :  %6s %4d %4d\n", TBL_OPSTR[hf], theta, theta_prev )
        $stderr.printf( "  GF PSI     PSI_PREV :  %6s %4d %4d\n", TBL_WFSTR[gf], psi,   psi_prev   )
        $stderr.printf( "  GF PHI     PHI_PREV :  %6s %4d %4d\n", TBL_WFSTR[ff], phi,   phi_prev   )
    #rescue InvalidLkValue
    #    rtn_code = 1
    #    $stderr.print "HMLKVA: ERROR CONDITION ?\n"
    #rescue InvalidLkValueAttention
    #    rtn_code = 2
    #    $stderr.print "HMLKVA: ATTENTION ?\n"
    #rescue InvalidBySystemError
    #    rtn_code = 10
    #    $stderr.print "HMLKVA: SYSTEM ERROR ?\n"
    end

    return w

end

############################################################
############################################################
############################################################

def matrix_el( gf, hf, ff )
    rtn_val = 0

#    p 'gf', gf
#    p 'hf', hf
#    p 'ff', ff
#    p 'FULL', FULL
#    p 'EMPTY', EMPTY
#    p 'CAS', CAS
#    p 'CAT', CAT

    begin
        #WHEN ( GF = EMPTY & FF = EMPTY )                                         
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( UNITY ) RTN_VAL = 1 ;              /* 84-09-03 */              
        #   OTHER ;                                                               
        #   END ;                                                                 
        if ( gf == EMPTY ) && ( ff == EMPTY ) then
            if hf == UNITY then
                rtn_val = 1
            elsif hf == ACS then                     # added by honda
                rtn_val = - Math.sqrt(2)
            else
                raise InvalidCombination
            end

        #WHEN ( ( GF = UP ! GF = DOWN ) & ( FF = UP ! FF = DOWN ) )               
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( UNITY ) RTN_VAL = SQRT ( TWO ) ;   /* 84-09-03 */              
        #   WHEN ( CAS   ) RTN_VAL = - 1 ;            /* 84-09-03 */              
        #   WHEN ( CAT   ) RTN_VAL = SQRT ( THREE ) ; /* 84-09-03 */              
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( gf == UP || gf == DOWN ) && ( ff == UP || ff == DOWN ) then
            if hf == UNITY then
                rtn_val = Math.sqrt( TWO )
            elsif hf == CAS then
                rtn_val = -1
            elsif hf == ACS then                     # added by honda
                rtn_val = -1
            elsif hf == CAT then
                rtn_val = Math.sqrt( THREE )
            elsif hf == ACT then                     # added by honda
#               rtn_val = - Math.sqrt( THREE )
                rtn_val = + Math.sqrt( THREE )       # tashi
            else
                raise InvalidCombination
            end

        #WHEN ( GF = FULL  & FF = FULL )                                         
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( UNITY ) RTN_VAL = 1 ;              /* 84-09-03 */              
        #   WHEN ( CAS   ) RTN_VAL = - SQRT ( TWO ) ; /* 84-09-03 */              
        #   WHEN ( CAT   ) RTN_VAL = 0 ;              /* 84-09-03 */              
        #   /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/             
        #   WHEN ( CCAA  ) RTN_VAL = 1 ;              /* not used */              
        #   /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/             
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( gf == FULL ) && ( ff == FULL ) then
            if hf == UNITY then
                rtn_val = 1
            elsif hf == CAS then
                rtn_val = - Math.sqrt( TWO )
                return rtn_val
            elsif hf == CAT then
                rtn_val = 0
            elsif hf == CCAA then          # not used
                rtn_val = 1
            else
                raise InvalidCombination
            end

        #WHEN ( GF = EMPTY & ( FF = UP ! FF = DOWN ) )                            
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( A     ) RTN_VAL = + SQRT ( TWO ) ; /* 84-09-03 */              
        #/* WHEN ( A     ) RTN_VAL = - SQRT ( TWO ) */                            
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( gf == EMPTY ) && ( ff == UP || ff == DOWN ) then
            if hf == A then
                rtn_val = + Math.sqrt( TWO )
            ###if hf == A then
            ### rtn_val = - Math.sqrt( TWO )
            else
                raise InvalidCombination
            end

        #WHEN ( FF = EMPTY & ( GF = UP ! GF = DOWN ) )                            
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( C     ) RTN_VAL = + SQRT ( TWO ) ; /* 84-09-03 */              
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( ff == EMPTY ) && ( gf == UP || gf == DOWN ) then
            if hf == C then
                rtn_val = + Math.sqrt( TWO )
            else
                raise InvalidCombination
            end

        #WHEN ( FF = FULL  & ( GF = UP ! GF = DOWN ) )                            
        #   SELECT ( HF ) ;                                                       
        #/* WHEN ( A     ) RTN_VAL = + SQRT ( TWO ) */                            
        #   WHEN ( A     ) RTN_VAL = - SQRT ( TWO ) ; /* 84-09-03 */              
        #   WHEN ( CAA   ) RTN_VAL = - 1 ;      /* C * AAS / 2  */                
        # /*            01-04   by SAS at HUCC    Sep.   6, 1984  */         
        #   OTHER ;                                   /* 84-09-03 */              
        #   END ;                                                                 
        elsif ( ff == FULL ) && ( gf == UP || gf == DOWN ) then
            if hf == A then
                rtn_val = + Math.sqrt( TWO )               ### noro/sas before 84-09-03 phase  (14-08-23)
            ### rtn_val = - Math.sqrt( TWO )               ### sas after 84-09-03 phase        (14-08-23)
            elsif hf == CAA then
                rtn_val = - 1                 ### C * AAS / 2
            else
                raise InvalidCombination
            end

        #WHEN ( GF = FULL  & ( FF = UP ! FF = DOWN ) )                            
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( C     ) RTN_VAL = + SQRT ( TWO ) ; /* 84-09-03 */              
        #   WHEN ( CCA   ) RTN_VAL = + 1 ;      /* CCS * A / 2  */                
        #   OTHER ;                                   /* 84-09-03 */              
        #   END ;                                                                 
        elsif ( gf == FULL ) && ( ff == UP || ff == DOWN ) then
            if hf == C then
#               rtn_val = + Math.sqrt( TWO )   #### tashi 2014-11-15
                rtn_val = - Math.sqrt( TWO )
            elsif hf == CCA then
                rtn_val = + 1                 ### CCS * A / 2
            else
                raise InvalidCombination
            end

        #WHEN ( GF = EMPTY & FF = FULL  )                                         
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( AAS   ) RTN_VAL = - ONE / SQRT ( TWO ) ;                       
        # /*               01-04   by SAS at HUCC    Sep.   6, 1984  */         
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( gf == EMPTY ) && ( ff == FULL ) then
            if hf == AAS then
            ### rtn_val = - ONE / Math.sqrt( TWO )         ### kamuy value (14-08-23)
                rtn_val = + Math.sqrt( TWO )               ### modified    (14-08-23)
            else
                raise InvalidCombination
            end

        #WHEN ( GF = FULL  & FF = EMPTY )                                         
        #   SELECT ( HF ) ;                                                       
        #   WHEN ( CCS   ) RTN_VAL = + ONE / SQRT ( TWO ) ;                       
        # /*               01-04   by SAS at HUCC    Sep.   6, 1984  */         
        #   OTHER ;                                                               
        #   END ;                                                                 
        elsif ( gf == FULL ) && ( ff == EMPTY ) then
            if hf == CCS then
            ### rtn_val = + ONE / Math.sqrt( TWO )
                rtn_val = + Math.sqrt( TWO )               ### modified (14-08-23)
            else
                raise InvalidCombination
            end

        end

        #OTHER  PUT SKIP LIST                                                     
        #   ('INVALID COMBINATION OF GF AND HF IN HMLSET ' ) ;                    
        #END ;                                                                    
    rescue InvalidCombination
        $stderr.printf( "INVALID COMBINATION OF GF AND HF IN HMLSET : gf,hf,ff = %s %s %s\n",
                        TBL_WFSTR[gf], TBL_OPSTR[hf], TBL_WFSTR[ff] )
    end

    return rtn_val

    #IF @@SUB@ >= 4 THEN PUT SKIP DATA ( RTN_VAL ) ;                             
    #RETURN ( RTN_VAL ) ;                                                        
    #END MATRIX_EL ;                                                             
    #                                                                            
    #%  INCLUDE  ERRSUB ;                                                        
    #END @SUB@ ;                                                                 
end

Debug = false
DebugDebug = false
Thresh_zero = 1.0E-7
Seg_walk = [[],
            [[EMPTY, EMPTY], [UP, UP], [DOWN, DOWN], [FULL, FULL], [UP, DOWN], [DOWN, UP]], # UNITY
            [[EMPTY, UP], [EMPTY, DOWN], [UP, FULL], [DOWN, FULL]],                         # A
            [[UP, EMPTY], [DOWN, EMPTY], [FULL, UP], [FULL, DOWN]],                         # C
            [[EMPTY, FULL]],                                                                # AAS
            [[UP, UP], [DOWN, DOWN], [FULL, FULL], [UP, DOWN], [DOWN, UP]],                 # CAS
            [[UP, UP], [DOWN, DOWN], [UP, DOWN], [DOWN, UP]],                               # CAT
            [[FULL, EMPTY]],                                                                # CCS
            [[UP, FULL], [DOWN, FULL]],                                                     # CAA
            [[FULL, UP], [FULL, DOWN]],                                                     # CCA
            [[FULL, FULL]]]                                                                 # CCAA

Num_ca = [ 0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 4 ] # number of creation or anihilation operators in UNITY, ..., CCAA
Walk_nam = [ "", "e", "u", "d", "f" ]
Op_nam = ["", "UNITY", "A", "C", "AAS", "CAS", "CAT", "CCS", "CAA", "CCA", "CCAA"]
$maxterm = 4

class TensorProduct

  def initialize( _op, _theta )
    @op = _op.clone 
    @theta = _theta.clone
    @theta_prev = [1] 
    (0...@theta.size-1).each do |i|
      @theta_prev.push(@theta[i])
    end
  end

  def lenop
    return @op.length
  end

  def num_ca
    tmp = 0
    @op.each do |x|
      tmp += Num_ca[ x ]
    end
    tmp
  end

  def op(iop)  
    return @op[iop]
  end

  def hf(iop)
    o = op(iop)
    if ( o == UNITY ) || ( o == CAS ) || ( o == ACS    ) ||
       ( o == CCS   ) || ( o == AAS ) || ( o == CASCAS ) then
      return 1
    elsif ( o == A ) || ( o == C ) || ( o == CCA ) || ( o == CAA ) then
      return 2
    elsif ( o == CAT ) || ( o == ACT ) then
      return 3
    else
      $stderr.printf( "Error: $opstbl::hf: illegal opcode: %s\n", o )
      exit
    end
  end

  def theta(iop)
    return @theta[iop]
  end

  def theta_prev(iop)
    return @theta_prev[iop]
  end

  def show
    print "Operator = "
    tmp = []
    @op.each do |a|
      tmp.push( Op_nam[a] )
    end
    p tmp
    print "Theta = "
    p @theta
    print "Theta_prev = "
    p @theta_prev
  end
end

class LoopElement
  attr_accessor :bottom, :top, :weight_f, :weight_g, :ijkl, :value, :seq_num 
#
  def initialize(_bottom, _top, _weight_g, _weight_f, _ijkl, _value, _origin)
    @bottom = _bottom
    @top = _top
    @weight_f = _weight_f
    @weight_g = _weight_g
    @ijkl = _ijkl.clone
    @seq_num = nil
    @value = _value
    @origin = _origin
  end

  def top
    @top
  end

  def show
    printf(" %5d   %5d    %5d    %5d    %12.8f   %5d    %5d    ", @bottom, @top, @weight_g, @weight_f, @value, @seq_num, @origin)
    p @ijkl
  end
end

def calc_TensorOp_average( drt, h_op, origin )

    all_expr = []
    iop = 0                                 # current operator number
# op                                           TensorProduct
    kmax = drt.r_upperlevel( drt.nnode - 1 )
    ng = Array.new( size = kmax, val = nil )      # node of psi (k)
    nf = Array.new( size = kmax, val = nil )      # node of phi (k)
    wg = Array.new( size = kmax, val = nil )      # weight of psi (k)
    wf = Array.new( size = kmax, val = nil )      # weight of phi (k)
    mk = Array.new( size = kmax, val = nil )      # value Lk (k)
    ijkl = []

    if Debug == true
        print "\nTensor Product", "\n"
        h_op.show
        print "\n"
    end

    inode = 0
    while kmax - drt.r_upperlevel( inode ) > h_op.lenop - 1
        expr = []; record = {}
        k = drt.r_upperlevel( inode ) 
        nf[ k ] = inode; ng[ k ] = inode; theta = 1; mk[ k ] = 1
        wg[ k ] = 0; wf[ k ] = 0
        hf = h_op.op( iop )
        if Debug == true
            print "===== from node ", ng[ k ], " with op ", Op_nam[h_op.op( iop )], "\n"
        end
        Seg_walk[ hf ].each do | walk |
            if Debug == true
                print "=====     walk : ", Walk_nam[walk[0]], ", ", Walk_nam[walk[1]], "\n"
            end
            theta = h_op.theta( iop )
            theta_prev = h_op.theta_prev( iop )
            tree_search( drt, inode, walk, k, ng, nf, wg, wf, 
                         iop, h_op, hf, theta, theta_prev, mk, ijkl, origin, expr, record )
        end
#       all_expr.push( expr.sort{ |x, y| x.top - y.top } )
        all_expr.push( expr.sort{ |x, y| x.top - y.top }.uniq )
        inode += 1
    end

    return all_expr.flatten
end

def tree_search( drt, bottom, walk, k, ng, nf, wg, wf, 
                 iop, h_op, hf, theta, theta_prev, mk, ijkl, origin, expr, record ) 

    if DebugDebug == true    
        print "===> tree_search \n"
        print "     k = ", k, " iop = ", iop, "\n" 
        print "     ng, nf : ", ng, ",  ", nf
    end
    next_ng = drt.r_uppercon( ng[ k ], walk[ 0 ] - 1 )
    next_nf = drt.r_uppercon( nf[ k ], walk[ 1 ] - 1 )
    if DebugDebug == true    
        print ";    next_ng, nf : ", next_ng, ", ", next_nf
        print ";    hf : ", hf, "\n"
    end
    num_ca = 0
    if next_ng != nil && next_nf != nil
        k += 1
        ng[ k ] = next_ng
        nf[ k ] = next_nf
        wg[ k ] = wg[ k -1 ] + drt.r_upperarcw( ng[ k - 1 ], walk[ 0 ] - 1 )
        wf[ k ] = wf[ k -1 ] + drt.r_upperarcw( nf[ k - 1 ], walk[ 1 ] - 1 )
        num_ca = Num_ca[ hf ]
        num_ca.times do |i|
            ijkl.push( k )
        end 
        psi = drt.r_upperspin( ng[ k ] ); phi = drt.r_upperspin( nf[ k ] )
        mk[ k ] = mk[ k - 1 ] * hmlkva( walk[ 0 ], walk[ 1 ], hf, psi, phi, theta, theta_prev )

        if DebugDebug == true    
            print "     g, f : ", Walk_nam[walk[0]], ", ", Walk_nam[walk[1]]
            print "    hf : ", hf
            print "    psi, phi : ", psi, ", ", phi
            print "    theta, theta_prev : ", theta, ", ", theta_prev, "\n"
            print "     mk_value = ", mk[k], "\n \n"
        end

        if mk[ k ].abs < Thresh_zero 
            k -= 1
            num_ca.times do |i|
                ijkl.pop
            end
            if DebugDebug == true    
                p "<=== tree_search" 
            end
            return
        end

        if iop == h_op.lenop - 1 && ng[ k ]  == nf[ k ]
            if wg[ k ] > wf[ k ]
                key = [ bottom, ng[ k ], wg[ k ], wf[ k ], ijkl ].flatten
            else
                key = [ bottom, ng[ k ], wf[ k ], wg[ k ], ijkl ].flatten
            end
            if !record.key?( key )
                tmp = LoopElement.new( bottom, ng[ k ], wg[ k ], wf[ k ], ijkl, mk[ k ], origin )
                expr.push( tmp )
                record.store( key, true )
                if Debug == true
                    print "Loop closed \n"
                    print "g = "
                    (1..k).each do |it|
                         printf( " %3d ", ng[it] )
                    end
                    print "\nf = "
                    (1..k).each do |it|
                         printf( " %3d ", nf[it] )
                    end
                    print "\n"
                    tmp.show
                end
            end
        else
            if ng.size - k > h_op.lenop - iop - 1
                if DebugDebug == true    
                    print "    UNITY \n"
                end
                hf = UNITY
                theta = theta
                theta_prev = theta
                Seg_walk[UNITY].each do |wk|
                    if DebugDebug == true    
                        print "    walk : ", Walk_nam[wk[0]], ", ", Walk_nam[wk[1]], "\n"
                    end
                    tree_search( drt, bottom, wk, k, ng, nf, wg, wf, 
                                 iop, h_op, hf, theta, theta_prev, mk, ijkl, origin, expr, record ) 
                end
            end
#
            if h_op.lenop - 1 > iop
                iop += 1
                hf = h_op.op( iop )
                theta_prev = theta
                theta = h_op.theta( iop )
                if DebugDebug == true    
                    print "    Next Operator ", hf, "\n"
                end
                Seg_walk[ hf ].each do |wk|
                if DebugDebug == true    
                        print "    walk : ", Walk_nam[wk[0]], ", ", Walk_nam[wk[1]], "\n"
                end
                     tree_search( drt, bottom, wk, k, ng, nf, wg, wf, 
                                 iop, h_op, hf, theta, theta_prev, mk, ijkl, origin, expr, record ) 
                end
                iop -= 1
            end
            k -= 1
        end
    end
    num_ca.times do |i|
        ijkl.pop
    end 
    if DebugDebug == true    
        p "<=== tree_search"
    end
    return
end

class BrooksCases
    attr_reader :drt, :norb, :expr, :expr_one
    def initialize( _drt )
        @drt = _drt.clone
        csf_list = @drt.mk_csf_list
        @norb = _drt.norb
        @num_coefficient = 0
        expr_all = self.all
        @expr = Expression.new( [@norb, csf_list, expr_all] ) 
        tmp_one = []; one_size = @norb * ( @norb + 1 ) / 2
        expr_all.each do | x |
            if x[0][1] < one_size
                tmp_one.push( x )
            end
        end
        @expr_one = Expression.new( [@norb, csf_list, tmp_one] )
    end
 
    def triangle( a, b )
        if a > b
            x = a; y = b 
        else
            x = b; y = a
        end
        x * ( x + 1 ) / 2 + y
    end

    def tri( p, q )
        if p > q
            x = p; y = q 
        else
            x = q; y = p
        end
        x * ( x - 1 ) / 2 + y
    end

    def expand_expr( expr_in )
        expr_out = []
        expr_in.each do | t |
            lower_weight = @drt.get_lower_path_weights( t[1].bottom )
            upper_weight = @drt.get_upper_path_weights( t[1].top )
            lower_weight.each do | x |
                upper_weight.each do | y  |
                    key = [ self.triangle( x + y + t[1].weight_g, x + y + t[1].weight_f ), t[1].seq_num ]
                    expr_out.push( [ key, t[1] ] )
                end
            end
        end
        expr_out.sort{ |x, y| x[0][0] - y[0][0] }
    end

    def add_expr( expr_in, coef, a )
        expr_out = {}
        n_oel = @norb * ( @norb + 1 ) / 2
        expr_in.each_with_index do |t, i|
            t.each do |x|
                @num_coefficient += 1
                w = @drt.get_lower_path_weights( x.bottom )[0] + @drt.get_upper_path_weights( x.top )[0]

                if a.size > 2
                   tmp = [ x.ijkl[ a[0] ], x.ijkl[ a[1] ], x.ijkl[ a[2] ], x.ijkl[ a[3] ] ]
                   x.seq_num = tri( tri(tmp[ 0 ],tmp[ 1 ] ),  tri(tmp[ 2 ],tmp[ 3 ] ) ) - 1 + n_oel
                else
                   tmp = [ x.ijkl[ a[0] ], x.ijkl[ a[1] ] ]
                   x.seq_num = tri(tmp[ 0 ],tmp[ 1 ] ) - 1
                end
                key = [ self.triangle( w + x.weight_g, w + x.weight_f ), x.seq_num ]

                if expr_out.key?( key )
                    y = expr_out[ key ]
                    y.value += coef[i] * x.value
                else
                    y = x.clone
                    y.value *= coef[i]
                    y.ijkl = tmp.clone
                end
                if y.value.abs > 0.00001
                    expr_out.store( key, y )
                else
                    expr_out.delete( key )
                end
            end
        end
        expr_out
    end

    def get_num_coefficient
        @num_coefficient
    end

    def case01
        op = [ TensorProduct.new( [ CAS ], [ 1 ] ) ]
        expr = []
        expr.push( calc_TensorOp_average( @drt, op[0], 0 ) )
        expr_ll = add_expr( expr, [ -Sqrt2 ], [ 0, 1 ] )
        return expr_ll
    end

    def case02
        op = [ TensorProduct.new( [ C, A ], [ 2, 1 ] ) ]
        expr = []
        expr.push( calc_TensorOp_average( @drt, op[0], 0 ) )
        expr_ij = add_expr( expr, [ -Sqrt2 ], [ 0, 1 ] )
        return expr_ij
    end

    def case03
        op = [ TensorProduct.new( [ A, C ], [ 2, 1 ] ) ]
        expr = []
        expr = expr.push( calc_TensorOp_average( @drt, op[0], 0 ) )
        expr_ij = add_expr( expr, [ -Sqrt2 ], [ 0, 1 ] )
        return expr_ij
    end

    def case1
        op = [ TensorProduct.new( [ A, C, C, A ], [ 2, 1, 2, 1 ] ), 
               TensorProduct.new( [ A, C, C, A ], [ 2, 3, 2, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 1 ) )
        end

        expr_ikjl = add_expr( expr, [ -1, -Sqrt3 ], [ 0, 2, 1, 3 ] )
        expr_ijkl = add_expr( [expr[0]], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_ikjl, expr_ijkl
    end

    def case2
        op = [ TensorProduct.new( [ C, A, C, A ], [ 2, 1, 2, 1 ] ), 
               TensorProduct.new( [ C, A, C, A ], [ 2, 3, 2, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 2 ) )
        end

        expr_iljk = add_expr( expr, [ -1, Sqrt3 ], [ 0, 3, 1, 2 ] )
        expr_ijkl = add_expr( [expr[0]], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_iljk, expr_ijkl
    end

    def case3
        op = [ TensorProduct.new( [ C, C, A, A ], [ 2, 1, 2, 1 ] ), 
               TensorProduct.new( [ C, C, A, A ], [ 2, 3, 2, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 3 ) )
        end

        expr_ikjl = add_expr( expr, [ -1, -Sqrt3 ], [ 0, 2, 1, 3 ] )
        expr_iljk = add_expr( expr, [ -1, Sqrt3 ], [ 0, 3, 1, 2 ] )

        return expr_ikjl, expr_iljk
    end

    def case4
        op = TensorProduct.new( [ A, CCS, A ], [ 2, 2, 1 ] )

        expr_ijjk = calc_TensorOp_average( @drt, op, 4 )
#       expr_ijjk = add_expr( [expr_ijjk], [ -1 ], [ 0, 1, 2, 3 ] )
        expr_ijjk = add_expr( [expr_ijjk], [ 1 ], [ 0, 1, 2, 3 ] )

        return expr_ijjk
    end

    def case5
        op = [ TensorProduct.new( [ C, CAS, A ], [ 2, 2, 1 ] ), 
               TensorProduct.new( [ C, CAT, A ], [ 2, 2, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 5 ) )
        end

        expr_ijjl = add_expr( expr, [ -1, Sqrt3 ], [ 0, 1, 2, 3 ] )
        expr_jjil = add_expr( [expr[0]], [ 2 ], [ 0, 3, 1, 2 ] )

        return expr_ijjl, expr_jjil
    end

    def case6
        op = TensorProduct.new([ CCS, A, A ], [ 1, 2, 1 ] )
        expr_ilik = calc_TensorOp_average( @drt, op, 6 )

#       expr_ilik = add_expr( [expr_ilik], [ -1 ], [ 0, 3, 1, 2 ] )
        expr_ilik = add_expr( [expr_ilik], [ 1 ], [ 0, 3, 1, 2 ] )

        return expr_ilik
    end

    def case7
        op = [ TensorProduct.new( [ CAS, C, A ], [ 1, 2, 1 ] ), 
               TensorProduct.new( [ CAT, C, A ], [ 3, 2, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 7 ) )
        end

        expr_ikil = add_expr( expr, [ -1, -Sqrt3 ], [ 0, 2, 1, 3 ] )
        expr_iikl = add_expr( [expr[0]], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_ikil, expr_iikl
    end

    def case8
        op = TensorProduct.new([ C, C, AAS ], [ 2, 1, 1 ] )
        expr_iljl = calc_TensorOp_average( @drt, op, 8 )

        expr_iljl = add_expr( [expr_iljl], [ -1 ], [ 0, 3, 1, 2 ] )
#       expr_iljl = add_expr( [expr_iljl], [  1 ], [ 0, 3, 1, 2 ] )

        return expr_iljl
    end
 
    def case9
        op = [ TensorProduct.new( [ A, C, CAS ], [ 2, 1, 1 ] ), 
               TensorProduct.new( [ A, C, CAT ], [ 2, 3, 1 ] ) ]
        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 9 ) )
        end

        expr_iljl = add_expr( expr, [ -1, Sqrt3 ], [ 0, 2, 1, 3 ] )

        return expr_iljl
    end

    def case10
        op = [ TensorProduct.new( [ C, A, CAS ], [ 2, 1, 1 ] ), 
               TensorProduct.new( [ C, A, CAT ], [ 2, 3, 1 ] ) ]

        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 10 ) )
        end

        expr_iljl = add_expr( expr, [ -1, -Sqrt3 ], [ 0, 2, 1, 3 ] )
        expr_ijll = add_expr( [ expr[ 0 ] ], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_iljl, expr_ijll
    end

    def case11
        op = [ TensorProduct.new( [ C, CAA ], [ 2, 1 ] ), 
               TensorProduct.new( [ CCA, A ], [ 2, 1 ] ) ]

        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 11 ) )
        end

#       expr_illl = add_expr( [ expr[ 0 ] ], [ -2 ], [ 0, 1, 2, 3 ] )
#       expr_iiil = add_expr( [ expr[ 1 ] ], [ -2 ], [ 0, 1, 2, 3 ] )
        expr_illl = add_expr( [ expr[ 0 ] ], [ 2 ], [ 0, 1, 2, 3 ] )
        expr_iiil = add_expr( [ expr[ 1 ] ], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_illl, expr_iiil
    end

    def case12
        op = TensorProduct.new( [ CCS, AAS ], [ 1, 1 ] )
        expr = calc_TensorOp_average( @drt, op, 12 )

#       expr_ilil = add_expr( [ expr ], [ -0.5 ], [ 0, 2, 1, 3 ] )
        expr_ilil = add_expr( [ expr ], [ 0.5 ], [ 0, 2, 1, 3 ] )

        return expr_ilil
    end

    def case13
        op = [ TensorProduct.new( [ CAS, CAS ], [ 1, 1 ] ), 
               TensorProduct.new( [ CAT, CAT ], [ 3, 1 ] ) ]

        expr = []
        op.each do |x|
            expr.push( calc_TensorOp_average( @drt, x, 13 ) )
        end

        expr_ilil = add_expr( expr, [ -1, Sqrt3 ], [ 0, 2, 1, 3 ] )
        expr_iill = add_expr( [ expr[ 0 ] ], [ 2 ], [ 0, 1, 2, 3 ] )

        return expr_ilil, expr_iill
    end

    def case14
        op = TensorProduct.new( [ CCAA ], [ 1 ] )
        expr = calc_TensorOp_average( @drt, op, 14 )

        expr_llll = add_expr( [ expr ], [ 1.0 ], [ 0, 1, 2, 3 ] )

        return expr_llll
    end

    def all
        run_flag = Array.new( 15, true )
  #     run_flag = Array.new( 15, false )
        run_flag[ 9 ] = false

        expr = []
        # case01
        if run_flag[ 0 ]
            expr.push( self.case01.to_a.sort{ |a, b| a[0] <=> b[0] } )
        # case02
            expr.push( self.case02.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case1
        if run_flag[ 1 ]
            tmp1, tmp2 = self.case1
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case2
        if run_flag[ 2 ]
            tmp1, tmp2 = self.case2
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case3
        if run_flag[ 3 ]
            tmp1, tmp2 = self.case3
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case4
        if run_flag[ 4 ]
            expr.push( self.case4.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case5
        if run_flag[ 5 ]
            tmp1, tmp2 = self.case5
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case6
        if run_flag[ 6 ]
            expr.push( self.case6.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case7
        if run_flag[ 7 ]
            tmp1, tmp2 = self.case7
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case8
        if run_flag[ 8 ]
            expr.push( self.case8.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case9
        if run_flag[ 9 ]
            expr.push( self.case9.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case10
        if run_flag[ 10 ]
            tmp1, tmp2 = self.case10
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case11
        if run_flag[ 11 ]
            tmp1, tmp2 = self.case11
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case12
        if run_flag[ 12 ]
            expr.push( self.case12.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case13
        if run_flag[ 13 ]
            tmp1, tmp2 = self.case13
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
            expr.push( tmp2.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        # case14
        if run_flag[ 14 ]
            tmp1 = self.case14
            expr.push( tmp1.to_a.sort{ |a, b| a[0] <=> b[0] } )
        end

        a = []
        expr.each do | x |
            x.each do | y |
                a.push( y )
            end
        end        

        expr_sorted = a.sort{ |x, y| x[0][0] - y[0][0] }
        printf("\nNumber of records of sorted exression  %8d \n", expr_sorted.size )

#        printf( "\nSummary of sorted exression \n" )
#        expr_sorted.each do |x|
#            printf("[ %5d, %5d ]", x[0][0], x[0][1])
#            x[1].show
#        end 
#        printf( "\n" )

        expr_expand = self.expand_expr( expr_sorted )
#       printf("\nNumber of records of expanded exression  %8d \n", expr_expand.size )
#       printf( "\nSummary of expanded exression \n" )
#       expr_expand.each do |x|
#           printf("[ %5d, %5d ]", x[0][0], x[0][1])
#           x[1].show
#       end 
#       printf( "\n" )

        return expr_expand
    end
end

def show( a )
    a.each do |x|
       x.show
    end
end

class Expression
    attr_accessor :nexpr, :norb, :ij, :pqrs, :coef, :csf_list
    def initialize( obj = nil )
        if obj.nil? then
            @nexpr = 0
            @norb = 0
            @csf_list = []
            @ij = []
            @pqrs = []
            @coef = []
        end

        if ( obj.instance_of?( Array ) && obj.size == 3 ) then   #
            @norb = obj[ 0 ]
            @csf_list = obj[ 1 ].clone
            @ij = []
            @pqrs = []
            @coef = []
            obj[2].each do | x |
                @ij.push( x[0][0] + 1 )
                @pqrs.push( x[0][1] )
                @coef.push( x[1].value )
            end
            @nexpr = @ij.length
        elsif ( obj.instance_of?( Expression ) ) then
            @norb = obj.norb
            @csf_list = obj.csf_list.clone
            @ij = obj.ij.clone
            @pqrs = obj.pqrs.clone
            @coef = obj.coef.clone 
            @nexpr = @ij.length
        elsif ( obj.instance_of?( Array ) && obj.size == 5 ) then
            @norb = obj[ 0 ]
            @csf_list = obj[ 1 ]
            @ij = obj[ 2 ]
            @pqrs = obj[ 3 ]
            @coef = obj[ 4 ]
            @nexpr = @ij.length
        else
        end
    end

    def show
        printf( " CSF List \n" )
        @csf_list.each do | x |
            p x
        end
        printf( " Energy Expression \n" )
        printf( "      IJ       pqrs     Coef \n" )
        ( 0...@ij.size ).each do | i |
            printf( " %7d    %10d    %24.16e \n", @ij[i], @pqrs[i], @coef[i] )
        end
    end
end
