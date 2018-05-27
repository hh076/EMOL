# -*- coding: utf-8 -*-
#!/usr/bin/ruby
#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative 'emol_drt_soci'
require_relative 'molinfo'
require_relative 'emol_trans'
require_relative 'emol_expression'
require_relative 'emol_hmatrix'
require_relative 'liu_v12'
require_relative 'emol_futil'
require_relative 'emol_util'
require_relative 'emol_consts'

def pair( p, q )
    if ( p > q )
        p * ( p - 1 ) / 2 + q
    else
        q * ( q - 1 ) / 2 + p
    end
end

class H_matrix_casscf < H_matrix
    def initialize( molinfo, expr, h_mo, eris_mo, core_energy )
        @n = expr.csf_list.size
        nstate = molinfo.nstate
        thresh = EmolConsts::THR_CONV_WF
#                   H 行列の生成
        mkhmat( expr, h_mo, eris_mo )
#                   H 行列の Liu によるアルゴリズムによる対角化
        eigen_val, @eigen_vector = Liu.new( self, nstate, thresh ).solve
#                   凍結殼エネルギーの付加
        @eigen_value = GSL::Vector.alloc( nstate )
        (0...nstate).each do | k |
            @eigen_value[ k ] = eigen_val[ k ] + core_energy
            printf( "     %5d   %20.15f \n", k, @eigen_value[ k ] )
        end
    end
end

class MO_classification_mcscf
   attr_reader :norb, :ncore, :nactive, :nexternal, :core_range, :act_range, :ext_range, :ncore_act, :ncore_ext, :nact_ext, :nrotations
   def initialize( _ncore, _nactive, _nexternal )
#           _ncore     : 閉殼軌道数
#           _nact      : 自由な個数の電子占有可能な軌道数
#           _nexternal : 空軌道数
      @ncore = _ncore
      @nactive  = _nactive
      @nexternal  = _nexternal
      @norb  = @ncore + @nactive + @nexternal
      @core_range = 0...@ncore
      @act_range = @ncore...@ncore + @nactive
      @ext_range = @ncore + @nactive...@norb
      @ncore_act = @ncore * @nactive
      @ncore_ext = @ncore * @nexternal
      @nact_ext  = @nactive * @nexternal
      @nrotations = @ncore_act + @ncore_ext + @nact_ext
      it = -1
      printf( "******Orbital Rotation \n" )
      for i in @core_range
          for t in @act_range
              it += 1
              printf( "%5d   %5d   %5d\n", it, i, t )
          end
      end
      for i in @core_range
          for a in @ext_range
              it += 1
              printf( "%5d   %5d   %5d\n", it, i, a )
          end
      end
      for t in @act_range
          for a in @ext_range
              it += 1
              printf( "%5d   %5d   %5d\n", it, t, a )
          end
      end
   end
end

class Transformation_mcscf < Transformation
    def initialize( obj = nil )
        if ( obj.instance_of?( Array ) ) then
            mk_moint( obj[ 0 ], obj[ 1 ], obj[ 2 ] )
        end
    end

    def mk_moint( obj0, obj1, cao )     
        load obj0
        molinfo = Data_molinfo::molinfo_data

        load obj1
        s = Data_molint::s_data
        h = Data_molint::h_data
        nz_g = Data_molint::nz_g_data

#       load obj[ 2 ]
#       cao = Data_rhf::cao_data

        nbf = molinfo.get_nbf
        npq = nbf * (nbf + 1) / 2
        n_frozen = molinfo.get_n_frozen
        nob = molinfo.get_nob

        dao = Density.new( cao, n_frozen, 2 ) 
        g = calc_gmat_nz_eris( nz_g, dao.get_body, nbf )
        f = h + g
        @core_energy = (dao.get_body * (h + f)).trace / 2.0 + molinfo.get_e_nucl
#
        print "\n \n ******T r a n s f o r m a t i o n \n \n"
        printf( "nbf, n_frozen, nob : %5d  %5d  %5d  \n", nbf, n_frozen, nob )
        printf( "frozen core energy : %14.8f \n", @core_energy) 

        cao_trans = GSL::Matrix.alloc( nbf, nob )
        for i in 0...nob
            for j in 0...nbf
                cao_trans[ j, i ] = cao[ j, n_frozen + i ]
            end
        end

#           EmolUtil::print_2dim_ary(cao.to_a, "Original MO matrix", 8, "%12.6f")
        EmolUtil::print_2dim_ary(cao_trans.to_a, "Transformation matrix", 8, "%12.6f")

        s_work, h_work, eris_work, cao_1dim = two_one_dimCopy( nbf, nob, s, f, nz_g, cao_trans )
#           一電子積分は凍結殼との相互作用を含む Fock を変換する ( h ではなく f )
        cpq = GSL::Vector.alloc( nbf * nob )
        @s_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
        @h_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
#            Inttrans::trnpp( 1, 1, nbf, nob, s_work, s_mo, 1, cao_1dim, cpq )
#            Inttrans::trnpp( 1, 1, nbf, nob, h_work, h_mo, 1, cao_1dim, cpq )
        Inttrans::trnpp( 1, 1, nbf, nob, s_work, @s_mo, 1, cao_1dim, cpq )
        Inttrans::trnpp( 1, 1, nbf, nob, h_work, @h_mo, 1, cao_1dim, cpq )

#            EmolUtil::print_ary_tri(s_mo.to_a, "Overlap matrix", 10, "%12.6f")
#            EmolUtil::print_ary_tri(h_mo.to_a, "H matrix", 10, "%12.6f")
        EmolUtil::print_ary_tri(@s_mo.to_a, "Overlap matrix", 10, "%12.6f")
        EmolUtil::print_ary_tri(@h_mo.to_a, "H matrix", 10, "%12.6f")

        mi = 1
        mpq = nob * ( nob + 1 ) / 2
        for n in 1..npq
            Inttrans::trnpp(mi, mi, nbf, nob, eris_work, eris_work, 1, cao_1dim, cpq)
            mi += npq
        end
        Inttrans::trnsps(npq, mpq, eris_work)

        buf = GSL::Vector.alloc( npq )
        @eris_mo = GSL::Vector.alloc( mpq*(mpq + 1)/2 + 1 )
        m = 1; kin = 0
        for n in 1..mpq
            Inttrans::trnrr(m, n, nbf, nob, eris_work, buf, cao_1dim, cpq)
            m += npq
            for k in 0...n
                @eris_mo[ kin + k ] = buf[ k ]
            end
            kin += n
        end
#        EmolUtil::print_ary_tri(eris_mo.to_a, "Eris matrix", 10, "%12.6f")

        printf( "Total energy calculated from mo integrals : %14.8f \n", 
        @core_energy + close_shell_energy(h_mo, eris_mo, molinfo.get_nel - n_frozen * 2) )
    end
end

class SX
    def initialize( sx )
        @sx = sx.clone
        @n = sx.size[ 0 ]
    end

    def get_n
        @n
    end

    def init_vec( dummy )
        v = GSL::Matrix.alloc( @n, 1 )
        v[ 0, 0 ] = 1.0
        for i in 1...@n
            v[ i, 0 ] = 0.00001
        end
        v
    end

    def get_diag( p )
        @sx[ p, p ]
    end

    def get_hc( b )
        ab = GSL::Matrix.alloc( @n, 1 )
        for i in 0...@n
            sum = 0.0
            for k in 0...@n
                sum += @sx[ i, k ] * b[ k, 0 ]
            end
            ab[ i, 0 ] = sum
        end
        ab
    end
end

def extract_active_integral( mo, core_energy, h_mo, eris_mo )
#   CASCI 用に active space の MO 積分を拔き出す -> h_mo_active, eris_mo_active
#   
#   h に core からの寄与を加える -> h_mo_active
#   core_energy に core 部分の寄与を加える。
    for i in mo.core_range
        ii = pair( i + 1, i + 1 )
        core_energy += 2.0 * h_mo[ ii - 1 ] + eris_mo[ pair( ii, ii ) - 1 ]
        for j in 0...i
            ij = pair( i + 1, j + 1 )
            jj = pair( j + 1, j + 1 )
            core_energy += 4.0 * eris_mo[ pair( ii, jj ) - 1 ] - 2.0 * eris_mo[ pair( ij, ij ) - 1 ]
        end
    end

#   active 部分の積分を拔き出し、core からの寄与を h に加える
    h_mo_active = []
    eris_mo_active = []
    for i in mo.act_range
        for j in mo.act_range.first..i
            ij = pair( i + 1, j + 1 )
            tmp = h_mo[ pair( i + 1, j + 1 ) - 1 ]
            for it in mo.core_range
                ik = pair( i + 1, it + 1 )
                jk = pair( j + 1, it + 1 )
                kk = pair( it + 1, it + 1 )
                tmp += 2.0 * eris_mo[ pair( ij, kk ) - 1 ] - eris_mo[ pair( ik, jk ) - 1 ]
            end
            h_mo_active.push( tmp )
            for k in mo.act_range.first..i
                 for l in mo.act_range.first..k
                     kl = pair( k + 1, l + 1 )
                     if kl > ij
                         break
                     end
                     eris_mo_active.push( eris_mo[ pair( ij, kl ) - 1 ] )
                 end
            end
        end
    end
    return core_energy, h_mo_active, eris_mo_active
end

def mkdensity( mo, cas_expr, cf )
    one_size = mo.nactive * ( mo.nactive + 1 ) / 2
    two_size = one_size * ( one_size + 1 ) / 2
    dns1 = GSL::Vector.alloc( one_size )
    dns2 = GSL::Vector.alloc( two_size )
    dns1_unfold = GSL::Vector.alloc( one_size )
    dns2_unfold = GSL::Vector.alloc( two_size )
    prev_ij = 0; i = 1; j = 0
    for k in 0...cas_expr.ij.size
        if cas_expr.ij[ k ] > prev_ij
            if i == j then
                i += 1
                j = cas_expr.ij[ k ] - i * ( i - 1 ) / 2
            else
                j += cas_expr.ij[ k ] - prev_ij
            end
            prev_ij = cas_expr.ij[ k ]

            if i == j
                cij = cf[ i - 1, 0 ] * cf[ j - 1, 0 ]
            else
                cij = 2 * cf[ i - 1, 0 ] * cf[ j - 1, 0 ]
            end
        end
        if cas_expr.pqrs[ k ] >= one_size
            dns2[ cas_expr.pqrs[ k ] - one_size ] += cij * cas_expr.coef[ k ]
        else
            dns1[ cas_expr.pqrs[ k ] ] += cij * cas_expr.coef[ k ]
        end
    end

    for i in 0...mo.nactive
        for j in 0...i + 1
            ij = pair( i + 1, j + 1 ) 
            if i == j
                ij_fold = 0
            else
                ij_fold = 1
            end
            dns1_unfold[ ij - 1 ] = dns1[ ij - 1 ] / 2** ij_fold

            for k in 0...i + 1
                for l in 0...k + 1
                    kl = pair( k + 1, l + 1 )
                    if ij < kl
                        break
                    end
                    ijkl = pair( ij, kl )
                    if k == l
                        kl_fold = 0
                    else
                        kl_fold = 1
                    end
                    if ij == kl
                        ijkl_fold = 0
                    else
                        ijkl_fold = 1
                    end
                    dns2_unfold[ ijkl - 1 ] = dns2[ ijkl - 1 ] / 2.0**( ij_fold + kl_fold + ijkl_fold )
                end
            end
        end
    end           
    EmolUtil::print_ary_tri(dns1, "density1", 10, "%14.8f")
    EmolUtil::print_ary_tri(dns1_unfold, "density_unfold", 10, "%14.8f")
    EmolUtil::print_ary_tri(dns2, "density2", 10, "%14.8f")
    EmolUtil::print_ary_tri(dns2_unfold, "density2_unfold", 10, "%14.8f")
    return dns1, dns2, dns1_unfold, dns2_unfold
end

def check_density( mo, dns1, dns2, h_mo_active, eris_mo_active, core_energy_active )
    energy_check = 0.0
    for p in 0...mo.nactive
        for q in 0..p
            pq = pair( p + 1, q + 1 ) 
            energy_check += dns1[ pq - 1 ] * h_mo_active[ pq - 1 ]
            for r in 0..p
                for s in 0..r
                    rs = pair( r + 1, s + 1 )
                    if rs > pq
                        break
                    end
                    pqrs = pair( pq, rs )
                    energy_check += dns2[ pqrs - 1 ] * eris_mo_active[ pqrs - 1 ] 
                end
            end
        end
    end
    printf( "Energy calculated by density1, 2 and MO integrals  %20.15f ( check_density )\n", energy_check + core_energy_active )
end

def check_density_unfold( mo, dns1, dns2, h_mo_active, eris_mo_active, core_energy_active )
    energy_check = 0.0
    for p in 0...mo.nactive
        for q in 0...mo.nactive
            pq = pair( p + 1, q + 1 ) 
            energy_check += dns1[ pq - 1 ] * h_mo_active[ pq - 1 ]
            for r in 0... mo.nactive
                for s in 0...mo.nactive
                    rs = pair( r + 1, s + 1 )
                    pqrs = pair( pq, rs )
                    energy_check += dns2[ pqrs - 1 ] * eris_mo_active[ pqrs - 1 ] 
                end
            end
        end
    end
    printf( "Energy calculated by density1, 2 and MO integrals  %20.15f ( check_density )\n", energy_check + core_energy_active )
end

def mkfock( mo, h_mo, eris_mo, dns1 )
# Fock(inactive)...(2.13a), Fock(active)...(2.13b) の生成
    fi = h_mo.clone
    fa = GSL::Vector.alloc( fi.size )
    pq = 0
    for p in 0...mo.norb
        for q in 0..p
            pq += 1
            for k in mo.core_range
                kk = pair( k + 1, k + 1 )
                pk = pair( p + 1, k + 1 )
                qk = pair( q + 1, k + 1 )
                fi[ pq - 1 ] += 2.0 * eris_mo[ pair( pq, kk ) - 1 ] - eris_mo[ pair( pk, qk ) - 1 ] 
            end

            for v in mo.act_range
#                for x in mo.act_range.first..v
                for x in mo.act_range
                    vx_cas = pair( v - mo.ncore + 1, x - mo.ncore + 1 )
                    vx = pair( v + 1, x + 1 )
                    pv = pair( p + 1, v + 1 )
                    qx = pair( q + 1, x + 1 )
#                   printf( "pq, vx, pv, qx : %5d %5d %5d %5d \n", pq, vx, pv, qx )
#                   printf( "pqvx, pvqx : %5d %5d \n", pair(pq, vx), pair(pv, qx) )
                    fa[ pq - 1 ] += dns1[ vx_cas - 1 ] * ( eris_mo[ pair( pq, vx ) - 1 ] - 0.5 * eris_mo[ pair( pv, qx ) - 1 ] )
                end
            end
        end
    end
    f = fi + fa
#    EmolUtil::print_ary_tri(fi, "FockI", 10, "%12.8f")
#    EmolUtil::print_ary_tri(fa, "FockA", 10, "%12.8f")
    EmolUtil::print_ary_tri(f, "Fock", 10, "%12.6f")
    return fi, fa, f
end

def mkBLB( mo, dns1, dns2, fi, fa, f, eris_mo )
# Roos, Siegbahn
# dns1, dns2 ... unfolded density matrix
    printf( "Number of orbital rotations: %5d ( %5d, %5d, %5d ) \n", mo.nrotations, mo.ncore_act, mo.ncore_ext, mo.nact_ext )
 
    it = 0   #   it = 0 ... reference CAS state
    sx_addr = []
    sx_addr.push( 1 )                             #  core -> act ( i->t )
    sx_addr.push( sx_addr[ 0 ] + mo.ncore_act )   #  core -> ext ( i->a ) 
    sx_addr.push( sx_addr[ 1 ] + mo.ncore_ext )   #  act  -> ext ( t->a )
    w_sx = GSL::Matrix.alloc( mo.nrotations + 1, mo.nrotations + 1 )
    w_sx[ 0, 0 ] = 0.0

# <0| H | i -> t >
#    p "<0| H | i -> t >"
    for i in mo.core_range
        for t in mo.act_range
#            printf( "%3d -> %3d : f %12.6f *** \n", i+1, t+1, f[ pair( i + 1, t + 1 ) - 1 ] )
            it += 1
            sum = 2.0 * fa[ pair( i + 1, t + 1 ) - 1 ]
#            printf( "1st term :  %12.6f\n", sum )    
            for u in mo.act_range
                tu_cas = pair( u - mo.ncore + 1, t - mo.ncore + 1 )
                iu = pair( i + 1, u + 1 )
                if t == u
                    sum += ( 2.0 - dns1[ tu_cas - 1 ] ) * fi[ iu - 1 ]
                else
                    sum += - dns1[ tu_cas - 1 ] * fi[ iu - 1 ]
                end
            end
#            printf( "1st + 2nd terms :  %12.6f\n", sum )    

            for u in mo.act_range
                iu = pair( i + 1, u + 1 )
                tu_cas = pair( u - mo.ncore + 1, t - mo.ncore + 1 )
                for v in  mo.act_range
                    for x in mo.act_range
                        vx = pair( v + 1, x + 1 )
                        vx_cas = pair( v - mo.ncore + 1, x - mo.ncore + 1 )
                        sum -= 2.0 * dns2[ pair( tu_cas, vx_cas ) - 1 ] * eris_mo[ pair( iu, vx ) - 1 ]
                    end
                end
            end
#            printf( "1st + 2nd + 3rd terms :  %12.6f\n", sum )    

            mt = 1.0 / Math::sqrt(2.0 - dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1) - 1 ])
            w_sx[ 0, it ] = sum * mt
            w_sx[ it, 0 ] = w_sx[ 0, it ]
#            printf( "0 - %3d :  %14.6e \n", it, w_sx[ 0, it ] )
#            printf( "%3d -> %3d :  %14.6e ( %14.6e )\n", i + 1, t + 1, w_sx[ 0, it ], f[ pair( i + 1, t + 1 ) - 1 ] )
        end
    end

# <0| H | i -> a >
#    p "<0| H | i -> a >"
    for i in mo.core_range
        for a in mo.ext_range
            it += 1
            w_sx[ 0, it ] = Math::sqrt( 2.0 ) * f[ pair( i + 1, a + 1 ) - 1 ]
            w_sx[ it, 0 ] = w_sx[ 0, it ]
#            printf( "0 - %3d :  %14.6e \n", it, w_sx[ 0, it ] )
        end
    end

# <0| H | t -> a >
#    p "<0| H | t -> a >"
    for t in mo.act_range
        for a in mo.ext_range
#           printf( "%3d -> %3d :  *** \n", t+1, a+1 )
            it += 1
            sum = 0.0
            for u in mo.act_range
                au = pair( u + 1, a + 1 )
                tu_cas = pair( u - mo.ncore + 1, t - mo.ncore + 1 )
                sum += fi[ au - 1 ] * dns1[ tu_cas - 1 ]
            end

            for u in mo.act_range
                au = pair( u + 1, a + 1 )
                tu_cas = pair( u - mo.ncore + 1, t - mo.ncore + 1 )
                for v in  mo.act_range
                    for x in mo.act_range
                        vx = pair( v + 1, x + 1 )
                        vx_cas = pair( v - mo.ncore + 1, x - mo.ncore + 1 )
                        sum += 2.0 * dns2[ pair( tu_cas, vx_cas ) - 1 ] * eris_mo[ pair( au, vx ) - 1 ]
                    end
                end
            end
            w_sx[ 0, it ] = sum / Math::sqrt( dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1) - 1  ] )
            w_sx[ it, 0 ] = w_sx[ 0, it ]
#            printf( "0 - %3d :  %14.6e \n", it, w_sx[ 0, it ] )
        end
    end

# < i -> t | H | j -> u > (4a)
    it = sx_addr[ 0 ] - 1 #   it = 0
    for i in mo.core_range
        for t in mo.act_range
            nt = dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
            root_mt = Math::sqrt( 2.0 - nt )
            it += 1
            ju = sx_addr[ 0 ] - 1 #  ju = 0
            for j in mo.core_range
                for u in mo.act_range
                    ju += 1
                    if it < ju
                        break   #####################################
                    end
                    ij = pair( i + 1, j + 1 ) - 1
                    tu = pair( t + 1, u + 1 ) - 1
                    nu = dns1[ pair( u - mo.ncore + 1, u - mo.ncore + 1 ) - 1 ]
                    root_mu = Math::sqrt( 2.0 - nu )
                    tmp = 0.0
                    if ( i == j )
                        if ( t == u )
                            tmp = f[ tu ] - f[ ij ]
#                            tmp = f[ tu ] - f[ ij ] + 0.5 * nt * eris_mo[ pair( tu + 1, tu + 1 ) - 1 ] 
#                            printf( "tu  %3d  f[tu]  %14.6f  ij  %3d  f[ij]  %14.6f  nt %14.6f \n", tu, f[tu], ij, f[ij], nt )
#                            printf( "pair( tu, tu )  %5d  eris %14.6f \n", pair( tu + 1, tu + 1 ), eris_mo[ pair(tu+1,tu+1) - 1 ] )
#                            printf( "tmp %14.6f \n", tmp )
                        else
                            tt = pair( t + 1, t + 1 ) - 1
                            uu = pair( u + 1, u + 1 ) - 1
                            dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
                            tmp = root_mt * root_mu * 0.5 * f[ tu ]
#                            tmp = root_mt * root_mu * ( 0.5 * f[ tu ] + 0.25 * ( nt * eris_mo[ pair( tt + 1, tu + 1 ) - 1 ] + nu * eris_mo[ pair( tu + 1, uu + 1 ) - 1 ] ) )
                        end
                    elsif ( t == u ) 
                        tmp = - f[ ij ]
                    end
                    w_sx[ it, ju ] = tmp
                    w_sx[ ju, it ] = tmp                  
#                    printf( "it - ju  %3d  %3d :  %14.6e \n", it, ju, w_sx[ it,  ju ] )
                end
            end
        end
    end
# < i -> t | H | j -> a > (4b)
    it = sx_addr[ 0 ] - 1 #   it = 0
    for i in mo.core_range
        for t in mo.act_range
            nt = dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
            root_mt2 = Math::sqrt( ( 2.0 - nt ) / 2.0 )
            tt = pair( t + 1, t + 1 ) - 1
            it += 1
            ja = sx_addr[ 1 ] - 1
            for j in mo.core_range
                for a in mo.ext_range
                    ja += 1
                    at = pair( a + 1, t + 1 ) - 1
                    tmp = 0.0
#                    tmp = 0.5 * nt * eris_mo[ pair( at + 1, tt + 1 ) - 1 ]
                    if ( i == j )
                        tmp += f[ at ]
                    end
                    w_sx[ it, ja ] = root_mt2 * tmp
                    w_sx[ ja, it ] = w_sx[ it, ja ]
#                    printf( "it - ja  %3d  %3d :  %14.6e \n", it, ja, w_sx[ it, ja ] )
                end
            end
        end
    end
# < i -> t | H | j -> a > = 0 (4c)
# < i -> t | H | j -> u > (4d)
    ia = sx_addr[ 1 ] - 1
    for i in mo.core_range
        for a in mo.ext_range
            ia += 1
            jb = sx_addr[ 1 ] - 1
            for j in mo.core_range
                for b in mo.ext_range
                    jb += 1
                    if ( ia < jb )
                        break
                    end
                    ij = pair( i + 1, j + 1 ) - 1
                    ab = pair( a + 1, b + 1 ) - 1
                    tmp = 0.0
                    if ( i == j )
                        tmp = f[ ab ]
                    end
                    if ( a == b )
                        tmp += - f[ ij ]
                    end
                    w_sx[ ia, jb ] = tmp
                    w_sx[ jb, ia ] = tmp
#                    printf( "ia - jb  %3d  %3d :  %14.6e \n", ia, jb, w_sx[ ia, jb ] )
                end
            end
        end
    end 
# < i -> a | H | t -> b > (4e)
    ia = sx_addr[ 1 ] - 1
    for i in mo.core_range
        for a in mo.ext_range
            ia += 1
            tb = sx_addr[ 2 ] - 1
            for t in mo.act_range
                nt = dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
                mt = 2.0 - nt
                it = pair( i + 1, t + 1 ) - 1
                tt = pair( t + 1, t + 1 ) - 1
                for b in mo.ext_range
                    tb += 1
#                printf( "it, ja, %5d, %5d\n", it, ja )
                    tmp = 0.0
                    if ( a == b )
                        tmp = - Math::sqrt( nt * 0.5 ) * f[ it ]
#                        tmp = - Math::sqrt( nt * 0.5 ) * f[ it ] + mt * 0.5 / Math::sqrt( 2 * nt ) * eris_mo[ pair( it + 1, tt + 1 ) - 1 ]
                    end
                    w_sx[ tb, ia ] = tmp
                    w_sx[ ia, tb ] = tmp
                end
            end
        end
    end 

# < t -> a | H | u -> b > (4f)
    ta = sx_addr[ 2 ] - 1
    for t in mo.act_range
        tt = pair( t + 1, t + 1 ) - 1
        for a in mo.ext_range
            nt = dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
            mt = 2.0 - nt
            ta += 1
            ub = sx_addr[ 2 ] - 1
            for u in mo.act_range
                uu = pair( u + 1, u + 1 ) - 1
                tu = pair( t + 1, u + 1 ) - 1
                for b in mo.ext_range
                    ub += 1
                    if ( ta < ub )
                        break
                    end
                    ab = pair( a + 1, b + 1 ) - 1
                    nu = dns1[ pair( u - mo.ncore + 1, u - mo.ncore + 1 ) - 1 ]
                    mu = 2.0 - nu
                    tmp = 0.0
                    if ( a == b )
                        if ( t == u )
                            tmp = f[ ab ] - f[ tu ]
#                            tmp = f[ ab ] - f[ tu ] + 0.5 * mt * eris_mo[ pair( tu + 1, tu + 1) - 1 ]
                        else
                            tmp = Math::sqrt( nt * nu ) * ( - 0.5 * f[ tu ] )
#                            tmp = Math::sqrt( nt * nu ) * ( - 0.5 * f[ tu ] + 0.25 * ( mt * eris_mo[ pair( tt + 1, tu + 1) - 1 ] + mu * eris_mo[ pair( tu + 1, uu + 1 ) - 1 ] ) )
                        end
                    elsif ( t == u ) 
                        tmp = f[ ab ]
                    end
#
                    w_sx[ ta, ub ] = tmp
                    w_sx[ ub, ta ] = tmp
#                    printf( "ta - ub  %3d  %3d :  %14.6e \n", ta, ub, w_sx[ ta, ub ] )
                end
            end
        end
    end 

#   Level shift on Diagonal elements following the algorithm of Jason2 by S. Yamamoto
    thresh_sx = 1.9; dsx_shift = 0.1
    sx_diag_min = thresh_sx
    for p in 1...mo.nrotations + 1 
        if sx_diag_min > w_sx[ p, p ]
            sx_diag_min = w_sx[ p, p ]
        end
    end
    shift = thresh_sx - sx_diag_min + dsx_shift
    printf( "shift = %14.8f\n", shift )
#
#    EmolUtil::print_2dim_ary_tri(w_sx.to_a, "BLB", 10, "%10.6f")
#    printf( "Debug print on Diagonal elements of w_sx\n" )
    for p in 1...mo.nrotations + 1 
        w_sx[ p, p ] = w_sx[ p, p ] + shift
#        printf( "%5d   %14.8f\n", p, w_sx[ p, p ] )
    end

#    printf( "Debug print on Off-diagonal elements ( column = 1 ) of w_sx\n" )
#    for p in 0...mo.nrotations + 1 
#        printf( "%14.8f, \n", w_sx[ 1, p ] )
#    end

#    printf( "Debug print on Off-diagonal elements ( column = 13 ) of w_sx\n" )
#    for p in 0...mo.nrotations + 1 
#        printf( "%14.8f, \n", w_sx[ 13, p ] )
#    end

    eval, evec = GSL::Eigen::symmv( w_sx )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )
#    printf( "\nEnergy lowering from BLB = %20.16f \n\n", eval[ 0 ] )
    sol = 0
#    for i in 1...eval.size
#        if evec[ 0, sol ].abs < evec[ 0, i ].abs
#            sol = i
#        end
#    end
    printf( "\nEnergy lowering from BLB\n" )
    printf( "%5d   %12.6f   %12.6f\n", sol, evec[ 0, sol ], eval[ sol ] )

#                   BLB 行列の Liu によるアルゴリズムによる対角化
#    printf( "BLB matrix diagnalization by Liu\n" )
#    sx = SX.new( w_sx )
#    thresh = EmolConsts::THR_CONV_WF
#    eval, evec = Liu.new( sx, 1, thresh ).solve
#    printf( "\nEnergy lowering from BLB = %20.16f \n\n", eval[ 0 ] )
#
    w = 0.0
    for p in 0...mo.nrotations + 1 
        w += evec[ p, sol ] * evec[ p, sol ]
#        if ( evec[ p, sol ].abs > 0.01 )
            printf( " %5d  %16.10f \n", p, evec[ p, sol ] )
#        end
    end
    printf( "  norm  %16.8e \n", Math::sqrt( w ) )
    return sol, eval, evec
end

def mkdensity_SX( mo, cao, sol, cvec, dns1, dns2 )
    a = GSL::Vector.alloc( mo.nrotations + 1 )
    a0 = cvec[ 0, sol ]
#                          (1.B.2)
    ait = GSL::Matrix.alloc( mo.ncore, mo.nactive )
    aia = GSL::Matrix.alloc( mo.ncore, mo.nexternal )
    ata = GSL::Matrix.alloc( mo.nactive, mo.nexternal )
#
    core_range = 0...mo.ncore
    act_range = 0...mo.nactive
    ext_range = 0...mo.nexternal
    iadr = 0
    tadr = mo.ncore
    aadr = tadr + mo.nactive

    nt = []; mt = []
    for t in act_range
        tmp = dns1[ pair( t + 1, t + 1 ) - 1 ]
        mt.push( 1.0 / Math::sqrt( 2.0 - tmp ) )
        nt.push( 1.0 / Math::sqrt( tmp ) )
    end

    k = 0
    for i in core_range
        for t in act_range
            k += 1
#            ait[ i , t ] = cvec[ k, sol ] / Math::sqrt( 2.0 - dns1[ pair( t + 1, t + 1 ) - 1 ] )
            ait[ i , t ] = cvec[ k, sol ] * mt[ t ]
        end
    end
#
    sqrt2 = Math::sqrt( 2.0 )                               
    for i in core_range
        for a in ext_range
            k += 1
            aia[ i, a ] = cvec[ k, sol ] / sqrt2
        end
    end
#
    for t in act_range
        for a in ext_range
            k += 1
#            ata[ t,  a ] = cvec[ k, sol ] / Math::sqrt( dns1[ pair( t + 1, t + 1 ) - 1 ] )
            ata[ t,  a ] = cvec[ k, sol ] * nt[ t ]
        end
    end

    EmolUtil::print_2dim_ary(ait.to_a, "ait", 10, "%14.8f")
    EmolUtil::print_2dim_ary(aia.to_a, "aia", 10, "%14.8f")
    EmolUtil::print_2dim_ary(ata.to_a, "ata", 10, "%14.8f")

    dsx = GSL::Matrix.alloc( mo.norb, mo.norb )
#    a0 = a0 / Math::sqrt( 1.00003777 )    
#    printf( "a0 = %f18.10 ", a0 )
#    p ""
#                            (1.B.3a)
    for i in core_range
        for j in 0...i + 1
            sum = 0.0
            if ( i == j )
                sum = 2.0
            end
            psum = 0.0
            for a in act_range
                psum += aia[ i, a ] * aia[ j, a ]
            end
            sum -= 2.0 * psum
            for t in act_range
                for u in act_range
                    tu_cas = pair( t + 1, u + 1 )
                    if t == u
                        sum += - ait[ i, t ] * ait[ j, u ] * ( 2.0 - dns1[ tu_cas - 1 ] )
                    else
                        sum += + ait[ i, t ] * ait[ j, u ] * dns1[ tu_cas - 1 ]
                    end
                end
            end
#            printf( "i, j, dsx  %5d  %5d  %10.5f\n", i, j, sum )
            dsx[ iadr + i, iadr + j ] = sum; dsx[ iadr + j, iadr + i ] = sum
        end
    end
#                            (1.B.3b)
    for i in core_range
        for t in act_range
            sum = 0.0
            for u in act_range
                tu_cas = pair( t + 1, u + 1 )
                if t == u
                    sum += ait[ i, u ] * ( 2.0 - dns1[ tu_cas - 1 ] )
                else
                    sum += - ait[ i, u ] * dns1[ tu_cas - 1 ]
                end
            end
            sum *= a0
            for u in act_range
                for a in ext_range
                    sum += - aia[ i, a ] * ata[ u, a ] * dns1[ tu_cas - 1 ]
                end
            end   
            dsx[ iadr + i, tadr + t ] = sum; dsx[ tadr + t, iadr + i ] = sum
#            printf( "i, t, dsx  %5d  %5d  %10.5f\n", i, t, sum )
        end
    end
#                            (1.B.3c)
    for i in core_range
        for a in ext_range
            dsx[ iadr + i, aadr + a ] = 2.0 * a0 * aia[ i, a ]
            dsx[ aadr + a, iadr + i ] = dsx[ iadr + i, aadr + a ]
#            printf( "i, a, dsx  %5d  %5d  %10.5f\n", i, a, dsx[ a + aadr, i + iadr ] )
        end
    end
#                            (1.B.3d)
    for t in act_range
        for u in 0...t + 1
            tu_cas = pair( t + 1, u + 1 )
            sum = dns1[ tu_cas - 1 ]
#            sum = d0
            for v in act_range
                for x in act_range
                    psum = 0.0
                    for i in core_range
########               3d 式の符号に誤りがあるので、JASON2 に従って修正した
#                        psum += ait[ i, v ] * ait[ i, x ]     # original
                         psum -= ait[ i, v ] * ait[ i, x ]     # JASON2
                    end
                    for a in ext_range
                         psum += ata[ v, a ] * ata[ x, a ]
                    end
                    vx = pair( v + 1, x + 1 )
                    tu = pair( t + 1, u + 1 )
                    sum += ( 2.0 * dns2[ pair( vx, tu ) - 1 ] - dns1[ vx - 1 ] * dns1[ tu - 1 ] ) * psum
                end
            end

            for i in core_range
                for v in act_range
                    tv = pair( t + 1, v + 1 )
                    uv = pair( u + 1, v + 1 )
                    if ( t == v )
                          sum += ait[ i, v ] * ait[ i, u ] * ( 1.0 - dns1[ tv - 1 ] )
                    else
                          sum += ait[ i, v ] * ait[ i, u ] * ( - dns1[ tv - 1 ] )
                    end
                    if ( u == v )
                        sum += ait[ i, v ] * ait[ i, t ] * ( 1.0 - dns1[ uv - 1 ] )
                    else
                        sum += ait[ i, v ] * ait[ i, t ] * ( - dns1[ uv - 1 ] )
                    end
                end
            end
            dsx[ tadr + t, tadr + u] = sum; dsx[ tadr + u, tadr + t ] = sum
#            printf( "t, u, dsx  %5d  %5d  %10.5f\n", t, u, sum )
        end
    end
#                            (1.B.3e)
    for t in act_range
        for a in ext_range
            sum = 0.0
            for u in act_range
                tu = pair( t + 1, u + 1 )
                sum += ata[ u, a ] * dns1[ tu - 1 ]
            end
            sum *= a0
            for i in core_range
                for u in act_range
                    tu = pair( t + 1, u + 1 )
                    if ( t == u )
                        sum += aia[ i, a ] * ait[ i, u ] * ( 2.0 - dns1[ tu - 1 ] )
                    else
                        sum += aia[ i, a ] * ait[ i, u ] * (     - dns1[ tu - 1 ] )
                    end
                end
            end
            dsx[ tadr + t, aadr + a ] = sum; dsx[ aadr + a, tadr + t ] = sum
        end
    end
#                            (1.B.3f)
    for a in ext_range
        for b in 0...a + 1
            sum = 0.0
            for i in core_range
                sum += aia[ i, a ] * aia[ i, b ]
            end
            sum *= 2.0
            for t in act_range
                for u in act_range
                    tu = pair( t + 1, u + 1 )
                    sum += ata[ t, a ] * ata[ u, b ] * dns1[ tu - 1 ]
                end
            end
            dsx[ aadr + a, aadr + b ] = sum; dsx[ aadr + b, aadr + a ] = sum
#            printf( "a, b, dsx  %5d  %5d  %10.5f\n", a, b, sum )
        end
    end
    EmolUtil::print_2dim_ary(dsx.to_a, "SXdensity", 7, "%16.10f")

#   check trace
    printf( "\n" )
    sum = 0.0
    for k in 0...mo.norb
        sum += dsx[ k, k ]
    end
    printf( "\n" )
    printf( "check trace of dsx : %12.8e\n", sum )    

    eval, evec = GSL::Eigen::symmv( dsx )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_DESC )
    printf( "check sum of occ : %12.8e\n", eval.sum )    
#   EmolUtil::print_2dim_ary(evec.to_a, "Natural orbitals", 10, "%14.8f")
#   EmolUtil::print_2dim_ary_with_value(evec.to_a, eval.to_a, "new MO", 10, "%14.8f")
    new_cmo = GSL::Matrix.alloc( evec.size[ 0 ], evec.size[ 0 ] )
    occ = GSL::Vector.alloc( evec.size[ 0 ] )
#                              core, active mo は previous mo の character 維持、それ以外の mo は占有数順
    flag = Array.new( mo.norb )
    for k in 0...mo.norb
        flag[ k ] = true
    end
#    for k in mo.core_range
    for k in mo.core_range.first...mo.act_range.last
        from = find_max( evec, k )
        occ[ k ] = eval[ from ]
        flag[ from ] = false
        for i in 0...evec.size[ 0 ]
            new_cmo[ i, k ] = evec[ i, from ]
        end
    end
    from = 0
#    for k in mo.act_range.first...mo.norb
    for k in mo.ext_range.first...mo.norb
        until flag[ from ]
            from += 1
        end
        occ[ k ] = eval[ from ]
        for i in 0...evec.size[ 0 ]
            new_cmo[ i, k ] = evec[ i, from ]
        end
        from += 1
    end
#    new_cmo = evec
#    occ = eval
    EmolUtil::print_2dim_ary_with_value(new_cmo.to_a, occ, "Natural Orbitals of superCI", 10, "%14.8f")

#                     no in ao representation
        nbf = cao.size[ 0 ]
        cao_active = GSL::Matrix.alloc(nbf, mo.norb)
        n_frozen = cao.size[ 1 ] - mo.norb
        for i in 0..mo.norb - 1
            for p in 0..nbf - 1
                cao_active[ p, i ] = cao[ p, n_frozen + i ]
            end
        end
#   EmolUtil::print_2dim_ary( cao_active.to_a, "Active MOs", 8, "%12.6f" )
        no_ao = cao_active * new_cmo
#   EmolUtil::print_2dim_ary( no_ao.to_a, "NOs-AOs", 8, "%12.6f" )

#                     add frozen core 
        if n_frozen > 0
            cao_new = cao.submatrix(nil, 0, n_frozen).horzcat( no_ao )
            fc = GSL::Vector.alloc( n_frozen )
            fc.set_all( 2.0 )
            occ = fc.connect( occ ).to_a
        else
            cao_new = no_ao
        end
  EmolUtil::print_2dim_ary_with_value( cao_new.to_a, occ, "Natural Orbitals of SuperCI in AO representation", 8, "%14.8f" )
############
    return cao_new
end

def find_max( a, k )
    max_clm = -1;     max_value = 1.0e-20
    for i in 0...a.size[ 0 ]
        if a[ k, i ].abs > max_value
            max_clm = i; max_value = a[ k, i ].abs
        end 
    end
    max_clm
end

def trans( nob, cmo, h_mo, eris_mo )
    cpq = GSL::Vector.alloc( nob * nob )
    new_h_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
    new_cmo_1dim = GSL::Vector.alloc( nob * nob )
    ip = -1
    for i in 0...nob
        for p in 0...nob
            ip += 1
            new_cmo_1dim[ ip ] = cmo[ p, i ] 
        end
    end
    Inttrans::trnpp( 1, 1, nob, nob, h_mo, new_h_mo, 1, new_cmo_1dim, cpq )
    #EmolUtil::print_ary_tri(new_h_mo.to_a, "H matrix", 10, "%12.6f")

    mi = 1
    nbf = nob
    mpq = nob * ( nob + 1 ) / 2
    npq = mpq
    eris_work = GSL::Vector.alloc( npq * npq )
    pqrs = -1
    for pq in 0...npq
        for rs in 0...pq + 1
            pqrs += 1
            eris_work[ pq * npq + rs ] = eris_mo[ pqrs ]
            eris_work[ rs * npq + pq ] = eris_mo[ pqrs ]
        end
    end

    for n in 1..npq
        Inttrans::trnpp(mi, mi, nbf, nob, eris_work, eris_work, 1, new_cmo_1dim, cpq)
        mi += npq
    end
    Inttrans::trnsps(npq, mpq, eris_work)

    res = []
    buf = GSL::Vector.alloc( npq )
    new_eris_mo = GSL::Vector.alloc( mpq*(mpq + 1)/2 + 1 )
    m = 1; kin = 0
    for n in 1..mpq
        Inttrans::trnrr(m, n, nbf, nob, eris_work, buf, new_cmo_1dim, cpq)
        m += npq
        for k in 0...n
            new_eris_mo[ kin + k ] = buf[ k ]
        end
        kin += n
    end
#    EmolUtil::print_ary_tri(new_eris_mo.to_a, "New Eris matrix", 10, "%12.6f")
    return new_h_mo, new_eris_mo
end

# 軌道、状態情報の取得
load "Data_molinfo"
molinfo = Data_molinfo::molinfo_data
nel = molinfo.nel
spin = molinfo.spin

n_frozen = molinfo.n_frozen
n_core = molinfo.n_core
n_active = molinfo.n_active
n_external = molinfo.n_external
nob = n_core + n_active + n_external
nve = nel - 2 * ( n_frozen + n_core )
mo = MO_classification_mcscf.new( n_core, n_active, n_external )
#val_range = ( n_frozen + n_core )...( n_frozen + n_core + n_active )
val_range = mo.act_range 

# CAS エネルギー表式の生成
br = BrooksCases.new( Second_order_CI.new( nve, 0, 0, nve, 0 ) )
cas_expr = br.expr
cas_expr_one = br.expr_one
#br.expr.show
#br.expr_one.show

load "Data_rhf"
cao = Data_rhf::cao_data
#trnint = Transformation.new( [ "Data_molinfo", "Data_molint", "Data_rhf" ] )
# 分子積分
#     core_energy には frozen_core のエネルギーを含む
#     one electron integrals (tri) 
#     frozen_core との相互作用を含む
#     two electron integrals (tri, tri, tri)

it = -1
max_iter = 40
energy = []
lowering = []
blb = []
lowering_BLB = 0.0
    printf( "===    it         Energy              E(it)-E(it-1)       BLBlowering\n" )
while it < max_iter
    molint = Transformation_mcscf.new( [ "Data_molinfo", "Data_molint", cao ] )
    core_energy = molint.core_energy
    h_mo = molint.h_mo
    eris_mo = molint.eris_mo

    it += 1
    core_energy_active, h_mo_active, eris_mo_active = extract_active_integral( mo, core_energy, h_mo, eris_mo )
#                      CASCI の実行 
    cas_hmat = H_matrix_casscf.new( molinfo, cas_expr, h_mo_active, eris_mo_active, core_energy_active )
    cf = cas_hmat.eigen_vector
    energy.push( cas_hmat.eigen_value[ 0 ] )
    if it > 0
        lowering.push( energy[ it ] - energy[ it - 1 ] )
    else
        lowering.push( energy[ it ] )
    end
#                      Density の生成
    dns1, dns2, dns1_unfold, dns2_unfold = mkdensity( mo, cas_expr, cf )
    check_density( mo, dns1, dns2, h_mo_active, eris_mo_active, core_energy_active )
    check_density_unfold( mo, dns1_unfold, dns2_unfold, h_mo_active, eris_mo_active, core_energy_active )
#                      Fock の生成
#    fi, fa, f = mkfock( mo, h_mo, eris_mo, dns1 )
    fi, fa, f = mkfock( mo, h_mo, eris_mo, dns1_unfold )
#                      BLB matrix elements の生成
#    eval, evec = mkBLB( mo, dns1, dns2, fi, fa, f, eris_mo )
    sol, eval, evec = mkBLB( mo, dns1_unfold, dns2_unfold, fi, fa, f, eris_mo )
    lowering_BLB = eval[ sol ]
    blb.push( lowering_BLB )
    printf( "=== %5d    %16.8f    %16.8f   %16.8f\n", it, energy[ it ], lowering[ it ], blb[ it ] )
    cao = mkdensity_SX( mo, cao, sol, evec, dns1_unfold, dns2_unfold )
#    EmolUtil::print_2dim_ary(new_cao.to_a, "debgug print on cao matrix", 10, "%12.6f")
end




