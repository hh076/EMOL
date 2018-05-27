# -*- coding: utf-8 -*-
#!/usr/bin/ruby
#printf "ruby version: %s\n", RUBY_VERSION
if ( RUBY_VERSION.to_f <= 1.8 ) then
  def require_relative( str )
    str2 = "./" + str
    require str2
  end
end

require_relative './emol_nl'
require_relative './emol_drt_soci'
require_relative './emol_molinfo'
require_relative './emol_trans'
require_relative './emol_expression'
require_relative './emol_hmatrix'
require_relative './emol_liu'
require_relative './emol_futil'
require_relative './emol_util'
require_relative './emol_consts'

def pair( p, q )
    if ( p > q )
        p * ( p - 1 ) / 2 + q
    else
        q * ( q - 1 ) / 2 + p
    end
end

class H_matrix_casscf < H_matrix
    def initialize( molinfo, expr, h_mo, eris_mo, core_energy )
        if EmolConsts::PRINT_LEVEL > 0 then
        	printf( "**** CASCI \n" )
	end
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
            if EmolConsts::PRINT_LEVEL > 0 then
                printf( "     %5d   %20.15f \n", k, @eigen_value[ k ] )
	    end
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
      self.show
    end

    def show
        orb_pair_tbl = GSL::Matrix.alloc( @norb, @norb )
        @orb_pair = []
        printf( "Number of MOs       %5d \n", @norb )
        printf( "   Inactive         %5d \n", @ncore )
        printf( "   Active           %5d \n", @nactive )
        printf( "   External         %5d \n\n", @nexternal )
        printf( "Non-Redundant Orbital Rotations\n" )
        it = 0
        for i in @core_range
            for t in @act_range
                it += 1
                orb_pair_tbl[ i, t ] = it
                orb_pair_tbl[ t, i ] = it
                @orb_pair.push( [i, t] )
            end
        end
        for i in @core_range
            for a in @ext_range
                it += 1
                orb_pair_tbl[ i, a ] = it
                orb_pair_tbl[ a, i ] = it
                @orb_pair.push( [i + 1, a + 1] )
            end
        end
        for t in @act_range
            for a in @ext_range
                it += 1
                orb_pair_tbl[ t, a ] = it
                orb_pair_tbl[ a, t ] = it
                @orb_pair.push( [t + 1, a + 1] )
            end
        end

#        @orb_pair.each_with_index do | x, i |
#            printf( "%5d ( %2d %2d ) \n", i, x[ 0 ], x[ 1 ]  )
#        end

        printf( "      " )
        for i in 0...@act_range.last
            printf(" %5d", i + 1 )
        end
        printf("\n")

        for i in 0...@norb
            printf(" %5d", i + 1 )
            for j in 0...@act_range.last
                printf(" %5d", orb_pair_tbl[ j, i ] )
            end
            printf("\n")
        end
    end

    def get_pair( it )
        @orb_pair[ it ]
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
        s = Data_int::s_data
        h = Data_int::h_data
        v = Data_int::v_data
        nz_g = Data_int::nz_g_data

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
        if EmolConsts::PRINT_LEVEL > 0 then
            print "\n**** Transformation\n"
            printf( "nbf, n_frozen, nob : %5d  %5d  %5d  \n", nbf, n_frozen, nob )
            printf( "frozen core energy : %14.8f \n", @core_energy) 
        end

        cao_trans = GSL::Matrix.alloc( nbf, nob )
        for i in 0...nob
            for j in 0...nbf
                cao_trans[ j, i ] = cao[ j, n_frozen + i ]
            end
        end
        s_work, t_work, v_work, h_work, eris_work, cao_1dim = two_one_dimCopy( nbf, nob, s, f, v, h, nz_g, cao_trans )
#           一電子積分は凍結殼との相互作用を含む Fock を変換する ( h ではなく f )
        cpq = GSL::Vector.alloc( nbf * nob )
        @s_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
        @h_mo = GSL::Vector.alloc( nob * ( nob + 1 ) / 2 )
        Inttrans::trnpp( 1, 1, nbf, nob, s_work, @s_mo, 1, cao_1dim, cpq )
        Inttrans::trnpp( 1, 1, nbf, nob, h_work, @h_mo, 1, cao_1dim, cpq )
        if EmolConsts::PRINT_LEVEL > 1 then
            EmolUtil::print_2dim_ary(cao_trans.to_a, "Transformation matrix", 8, 
                                     "%12.6f")
            EmolUtil::print_ary_tri(@s_mo.to_a, "Overlap matrix", 8, "%12.8f")
            EmolUtil::print_ary_tri(@h_mo.to_a, "H matrix", 8, "%12.8f")
        end

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
        if EmolConsts::PRINT_LEVEL > 2 then
            EmolUtil::print_ary_tri(eris_mo.to_a, "Eris matrix", 8, "%12.8f")
        end

        if EmolConsts::PRINT_LEVEL > 1 then
            printf( "Total energy calculated from mo integrals : %14.8f \n", 
            @core_energy + close_shell_energy(h_mo, eris_mo, molinfo.get_nel - n_frozen * 2) )
        end
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

def mkdensity_CAS( mo, cas_expr, cf )
    if EmolConsts::PRINT_LEVEL > 0 then
    	printf( "**** MkDensity_CAS \n" )
    end
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

    if EmolConsts::PRINT_LEVEL > 1 then           
        EmolUtil::print_ary_tri(dns1, "density1", 10, "%14.8f")
        EmolUtil::print_ary_tri(dns1_unfold, "density_unfold", 10, "%14.8f")
        EmolUtil::print_ary_tri(dns2, "density2", 10, "%14.8f")
        EmolUtil::print_ary_tri(dns2_unfold, "density2_unfold", 10, "%14.8f")
    end

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
#
    if EmolConsts::PRINT_LEVEL > 1 then           
    	printf( "**** Mkfock \n" )
    end
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
                for x in mo.act_range
                    vx_cas = pair( v - mo.ncore + 1, x - mo.ncore + 1 )
                    vx = pair( v + 1, x + 1 )
                    pv = pair( p + 1, v + 1 )
                    qx = pair( q + 1, x + 1 )
                    fa[ pq - 1 ] += dns1[ vx_cas - 1 ] * ( eris_mo[ pair( pq, vx ) - 1 ] - 0.5 * eris_mo[ pair( pv, qx ) - 1 ] )
                end
            end
        end
    end
    f = fi + fa

    if EmolConsts::PRINT_LEVEL > 1 then           
        EmolUtil::print_ary_tri(fi, "FockI", 8, "%12.8f")
        EmolUtil::print_ary_tri(fa, "FockA", 8, "%12.8f")
        EmolUtil::print_ary_tri(f, "Fock", 8, "%12.6f")
    end

    return fi, fa, f
end

def mkfock_general( mo, fi, f, eris_mo, dns1, dns2 )
#
#   Siegbahn(1981) eq(14)
#
#
    fg = GSL::Matrix.alloc( mo.ncore + mo.nactive, mo.norb )
    for i in mo.core_range
        for q in 0...mo.norb
            fg[ i, q ] = 2.0 * f[ pair( i + 1, q + 1 ) - 1 ]    
        end     
    end

    for t in mo.act_range
        for q in 0...mo.norb
            tmp = 0.0
            for u in mo.act_range
                tu_cas = pair( t - mo.ncore + 1, u - mo.ncore + 1 )
                qu = pair( q + 1, u + 1 )
                tmp += dns1[ tu_cas - 1 ] * fi[ qu - 1 ]
            end  
            tmp2 = 0.0        
            for u in mo.act_range
                tu_cas = pair( t - mo.ncore + 1, u - mo.ncore + 1 )
                qu = pair( q + 1, u + 1 )
                for v in mo.act_range
                    for x in mo.act_range
                        vx_cas = pair( v - mo.ncore + 1, x - mo.ncore + 1 )
                        vx = pair( v + 1, x + 1 )
                        tuvx_cas = pair( tu_cas, vx_cas )
                        quvx = pair( qu, vx )
                        tmp2 += dns2[ tuvx_cas - 1 ] * eris_mo[ quvx - 1 ]
                    end
                end    
            end
            fg[ t, q ] = tmp + 2.0 * tmp2
        end
    end
    if EmolConsts::PRINT_LEVEL > 1 then           
        EmolUtil::print_2dim_ary(fg.transpose.to_a, "FockGeneral", 8, "%12.8f")
    end
    return fg
end

def mkSX( mo, dns1, dns2, fi, fa, f, eris_mo, energy_shift )
# Roos, Siegbahn
# dns1, dns2 ... unfolded density matrix
#    printf( "Number of orbital rotations: %5d ( %5d, %5d, %5d ) \n", mo.nrotations, mo.ncore_act, mo.ncore_ext, mo.nact_ext )
#
   if EmolConsts::PRINT_LEVEL > 0 then           
   	printf( "**** MkSX \n" ) 
   end	
   it = 0   #   it = 0 ... reference CAS state
    sx_addr = []
    sx_addr.push( 1 )                             #  core -> act ( i->t )
    sx_addr.push( sx_addr[ 0 ] + mo.ncore_act )   #  core -> ext ( i->a ) 
    sx_addr.push( sx_addr[ 1 ] + mo.ncore_ext )   #  act  -> ext ( t->a )
    w_sx = GSL::Matrix.alloc( mo.nrotations + 1, mo.nrotations + 1 )
    w_sx[ 0, 0 ] = 0.0
    blb_max = EmolConsts::THR_ZERO

# <0| H | i -> t >
#    p "<0| H | i -> t >"
    for i in mo.core_range
        for t in mo.act_range
            it += 1
            sum = 2.0 * fa[ pair( i + 1, t + 1 ) - 1 ]
            for u in mo.act_range
                tu_cas = pair( u - mo.ncore + 1, t - mo.ncore + 1 )
                iu = pair( i + 1, u + 1 )
                if t == u
                    sum += ( 2.0 - dns1[ tu_cas - 1 ] ) * fi[ iu - 1 ]
                else
                    sum += - dns1[ tu_cas - 1 ] * fi[ iu - 1 ]
                end
            end

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

            mt = 1.0 / Math::sqrt(2.0 - dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1) - 1 ])
            w_sx[ 0, it ] = sum * mt
            w_sx[ it, 0 ] = w_sx[ 0, it ]
            if w_sx[ 0, it ].abs > blb_max
                blb_max = w_sx[ 0, it ].abs
            end
        end
    end

# <0| H | i -> a >
    for i in mo.core_range
        for a in mo.ext_range
            it += 1
            w_sx[ 0, it ] = Math::sqrt( 2.0 ) * f[ pair( i + 1, a + 1 ) - 1 ]
            w_sx[ it, 0 ] = w_sx[ 0, it ]
            if w_sx[ 0, it ].abs > blb_max
                blb_max = w_sx[ 0, it ].abs
            end
        end
    end

# <0| H | t -> a >
    for t in mo.act_range
        for a in mo.ext_range
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
            if w_sx[ 0, it ].abs > blb_max
                blb_max = w_sx[ 0, it ].abs
            end
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
                        else
                            tt = pair( t + 1, t + 1 ) - 1
                            uu = pair( u + 1, u + 1 ) - 1
                            dns1[ pair( t - mo.ncore + 1, t - mo.ncore + 1 ) - 1 ]
                            tmp = root_mt * root_mu * 0.5 * f[ tu ]
                        end
                    elsif ( t == u ) 
                        tmp = - f[ ij ]
                    end
                    w_sx[ it, ju ] = tmp
                    w_sx[ ju, it ] = tmp                  
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
                    if ( i == j )
                        tmp += f[ at ]
                    end
                    w_sx[ it, ja ] = root_mt2 * tmp
                    w_sx[ ja, it ] = w_sx[ it, ja ]
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
                    tmp = 0.0
                    if ( a == b )
                        tmp = - Math::sqrt( nt * 0.5 ) * f[ it ]
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
                        else
                            tmp = Math::sqrt( nt * nu ) * ( - 0.5 * f[ tu ] )
                        end
                    elsif ( t == u ) 
                        tmp = f[ ab ]
                    end
                    w_sx[ ta, ub ] = tmp
                    w_sx[ ub, ta ] = tmp
                end
            end
        end
    end 

#   Level shift on Diagonal elements following the algorithm of Jason2 by S. Yamamoto
#   energy_shift = thresh_sx + dsx_shift
    thresh_sx = energy_shift
#   thresh_sx = 1.9; dsx_shift = 0.6

#   sx_diag_min = thresh_sx
    sx_diag_min = thresh_sx - 0.6
    p_min = 0
    for p in 1...mo.nrotations + 1 
        if sx_diag_min > w_sx[ p, p ] then
            sx_diag_min = w_sx[ p, p ]
            p_min = p
        end
    end
#   shift = thresh_sx - sx_diag_min + dsx_shift
    shift = thresh_sx - sx_diag_min
    for p in 1...mo.nrotations + 1 
        w_sx[ p, p ] = w_sx[ p, p ] + shift
    end

    if EmolConsts::PRINT_LEVEL > 1 then
        printf( "\n*** Energy Shift for SX Diagonal Elements \n" )
        printf( "min diagonal SX value   = %14.8f \n", sx_diag_min )
        printf( "min diagonal SX element = %5d \n", p_min )
        printf( "Energy shift        = %14.8f\n", shift )
    end

    if EmolConsts::PRINT_LEVEL > 2 then
        EmolUtil::print_2dim_ary_tri(w_sx.to_a, "SX", 8, "%14.8f")
    end

    eval, evec = GSL::Eigen::symmv( w_sx )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )

#    printf( "\nEnergy lowering from SX = %20.16f \n\n", eval[ 0 ] )
     sol = 0
     if EmolConsts::PRINT_LEVEL > 0
         printf( "Energy lowering = %12.8f \n\n", eval[ sol ] )
         printf( "   Pair      Coef\n" )
         for i in 0...evec.size[ 0 ]
             if i == 0 then
                 printf( "   Ref     %12.8f \n", evec[ sol, i ] )
             elsif evec[ sol, i ].abs > 0.001
                 pair = mo.get_pair( i - 1 )
                 printf( " (%2d-%2d)   %12.8f \n", pair[ 0 ] + 1, pair[ 1 ] + 1, evec[ sol, i ] )
             end
         end
    end
    return blb_max, shift, sol, eval, evec
end

def mkdensity_SX( mo, cao, sol, cvec, dns1, dns2 )
    if EmolConsts::PRINT_LEVEL > 0
	    printf( "**** MkDensty_SX\n" )
    end
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
            ata[ t,  a ] = cvec[ k, sol ] * nt[ t ]
        end
    end

    if EmolConsts::PRINT_LEVEL > 2 then
        EmolUtil::print_2dim_ary(ait.to_a, "ait", 10, "%14.8f")
        EmolUtil::print_2dim_ary(aia.to_a, "aia", 10, "%14.8f")
        EmolUtil::print_2dim_ary(ata.to_a, "ata", 10, "%14.8f")
    end

    dsx = GSL::Matrix.alloc( mo.norb, mo.norb )
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
        end
    end
#                            (1.B.3c)
    for i in core_range
        for a in ext_range
            dsx[ iadr + i, aadr + a ] = 2.0 * a0 * aia[ i, a ]
            dsx[ aadr + a, iadr + i ] = dsx[ iadr + i, aadr + a ]
        end
    end
#                            (1.B.3d)
    for t in act_range
        for u in 0...t + 1
            tu_cas = pair( t + 1, u + 1 )
            sum = dns1[ tu_cas - 1 ]
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
        end
    end

#   check trace
    if EmolConsts::PRINT_LEVEL > 1 then
        EmolUtil::print_2dim_ary(dsx.to_a, "SXdensity", 8, "%14.10f")
        printf( "\n" )
        sum = 0.0
        for k in 0...mo.norb
            sum += dsx[ k, k ]
        end
        printf( "\n" )
        printf( "check trace of dsx : %12.8e\n", sum )
    end    

    eval, evec = GSL::Eigen::symmv( dsx )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_DESC )

    if EmolConsts::PRINT_LEVEL > 1 then
        printf( "check sum of occ : %12.8e\n", eval.sum )    
    end
#    EmolUtil::print_2dim_ary(evec.to_a, "Natural orbitals", 8, "%12.6f")
#    EmolUtil::print_2dim_ary_with_value(evec.to_a, eval.to_a, "new MO", 8, "%12.6f")
    new_cmo = GSL::Matrix.alloc( evec.size[ 0 ], evec.size[ 0 ] )
    occ = GSL::Vector.alloc( evec.size[ 0 ] )
#                              core, active mo は previous mo の character 維持、それ以外の mo は占有数順
    flag = Array.new( mo.norb )
    for k in 0...mo.norb
        flag[ k ] = true
    end
    for k in mo.core_range
#    for k in mo.core_range.first...mo.act_range.last
        from = find_max( evec, k )
        occ[ k ] = eval[ from ]
        flag[ from ] = false
        for i in 0...evec.size[ 0 ]
            new_cmo[ i, k ] = evec[ i, from ]
        end
    end
    from = 0
    for k in mo.act_range.first...mo.norb
#    for k in mo.ext_range.first...mo.norb
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
    sqcdf = 0.0
    for i in mo.core_range
        sum = 0.0
        for k in mo.core_range
            sum += new_cmo[ i, k ] * new_cmo[ i, k ] 
        end
        sqcdf += occ[ i ] * ( 1.0 - sum )
    end
    for i in mo.act_range
        sum = 0.0
        for k in mo.act_range
            sum += new_cmo[ i, k ] * new_cmo[ i, k ] 
        end
        sqcdf += occ[ i ] * ( 1.0 - sum )
    end
#
    if EmolConsts::PRINT_LEVEL > 2
        EmolUtil::print_2dim_ary_with_value(new_cmo.to_a, occ, "Natural Orbitals in MO representation", 8, "%14.8f")
    end
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

    if EmolConsts::PRINT_LEVEL > 0
        EmolUtil::print_2dim_ary_with_value( cao_new.to_a, occ, "Natural Orbitals in AO representation", 8, "%14.8f" )
############
    end
    return occ, cao_new, sqcdf
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

def driver( mode )
	printf( "\n\n\n" )
	printf( "************************************************************* \n" )
	printf( "    M C S C F                        coded by Takeshi Noro \n" )
	#printf( "       Roos et al ( Chem.Phys. vol.48, 157 (1980 ) \n"  )
	#printf( "       Siegbahn et al ( Physicaf Scripta vol.21, 323 (1980 ) \n" )
	printf( "************************************************************* \n\n" )

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

	istate        = molinfo.get_istate
	thr_sqcdf     = molinfo.get_thresh_mcscf
	max_iteration = molinfo.get_max_iteration
	energy_shift  = molinfo.get_energy_shift

        if ( mode == "START" ) then
		load "Data_rhf"
		cao = Data_rhf::cao_data
	else
		load "Data_mcscf"
		cao = Data_mcscf::cao_data
	end

	printf( "State Information \n" )
	printf( "number of electrons          %5d \n", nel )
	printf( "number of active electrons   %5d \n", nel - n_core * 2 )
	printf( "spin (S*2)                   %5d \n\n", spin )

	printf( "target state                 %5d \n", istate )
	printf( "max mcscf iterations         %5d \n", max_iteration )
	printf( "SQCDF threshold              %12.6e \n", thr_sqcdf )
	printf( "energy shift                 %12.6e \n\n", energy_shift )

	printf( "Orbital Information \n" )
	mo = MO_classification_mcscf.new( n_core, n_active, n_external )
	#val_range = ( n_frozen + n_core )...( n_frozen + n_core + n_active )
	val_range = mo.act_range 

	# CAS エネルギー表式の生成
	br = BrooksCases.new( Second_order_CI.new( nve, 0, 0, nve, 0 ) )
	cas_expr = br.expr
	cas_expr_one = br.expr_one
	#br.expr.show
	#br.expr_one.show

	#load "Data_rhf"
	#cao = Data_rhf::cao_data
	#load "Data_mcscf"
	#cao = Data_mcscf::cao_data

	if ( mode == "START" ) then
		EmolUtil::print_2dim_ary(cao.to_a, "Initial MO from RHF", 8, "%12.8f")
	else
		EmolUtil::print_2dim_ary(cao.to_a, "Initial MO from the previous MCSCF", 8, "%12.8f")
	end

	#trnint = Transformation.new( [ "Data_molinfo", "Data_int", "Data_rhf" ] )
	# 分子積分
	#     core_energy には frozen_core のエネルギーを含む
	#     one electron integrals (tri) 
	#     frozen_core との相互作用を含む
	#     two electron integrals (tri, tri, tri)

	printf( "\nMCSCF iteration starts >>>> \n\n" )

	it = -1
	sqcdf = 1.0
	energy = []
	lowering = []
	printf( "===    it      Energy         E(it)-E(it-1)     maxBLB         SQCDF            LevelShift  Step\n" )
	while ( it < max_iteration && sqcdf > thr_sqcdf ) 
	    molint = Transformation_mcscf.new( [ "Data_molinfo", "Data_int", cao ] )
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
	    dns1, dns2, dns1_unfold, dns2_unfold = mkdensity_CAS( mo, cas_expr, cf )
	    if EmolConsts::PRINT_LEVEL > 1 then
	        check_density( mo, dns1, dns2, h_mo_active, eris_mo_active, core_energy_active )
	        check_density_unfold( mo, dns1_unfold, dns2_unfold, h_mo_active, eris_mo_active, core_energy_active )
	    end
#                      Fock の生成
	    fi, fa, f = mkfock( mo, h_mo, eris_mo, dns1_unfold )
	    #    EmolUtil::print_2dim_ary(new_cao.to_a, "debgug print on cao matrix", 10, "%12.6f")
	#  if sqcdf > THR_SQCDF then
	#                      SX matrix elements の生成
            blb_max, shift, sol, eval, evec = mkSX( mo, dns1_unfold, dns2_unfold, fi, fa, f, eris_mo, energy_shift )
            occ, cao, sqcdf = mkdensity_SX( mo, cao, sol, evec, dns1_unfold, dns2_unfold )
            process = "SX"
	    printf( "=== %5d    %12.8f     %12.8f    %12.8f     %10.6e  %10.5f     %s \n", it, energy[ it ], lowering[ it ], blb_max, sqcdf, shift, process )
	end

	printf( "\n\n ***** Results of CASSCF\n\n" )
	eval = cas_hmat.eigen_value
	evec = cas_hmat.eigen_vector
	csf_list = cas_expr_one.csf_list
	bs_name = molinfo.get_bs_name
	print_casci_result( istate, eval, evec, csf_list, n_core, n_active, n_external )
	EmolUtil::print_2dim_ary_with_value_moname( cao.to_a, occ.to_a, bs_name, "Natural Orbitals and Occupation Numbers ", 6, "%12.8f" )

	return istate, eval, cao, sqcdf
end

def print_casci_result( istate, eval, evec, csf_list, n_core, n_active, n_external )
    printf( "TOTAL ENERGY =  %16.10f \n\n", eval[ istate - 1 ] )
    printf( "    CSF       COEF       ACTIVE \n")
    printf( "    ---    ----------    ") 
    n_active.times do | i |
        printf("--")
    end
    printf("  \n")

    vec_print = {}
    vec = evec.column( istate - 1 )
    vec.to_a.each_with_index do |x, i|
        if x.abs > 0.001
            vec_print.store( i,  x )
        end
    end

    vec_print.sort{ |(k1, v1), (k2, v2)| v2.abs <=> v1.abs }.each do | x |
        printf("  %5d  ", x[ 0 ] + 1)
        printf(" %10.6f   ", x[ 1 ])
        a = csf_list[ x[ 0 ] ]
 
        k = 0
        n_active.times do | i |
            printf("%2s", a[ k ] )
            k += 1
        end
        printf("  \n")
    end
end
