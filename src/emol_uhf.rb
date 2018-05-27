#!/usr/bin/ruby

require_relative 'test_rexml'
require_relative 'emol_util'
require_relative 'emol_consts'
require_relative 'emol_int'
require_relative 'emol_rhf'
require_relative 'molinfo'
require 'gsl'

def calc_gmat_nz_eris_uhf ( nz_eris, da, db, nbf )
  ga = GSL::Matrix.alloc(nbf, nbf)
  gb = GSL::Matrix.alloc(nbf, nbf)

  nz_eris.each do |xx|
    i = xx.get_i; j = xx.get_j; k = xx.get_k; l = xx.get_l; 
    x = xx.get_value

    if i == j then; x *= 0.5; end 
    if k == l then; x *= 0.5; end 
    if i == k && j == l then; x *= 0.5; end

    tmp = 4 * x * ( da [ k, l ] + db [ k, l] ); ga [ i, j ] += tmp; gb [ i, j ] += tmp
    tmp = 4 * x * ( da [ i, j ] + db [ i, j] ); ga [ k, l ] += tmp; gb [ k, l ] += tmp
    ga [ i, k ] -= 2 * x * da [ j, l ]; gb [ i, k ] -= 2 * x * db [ j, l ]
    ga [ i, l ] -= 2 * x * da [ j, k ]; gb [ i, l ] -= 2 * x * db [ j, k ]
    ga [ j, k ] -= 2 * x * da [ i, l ]; gb [ j, k ] -= 2 * x * db [ i, l ]
    ga [ j, l ] -= 2 * x * da [ i, k ]; gb [ j, l ] -= 2 * x * db [ i, k ]
  end

  for i in 0..nbf-1 do
    for j in 0..i-1 do
      tmpa = ga [ j, i ] + ga [ i, j ]; tmpb = gb [ j, i ] + gb [ i, j ] 
      ga [ j, i ] = ga [ i, j ] = tmpa / 2.0; gb [ j, i ] = gb [ i, j ] = tmpb / 2.0
    end
  end
  return ga, gb
end

def uhf( molinfo, s, h, nz_g )

    nbf = molinfo.get_nbf
    nel = molinfo.get_nel
    nelb = (nel - molinfo.get_spin) / 2
    nela = nel - nelb
  
#
#   Canonical orbital transformation matrix : tm
#
    eval, evec = GSL::Eigen::symmv( s )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_DESC )
#   EmolUtil::print_2dim_ary_tri(s.to_a, "Overlap matrix", 10, "%12.6f")
#   EmolUtil::print_2dim_ary_tri(h.to_a, "Hcore matrix", 10, "%12.6f")

    num_thrown = 0; thrown = EmolConsts::THR_THROWN
    eval.each do | x |
      if x < thrown
        num_thrown += 1
      end
    end
    eval = eval.subvector( eval.size - num_thrown )
    evec = evec.submatrix( nil, 0, evec.size[0] - num_thrown )
    for i in 0..eval.size-1
      eval[i] = 1.0 / Math.sqrt( eval[i] )
    end
#   EmolUtil::print_2dim_ary_with_value(evec.to_a, eval.to_a, "Rectangular martrix", 10, "%12.6f")
    tm = evec * GSL::Matrix.diagonal( eval )
#                                                                      X = U s^{-1/2}   Eq(3.172)
#   EmolUtil::print_2dim_ary((tm.transpose * s * tm).to_a, "Check X(T)SX = 1", 10, "%12.6f")

#
#   make initial guess for mo from core Hamiltonian h
#
    eval, evec = GSL::Eigen::symmv( tm.transpose * h * tm )
    GSL::Eigen::symmv_sort( eval, evec, type=GSL::Eigen::SORT_VAL_ASC )
    caoa = caob = tm * evec
    EmolUtil::print_2dim_ary_with_value( caoa.to_a, eval.to_a, "Eigenvalues and Eigenvectors of Hcore", 8, "%12.6f" )
    daoa = Density.new( caoa, nela, 1.0 ); daob = Density.new( caob, nelb, 1.0 )
#
#   scf iteration
#
    diisa = DIIS.new; diisb = DIIS.new; 
    errora = GSL::Vector.alloc( [100.0, 100.0] ); errorb = GSL::Vector.alloc( [100.0, 100.0] )
    iter = -1; ga = nil; gb = nil;
    printf("\n  Iter:     Energy                     Error \n")
#   while (errora * errora.col + errorb * errorb.col)/2.0 > EmolConsts::THR_CONV
    while (errora * errora.col + errorb * errorb.col)/2.0 > molinfo.get_thresh_hf
      ga, gb = calc_gmat_nz_eris_uhf( nz_g, daoa.get_body, daob.get_body, nbf )
      fa = h + ga; fb = h + gb; 
      iter += 1
#
      total_energy = (daoa.get_body * (h + fa) + daob.get_body * (h + fb)).trace / 2.0 + molinfo.get_e_nucl
      if iter > 0 then
        printf( "   %3d     %20.16e     %20.16e \n", iter, total_energy, errora*errora.col+errorb*errorb.col )
      else
        printf( "   %3d     %20.16e  \n",            iter, total_energy )
      end

      fmoa = caoa.transpose * fa * caoa; fmob = caob.transpose * fb * caob
      erra = []
      for i in 0..nela-1
        for j in nela..nbf-1
          erra.push( fmoa[i, j] )
        end
      end
      errb = []
      for i in 0..nelb-1
        for j in nelb..nbf-1
          errb.push( fmob[i, j] )
        end
      end
      errora = GSL::Vector.alloc( erra ); errorb = GSL::Vector.alloc( errb )
      diisa.err_push( errora ); diisb.err_push( errorb )
      diisa.fock_push( fa ); diisb.fock_push( fb )
#
#     printf("dmat, gmat, fmat: \n")
#     ij = -1
#     for i in 0..g.size1-1
#       for j in 0..i
#         ij += 1
#         printf(" %5d  %5d  %6d  %20.15e  %20.15e  %20.15e  %20.15e \n", j, i, ij, d[j,i], h[j,i], g[j,i], f[j,i])
#       end
#     end
#
      evala, evec = GSL::Eigen::symmv( tm.transpose * diisa.get_averageFock * tm )
      GSL::Eigen::symmv_sort( evala, evec, type=GSL::Eigen::SORT_VAL_ASC )
      caoa = tm * evec
      daoa = Density.new( caoa, nela, 1.0 )
      evalb, evec = GSL::Eigen::symmv( tm.transpose * diisb.get_averageFock * tm )
      GSL::Eigen::symmv_sort( evalb, evec, type=GSL::Eigen::SORT_VAL_ASC )
      caob = tm * evec
      daob = Density.new( caob, nelb, 1.0 )

  end
  EmolUtil::print_2dim_ary_with_value(caoa.to_a, evala.to_a, "Eigenvectors of alpha spin", 8, "%12.6f")
  EmolUtil::print_2dim_ary_with_value(caob.to_a, evalb.to_a, "Eigenvectors of beta spin", 8, "%12.6f")

# return total_energy, fa, fb, ga, gb, caoa, caob, daoa.get_body, daob.get_body, tm
end

molinfo = Molinfo.new( "input" )
s, h, nz_g = emol_int( molinfo )
#uhf_energy, f, g, cao, dao, tm = uhf( molinfo, s, h, nz_g )
uhf( molinfo, s, h, nz_g )
