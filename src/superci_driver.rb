# coding: utf-8
require 'emol_nl'
require 'emol_superci'

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
