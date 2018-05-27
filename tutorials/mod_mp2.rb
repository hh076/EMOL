
require 'emol_molinfo'
require 'emol_util'
require 'emol_rhf'
require 'emol_trans'

load 'Data_molinfo'
load 'Data_rhf'
load 'Data_trnint'

molinfo  = Data_molinfo::molinfo_data
ene      = Data_rhf::mo_energies_data
ene_elec = Data_rhf::elec_energy_data
ene_hf   = Data_rhf::hf_energy_data
trnint   = Data_trnint::trnint_data
nao      = molinfo.nbf
nob      = molinfo.nob
nocc     = molinfo.nel / 2
ene_core = trnint.core_energy

eris_mo = trnint.eris_mo

def ij( i_, j_ )
    i = i_; j = j_
    if i_ < j_
        i = j_; j = i_
    end
    i * ( i + 1 ) / 2 + j
end

printf( "\n\nMP2 calculations:\n\n" )

emp2 = 0
for i in 0...nocc do
  for j in 0...nocc do
    for a in nocc...nob do
      ia     = ij( i, a )
      ja     = ij( j, a )
      for b in nocc...nob do
        ib   = ij( i, b )
        jb   = ij( j, b )
        iajb = ij( ia, jb )
        ibja = ij( ib, ja )
        eri1 = eris_mo[ iajb ]
        eri2 = eris_mo[ ibja ]
        den3 = ( ene[ i ] + ene[ j ] - ene[ a ] - ene[ b ] )
        emp2 += eri1 * ( 2 * eri1 - eri2 ) / den3
      end
    end
  end
end

printf( "E_HF:         %28.16e\n", ene_hf )
printf( "E_MP2:        %28.16e\n", emp2 )
printf( "E_HF + E_MP2: %28.16e\n", ene_hf + emp2 )

print "\n MP2 : completed.\n"
