#!/bin/sh

DIR=$PWD
RUBYLIB=$DIR/../src:$DIR/../modules:$RUBYLIB
EMOL_HOME="$DIR/.."
export RUBYLIB EMOL_HOME

echo 'ruby path:' `which ruby`

time (
        ruby mod_input.rb h2o-dz.inp
        ruby mod_int.rb
        ruby mod_rhf.rb
    )
time ruby mod_trn.rb

time ruby mod_symbolic_expr.rb
time ruby mod_hmatrix_symbolic.rb
time ruby util_symbolic2real_expr.rb
time ruby mod_no.rb

#ruby mod_expr.rb
#ruby mod_hmatrix.rb
#ruby mod_no.rb
