#!/bin/sh

DIR=$PWD
RUBYLIB=$DIR/../src:$DIR/../modules:$RUBYLIB
EMOL_HOME="$DIR/.."
export RUBYLIB EMOL_HOME

echo 'ruby path:' `which ruby`

time ruby mod_input.rb h2o-dz.inp
time ruby mod_int.rb
time ruby mod_rhf.rb

time ruby mod_superci.rb h2o-dz.inp

