#!/bin/sh

DIR=$PWD
RUBYLIB=$DIR/../src:$DIR/../modules:$RUBYLIB
EMOL_HOME="$DIR/.."
export RUBYLIB EMOL_HOME

echo 'ruby path:' `which ruby`

ruby mod_input.rb h2o-dz.inp
ruby mod_int.rb
ruby mod_rhf.rb

ruby mod_superci.rb h2o-dz.inp
