#!/bin/sh

DIR=$PWD
RUBYLIB=$DIR/../src:$DIR/../modules:$RUBYLIB
EMOL_HOME="$DIR/.."
export RUBYLIB EMOL_HOME

echo 'ruby path:' `which ruby`

ruby mod_superci_restart.rb
