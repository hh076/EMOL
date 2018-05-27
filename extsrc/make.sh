#!/bin/sh

MODULES='molint netlib trn wig symbolic'
OPT=
#RUBY=$HOME/ruby-2.4.1/bin/ruby
RUBY=/usr/local/bin/ruby
MAKE=make

if [ x$1 = 'xclean' ] ; then
    OPT=clean
fi

for M in $MODULES
do
    (
        echo $M:
        cd $M
	touch *.c *.f
        $RUBY extconf.rb
        $MAKE $OPT
    )
done
