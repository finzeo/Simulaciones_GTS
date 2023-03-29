#! /usr/bin/env bash

USFDIR=$1
## retv=`grep '^set.CXXOPTFLAGS' $USFDIR/cmake/csuserfuns.cmake`
grep '^set.CXXOPTFLAGS' $USFDIR/cmake/csuserfuns.cmake | perl -pe 's/.*\"(.*)\".*/$1/;'
