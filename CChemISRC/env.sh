#!/bin/bash
# SCRIPT=$(readlink -f $0) only for linux
# Works on both linux and Mac OS X.

function read_link() {
    local path=$1
        if [ -d $path ] ; then
        local abspath=$(cd $path; pwd)
        else
            local prefix=$(cd $(dirname -- $path) ; pwd)
            local suffix=$(basename $path)
            local abspath="$prefix/$suffix"
        fi
    echo $abspath
}
SCRIPT=$(read_link $0)
export  WORKDIR=`dirname $SCRIPT`
echo WORKDIR=$WORKDIR
export LIBCCHEMIDIR=$WORKDIR/cchemilib
export LIBCCHEMISRC=$LIBCCHEMIDIR/src
export LIBCCHEMI=$LIBCCHEMIDIR/lib/libcchemi.so
echo $LIBCCHEMISRC
