#!/bin/sh

VERSION=1.19.1

srcdir=${0%/version.sh}
if [ -e $srcdir/.git ]
then

    v=`cd $srcdir && git describe --always --match '[0-9].[0-9]*' --dirty`
    case $v in
        [0-9]*.[0-9]*) VERSION="$v" ;;
        [0-9a-f][0-9a-f]*) VERSION="$VERSION-1-g$v" ;;
    esac
fi

if [ "$1" = "numeric" ]
then
    v1=`expr "$VERSION" : '\([0-9]*\)'`
    v2=`expr "$VERSION" : '[0-9]*.\([0-9]*\)'`
    v3=`expr "$VERSION" : '[0-9]*.[0-9]*.\([0-9]*\)'`
    if [ -z "`expr "$VERSION" : '\([0-9.]*\)$'`" ]
    then
        VERSION="$v1.$v2.255"
    else
        VERSION="$v1.$v2${v3:+.}$v3"
    fi
fi

echo $VERSION
