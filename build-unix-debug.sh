#!/bin/sh
NB_PROC=`grep processor /proc/cpuinfo | wc -l`
NB_PROC=$(( $NB_PROC - 2))
if [ $NB_PROC -lt 1 ]
then
    NB_PROC=1
fi

rm -rf build-unix-debug bin

mkdir build-unix-debug
cd build-unix-debug 
cmake -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_INSTALL_PREFIX=$APP_DIR_SRC .. && make -j$NB_PROC install
cd ..
rm -r build-unix-debug
