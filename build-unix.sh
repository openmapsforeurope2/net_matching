#!/bin/sh
set -e

NB_PROC=`grep processor /proc/cpuinfo | wc -l`
NB_PROC=$(( $NB_PROC - 2))
if [ $NB_PROC -lt 1 ]
then
    NB_PROC=1
fi

cmake -S . -B $CONFIG \
      -G "Unix Makefiles" \
      -DCMAKE_BUILD_TYPE=$CONFIG \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
      -DCMAKE_MODULE_PATH=$CMAKE_MODULE_PATH \
      -DCMAKE_INSTALL_PREFIX=$APP_DIR_SRC
cmake --build $CONFIG -j$NB_PROC
cmake --install $CONFIG
rm -rf $CONFIG

