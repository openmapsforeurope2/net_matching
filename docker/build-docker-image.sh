#!/bin/sh
DOCKER_NAME=net_matching

GIT_BRANCH=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
GIT_BRANCH_LOWER=$(echo $GIT_BRANCH | tr '[:upper:]' '[:lower:]')

DOCKER_TAG=$(head -n 1 ./../VERSION)

if [ $GIT_BRANCH_LOWER = "main" ]
then
    DOCKER_TAG="latest"
fi

echo $GIT_BRANCH
echo $DOCKER_TAG

NB_PROC=`grep processor /proc/cpuinfo | wc -l`
NB_PROC=$(( $NB_PROC - 2))
if [ $NB_PROC -lt 1 ]
then
    NB_PROC=1
fi
echo $NB_PROC

docker build --no-cache --build-arg GIT_BRANCH=$GIT_BRANCH --build-arg NB_PROC=$NB_PROC -t $DOCKER_NAME:$DOCKER_TAG -f Dockerfile ./..
