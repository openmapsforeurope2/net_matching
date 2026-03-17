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

DOCKER_BUILDKIT=1 docker build --no-cache -t $DOCKER_NAME:$DOCKER_TAG -f Dockerfile ./..
