#!/bin/sh
DOCKER_NAME=net_matching_dev
DOCKER_TAG="latest"

docker build --no-cache -t $DOCKER_NAME:$DOCKER_TAG -f Dockerfile.dev ./..