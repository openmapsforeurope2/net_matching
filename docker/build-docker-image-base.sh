#!/bin/sh
PROJECT_NAME=net_matching_base
DOCKER_TAG="latest"

docker build --no-cache -t $PROJECT_NAME:$DOCKER_TAG -f Dockerfile.base ./..