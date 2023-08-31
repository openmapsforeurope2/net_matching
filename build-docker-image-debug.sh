#!/bin/sh
PROJECT_NAME=tn_matching
DOCKER_TAG="latest"

docker build --no-cache -t $PROJECT_NAME:$DOCKER_TAG -f Dockerfile.debug .