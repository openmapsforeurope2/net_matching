#!/bin/sh
PROJECT_NAME=cp_connection
DOCKER_TAG="latest"

docker build --no-cache -t $PROJECT_NAME:$DOCKER_TAG -f Dockerfile.debug .