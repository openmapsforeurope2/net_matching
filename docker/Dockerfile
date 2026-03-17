# =========================
# STAGE 1 : BUILD
# =========================
FROM libepg AS build

RUN apt-get update -y \
 && apt-get -y install git \
                       g++ \
                       gcc \ 
                       cmake \
                       make \
                       pax-utils \
 && rm -rf /var/lib/apt/lists/*

ENV APP_DIR_SRC=/usr/local/src/net_matching/                                              

WORKDIR $APP_DIR_SRC/../

COPY . $APP_DIR_SRC

WORKDIR $APP_DIR_SRC

RUN ./build-unix-release.sh

RUN mkdir -p /deps \
 && lddtree -l /usr/local/src/net_matching/bin/net_matching \
 | grep '^/' \
 | xargs -I{} cp -v {} /deps

# =========================
# STAGE 2
# =========================
FROM ubuntu:24.04

ENV APP_DIR_SRC=/usr/local/src/net_matching/

COPY --from=build $APP_DIR_SRC/bin    $APP_DIR_SRC/bin
COPY --from=build $APP_DIR_SRC/config $APP_DIR_SRC/config
COPY --from=build /deps               $APP_DIR_SRC/deps

ENV LD_LIBRARY_PATH=$APP_DIR_SRC/deps \
    PATH=$PATH:$APP_DIR_SRC/bin

WORKDIR $APP_DIR_SRC