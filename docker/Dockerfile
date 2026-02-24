FROM libepg AS build

ARG GIT_BRANCH
ARG NB_PROC

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

RUN git checkout $GIT_BRANCH \
 && mkdir build-unix-release \
 && cd build-unix-release && cmake .. -DCMAKE_INSTALL_PREFIX=$APP_DIR_SRC -DCMAKE_BUILD_TYPE=Release && make -j$NB_PROC install \
 && cd .. && rm -r build-unix-release

RUN mkdir -p /deps \
 && lddtree -l /usr/local/src/net_matching/bin/net_matching \
 | grep '^/' \
 | xargs -I{} cp -v {} /deps

# ──────────────────────────────
# Image finale minimale
# ──────────────────────────────
FROM debian:12

RUN apt-get update -y \
 && apt-get -y install libssl3 rclone

# Dependencies
RUN mkdir -p /usr/local/src/deps

COPY --from=build /deps /usr/local/src/deps

ENV LD_LIBRARY_PATH=/usr/local/src/deps

# App
ENV APP_DIR_SRC=/usr/local/src/net_matching/
ENV PATH=$PATH:$APP_DIR_SRC/bin

COPY --from=build $APP_DIR_SRC/bin $APP_DIR_SRC/bin
COPY --from=build $APP_DIR_SRC/config $APP_DIR_SRC/config

WORKDIR $APP_DIR_SRC