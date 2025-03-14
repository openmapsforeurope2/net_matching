FROM libepg

ARG GIT_BRANCH
ARG NB_PROC

ENV APP_DIR_SRC=/usr/local/src/net_matching/                                              

WORKDIR $APP_DIR_SRC/../

COPY . $APP_DIR_SRC

WORKDIR $APP_DIR_SRC

RUN ls

RUN git checkout $GIT_BRANCH

RUN mkdir build-unix-release
RUN cd build-unix-release && cmake .. -DCMAKE_INSTALL_PREFIX=$APP_DIR_SRC -DCMAKE_BUILD_TYPE=Release && make -j$NB_PROC install
RUN rm -r build-unix-release

ENV PATH=$PATH:$APP_DIR_SRC/bin
