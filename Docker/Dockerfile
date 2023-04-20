# The first stage will install build dependencies and compile the executable

FROM debian:bullseye AS builder

RUN set -ex;                                                                        \
    apt-get update;                                                                 \
    apt-get install -y g++ curl cmake libboost-all-dev libgsl-dev zlib1g-dev; \
    mkdir -p /usr/src;                                                              \
    cd /usr/src;                                                                    \
    curl -L https://github.com/hahnlab/CAGEE/releases/download/1.1/CAGEE1.1.tar.gz | tar -zxf -; \
    curl -L https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz | tar -zxf -;

RUN set -ex;              \
    cd /usr/src/eigen-3.4.0;  \
    mkdir build && cd build;  \
    cmake ..; make install

COPY . /usr/src/CAGEE

RUN set -ex;              \
    cd /usr/src/CAGEE;  \
    mkdir build && cd build;  \
    cmake ..; make; make install

# The second stage will install the runtime dependencies only and copy
# the compiled executables

FROM debian:bullseye AS runtime

RUN set -ex;         \
    apt-get update;  \
    apt-get install -y libboost-program-options-dev libomp-dev libgsl-dev

COPY --from=builder /usr/local/bin /usr/local/bin

ENTRYPOINT ["/usr/local/bin/cagee"]

