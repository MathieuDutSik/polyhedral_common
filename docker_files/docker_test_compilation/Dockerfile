# This dockerfile allows you to access to functionality of C++ polyhedral
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker".
#    Do NOT rename the file.
# ---Compile the dockerfile with the command
#    docker build -t test_compilation:build-env .
# ---For building with no-cache option, use
#    docker build --no-cache -t test_compilation:build-env .
# ---Run the docker and access a bash by doing
#    docker run -i -t test_compilation:build-env /bin/bash
#
# Ask Mathieu Dutour Sikiric (mathieu.dutour@gmail.com) for questions.

FROM ubuntu:26.04

# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
LABEL maintainer="mathieu.dutour@gmail.com"


# For unclear reason, all ubuntu installs must start by the following
RUN apt-get update -y


# This specific package requires to be made non-interactive.
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata


# Install classic Ubuntu packages
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe g++ gcc gfortran clang git cmake

# openmpi
RUN apt-get install -y openmpi-bin libopenmpi-dev


# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lgmp -lgmpxx"


# BOOST
RUN apt-get install -y libboost-dev libboost-mpi-dev libboost-serialization-dev
ENV BOOST_INCDIR "/usr/include"
ENV BOOST_LINK "-lboost_serialization"
ENV MPI_LINK_CPP "-lboost_mpi"


# downloading eigen (no compiling needed)
RUN mkdir -p /opt
RUN cd /opt && git clone https://gitlab.com/libeigen/eigen.git
ENV EIGEN_PATH /opt/eigen

# cloning polyhedral code
RUN mkdir -p GIT
RUN cd GIT && git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive


# bliss
RUN cd GIT && git clone https://github.com/MathieuDutSik/bliss
RUN cd GIT/bliss && mkdir build && cd build && cmake .. && make
ENV LIBBLISS_INCDIR /GIT/bliss/src
ENV LIBBLISS_LINK "-L/GIT/bliss/build -lbliss -Wl,-rpath,/GIT/bliss/build"


# nauty/traces
ENV NAUTY_PATH /opt/nauty
ENV NAUTY_INCLUDE "-I$NAUTY_PATH/include"
ENV NAUTY_LINK "$NAUTY_PATH/lib/libnauty.a"
RUN cd GIT && git clone --depth 1 --branch traces-grpsize-factors https://github.com/MathieuDutSik/nauty
RUN cd GIT/nauty && ./configure --prefix=$NAUTY_PATH && make && make install
