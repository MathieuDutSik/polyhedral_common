# This dockerfile allows you to access to functionality of the Indefinite forms
# domain for debugging.
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker".
#    Do NOT rename the file.
# ---Compile the dockerfile with the command
#    docker build -t indefinite:build-env .
# ---For building with no-cache option, use
#    docker build --no-cache -t indefinite:build-env .
# ---Run the docker and access a bash by doing
#    docker run -i -t indefinite:build-env /bin/bash
#
# Ask Mathieu Dutour Sikiric (mathieu.dutour@gmail.com) for questions.

FROM ubuntu:24.04

# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
LABEL maintainer="mathieu.dutour@gmail.com"

# For unclear reason, all ubuntu installs must start by the following
RUN apt-get update -y


# This specific package requires to be made non-interactive.
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata


# Install classic Ubuntu packages
RUN apt-get install -y libgdbm-dev libsnappy-dev pkg-config rlwrap
RUN apt-get install -y screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe valgrind
RUN apt-get install -y g++ gcc
RUN apt-get install -y git cmake make

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

# nauty/traces
ENV NAUTY_PATH /opt/nauty
ENV NAUTY_INCLUDE "-I$NAUTY_PATH/include"
ENV NAUTY_LINK "-L$NAUTY_PATH/lib -lnauty_shared -Wl,-rpath,$NAUTY_PATH/lib"
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty && mkdir build
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build && cmake -DCMAKE_INSTALL_PREFIX:PATH=$NAUTY_PATH ..
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build && make all install


# Compiling polyhedral functionality
RUN cd GIT/polyhedral_common/src_indefinite && make -f Makefile_indef_ci
