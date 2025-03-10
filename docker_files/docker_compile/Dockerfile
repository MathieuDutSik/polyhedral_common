# This dockerfile allows you to access to functionality of C++ polyhedral
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker".
#    Do NOT rename the file.
# ---Compile the dockerfile with the command
#    docker build -t polyhedralcpp:build-env .
# ---For building with no-cache option, use
#    docker build --no-cache -t polyhedralcpp:build-env .
# ---Run the docker and access a bash by doing
#    docker run -i -t polyhedralcpp:build-env /bin/bash
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
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe
RUN apt-get install -y g++ gcc gfortran clang
RUN apt-get install -y git cmake


# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lgmp -lgmpxx"


# BOOST
RUN apt-get install -y libboost-dev libboost-mpi-dev libboost-serialization-dev
ENV BOOST_INCDIR "/usr/include"
ENV BOOST_LINK "-lboost_serialization"


# downloading eigen (no compiling needed)
RUN mkdir -p /opt
RUN cd /opt && git clone https://gitlab.com/libeigen/eigen.git
ENV EIGEN_PATH /opt/eigen


# GLPK
RUN apt-get install -y libglpk-dev
ENV GLPK_PATH /usr
ENV GLPK_INCLUDE "-I$GLPK_PATH/include"
ENV GLPK_LINK="-L$GLPK_PATH/lib -lglpk -Wl,-rpath,$GLPK_PATH/lib"


# CDDLIB
RUN apt-get install -y autoconf autotools-dev libtool
RUN apt-get install -y texlive-latex-base # needed for the manual to be compiled
RUN git clone https://github.com/cddlib/cddlib
RUN cd cddlib && ./bootstrap && ./configure --prefix=/opt/cddlib && make && make install
ENV CDDLIB_PATH /opt/cddlib
ENV CDDLIB_INCLUDE "-I$CDDLIB_PATH/include/cddlib"
ENV CDDLIB_DOUBLE_LINK "-L$CDDLIB_PATH/lib  -lcdd -Wl,-rpath,$CDDLIB_PATH/lib"
ENV CDDLIB_GMP_LINK "-L$CDDLIB_PATH/lib  -lcddgmp -Wl,-rpath,$CDDLIB_PATH/lib"


# openmpi
RUN apt-get install -y openmpi-bin libopenmpi-dev
ENV MPI_INCLUDE_CPP "-I/usr/lib/x86_64-linux-gnu/openmpi/include/"
ENV MPI_LINK_CPP "-lmpi -lboost_mpi -lboost_serialization -lmpi_cxx"


# nauty/traces
ENV NAUTY_PATH /opt/nauty
ENV NAUTY_INCLUDE "-I$NAUTY_PATH/include"
ENV NAUTY_LINK "-L$NAUTY_PATH/lib -lnauty_shared -Wl,-rpath,$NAUTY_PATH/lib"
RUN mkdir -p GIT
RUN cd GIT && git clone https://github.com/MathieuDutSik/nauty
RUN cd GIT/nauty && ./configure --prefix=/opt/nauty
RUN cd GIT/nauty && make && make install


# cloning polyhedral code
RUN cd GIT && git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive


# bliss
ENV LIBBLISS_PATH /opt/bliss
ENV LIBBLISS_INCLUDE "-I$LIBBLISS_PATH/include"
ENV LIBBLISS_LINK "-L$LIBBLISS_PATH/lib -lbliss -Wl,-rpath,$LIBBLISS_PATH/lib"
RUN mkdir -p /opt/bliss/include
RUN mkdir -p /opt/bliss/lib
RUN cd /opt/bliss/lib && unzip /GIT/polyhedral_common/basic_common_cpp/ExternalLib/bliss-0.77.zip
RUN cd /opt/bliss/lib && cmake /opt/bliss/lib/bliss-0.77 && make
RUN cd /opt/bliss/include && cp /opt/bliss/lib/bliss-0.77/src/*.hh .


# Compiling polyhedral functionality
RUN cd GIT/polyhedral_common && ./compile_docker.sh
#RUN cd GIT/polyhedral_common && ./compile_dualdesc.sh
