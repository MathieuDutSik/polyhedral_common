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

# downloading netcdf functionality and compiling it
RUN apt-get install -y libcurl4-openssl-dev m4
RUN wget https://www.gfd-dennou.org/arch/netcdf/unidata-mirror/netcdf-c-4.8.0.tar.gz
RUN wget https://fossies.org/linux/misc/netcdf-cxx4-4.3.1.tar.gz
RUN tar -vxzf netcdf-c-4.8.0.tar.gz && cd netcdf-c-4.8.0 && export LDFLAGS="-L/opt/hdf5-1.12.0/lib -lhdf5_hl -lhdf5 -Wl,-rpath,/opt/hd# f5-1.12.0/lib" && export CPPFLAGS="-I/opt/hdf5-1.12.0/include" && ./configure --prefix=/opt/netcdf-4.8.0 && make && make install
RUN tar -vxzf netcdf-cxx4-4.3.1.tar.gz && cd netcdf-cxx4-4.3.1 && export LDFLAGS="-L/opt/netcdf-4.8.0/lib -lnetcdf -Wl,-rpath,/opt/netcdf-4.8.0/lib" && export CPPFLAGS="-I/opt/netcdf-4.8.0/include -I/opt/hdf5-1.12.0/include" && ./configure --prefix=/opt/netcdf-cxx4-4.3.1 && make && make install
ENV NETCDF_CXX_PATH /opt/netcdf-cxx4-4.3.1
ENV NETCDF_C_PATH /opt/netcdf-4.8.0
ENV NETCDF_CXX_ALLINC "-I$NETCDF_CXX_PATH/include -I$NETCDF_C_PATH/include"
ENV NETCDF_CXX_LINK "-L$NETCDF_CXX_PATH/lib -lnetcdf_c++4 -Wl,-rpath,$NETCDF_CXX_PATH/lib"


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


# mpfrc++
ENV MPFR_INCDIR /usr/include
ENV MPFR_LINK "-lmpfr"
ENV MPREAL_INCLUDE "-I/opt/mpfrc++-3.6.8"
RUN cp -r GIT/polyhedral_common/basic_common_cpp/ExternalLib/mpfrc++-3.6.8 /opt/mpfrc++-3.6.8


# openmpi
RUN apt-get install -y openmpi-bin libopenmpi-dev
ENV MPI_INCLUDE_CPP "-I/usr/lib/x86_64-linux-gnu/openmpi/include/"
ENV MPI_LINK_CPP "-lmpi -lboost_mpi -lboost_serialization -lmpi_cxx"


# Compiling polyhedral functionality
RUN cd GIT/polyhedral_common && ./compile.sh
#RUN cd GIT/polyhedral_common && ./compile_dualdesc.sh
