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

FROM ubuntu:22.04

# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
LABEL maintainer="mathieu.dutour@gmail.com"


# For unclear reason, all ubuntu installs must start by the following
RUN apt-get update -y


# This specific package requires to be made non-interactive.
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata


# Install classic Ubuntu packages
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe
RUN apt-get install -y g++-11 g++ gcc gfortran clang
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


# Threading Building Blocks
RUN apt-get install -y libtbb-dev
ENV TBB_INCDIR=/usr/include
ENV TBB_LINK="-ltbb"


# Gnu Scientific library
RUN apt-get install -y libgsl-dev
ENV GSL_INCDIR=/usr/include
ENV GSL_C_LINK="-lgsl -lgslcblas"


# GLPK
RUN apt-get install -y libglpk-dev
ENV GLPK_PATH /usr
ENV GLPK_INCLUDE "-I$GLPK_PATH/include"
ENV GLPK_LINK="-L$GLPK_PATH/lib -lglpk -Wl,-rpath,$GLPK_PATH/lib"


# flint
RUN apt-get install -y libflint-dev libntl-dev libmpfr-dev
ENV FLINT_INCDIR /include
ENV FLINT_LINK -lflint


# CDDLIB
RUN apt-get install -y autoconf autotools-dev libtool
RUN apt-get install -y texlive-latex-base # needed for the manual to be compiled
RUN git clone https://github.com/cddlib/cddlib
RUN cd cddlib && ./bootstrap && ./configure --prefix=/opt/cddlib && make && make install
ENV CDDLIB_PATH /opt/cddlib
ENV CDDLIB_INCLUDE "-I$CDDLIB_PATH/include/cddlib"
ENV CDDLIB_DOUBLE_LINK "-L$CDDLIB_PATH/lib  -lcdd -Wl,-rpath,$CDDLIB_PATH/lib"
ENV CDDLIB_GMP_LINK "-L$CDDLIB_PATH/lib  -lcddgmp -Wl,-rpath,$CDDLIB_PATH/lib"


# PPL
RUN wget https://www.bugseng.com/products/ppl/download/ftp/releases/1.2/ppl-1.2.tar.gz
RUN tar -vxzf ppl-1.2.tar.gz && cd ppl-1.2 && ./configure && make

# Downloading PPL and compiling it
RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.gz
RUN tar -vxzf hdf5-1.12.0.tar.gz && cd hdf5-1.12.0 && ./configure --with-default-api-version=v18 --prefix=/opt/hdf5-1.12.0 && make && make install
RUN mkdir -p /bin && cd /bin && ln -s ../ppl-1.2/demos/ppl_lcdd/.libs/ppl_lcdd

# downloading netcdf functionality and compiling it
RUN apt-get install -y libcurl4-openssl-dev m4
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-c-4.8.0.tar.gz
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-fortran-4.5.3.tar.gz
RUN wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-cxx4-4.3.1.tar.gz
RUN tar -vxzf netcdf-c-4.8.0.tar.gz && cd netcdf-c-4.8.0 && export LDFLAGS="-L/opt/hdf5-1.12.0/lib -lhdf5_hl -lhdf5 -Wl,-rpath,/opt/hdf5-1.12.0/lib" && export CPPFLAGS="-I/opt/hdf5-1.12.0/include" && ./configure --prefix=/opt/netcdf-4.8.0 && make && make install
RUN tar -vxzf netcdf-fortran-4.5.3.tar.gz && cd netcdf-fortran-4.5.3 && export LDFLAGS="-L/opt/netcdf-4.8.0/lib -lnetcdf -Wl,-rpath,/opt/netcdf-4.8.0/lib" && export CPPFLAGS="-I/opt/netcdf-4.8.0/include -I/opt/hdf5-1.12.0/include" && ./configure --prefix=/opt/netcdf-fortran-4.5.3 && make && make install
RUN tar -vxzf netcdf-cxx4-4.3.1.tar.gz && cd netcdf-cxx4-4.3.1 && export LDFLAGS="-L/opt/netcdf-4.8.0/lib -lnetcdf -Wl,-rpath,/opt/netcdf-4.8.0/lib" && export CPPFLAGS="-I/opt/netcdf-4.8.0/include -I/opt/hdf5-1.12.0/include" && ./configure --prefix=/opt/netcdf-cxx4-4.3.1 && make && make install
ENV NETCDF_CXX_PATH /opt/netcdf-cxx4-4.3.1
ENV NETCDF_C_PATH /opt/netcdf-4.8.0
ENV NETCDF_CXX_ALLINC "-I$NETCDF_CXX_PATH/include -I$NETCDF_C_PATH/include"
ENV NETCDF_CXX_LINK "-L$NETCDF_CXX_PATH/lib -lnetcdf_c++4 -Wl,-rpath,$NETCDF_CXX_PATH/lib"




# LINBOX
RUN wget https://raw.githubusercontent.com/linbox-team/linbox/master/linbox-auto-install.sh
RUN chmod +x linbox-auto-install.sh
RUN ./linbox-auto-install.sh --enable-openblas=yes
RUN mkdir /opt/linbox
RUN mv /tmp/lib /opt/linbox
RUN mv /tmp/include /opt/linbox
RUN mv /tmp/bin /opt/linbox
ENV LINBOX_PATH /opt/linbox
ENV LINBOX_INCLUDE "-DHAVE_CONFIG_H -DDISABLE_COMMENTATOR -DNDEBUG -UDEBUG -fopenmp -I$LINBOX_PATH/include"
ENV LINBOX_LINK "-fopenmp -lntl -lmpfr -lflint -L$LINBOX_PATH/lib -lopenblas -lpthread -lgfortran -lgivaro -lgmpxx -lgmp -Wl,-rpath,$CDDLIB_PATH/lib"


# cloning polyhedral code
RUN mkdir -p GIT
RUN cd GIT && git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive


# libshort
ENV LIBSHORT_PATH /opt/libshort
ENV LIBSHORT_INCDIR "$LIBSHORT_PATH/include"
ENV LIBSHORT_LINK "-L$LIBSHORT_PATH/lib -lshort -Wl,-rpath,$LIBSHORT_PATH/lib"
RUN cd GIT/polyhedral_common/ExternalLib/LibShort  && make -f Makefile_lib && make -f Makefile_lib install


# bliss
ENV LIBBLISS_PATH /opt/bliss
ENV LIBBLISS_INCDIR $LIBBLISS_PATH/include
ENV LIBBLISS_LINK "-L$LIBBLISS_PATH/lib -lbliss -Wl,-rpath,$LIBBLISS_PATH/lib"
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/bliss && mkdir build
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/bliss/build && cmake -DCMAKE_INSTALL_PREFIX:PATH=$LIBBLISS_PATH ..
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/bliss/build && make all install


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
