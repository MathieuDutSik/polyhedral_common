# This dockerfile allows you to access to functionality of C++ polyhedral
#
# It is used in the following way:
# ---Copy this "Dockerfile" to a directory, e.g. "docker".
#    Do NOT rename the file.
# ---Compile the dockerfile with the command
#    docker build -t standard_working_environment:build-env .
# ---Run the docker and access a bash by doing
#    docker run -i -t standard_working_environment:build-env /bin/bash
#
# Ask Mathieu Dutour Sikiric (mathieu.dutour@gmail.com) for questions.

FROM ubuntu

# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
MAINTAINER Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>


# For unclear reason, all ubuntu installs must start by the following
RUN apt-get update -y


# This specific package requires to be made non-interactive.
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata


# Install classic Ubuntu packages
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe
RUN apt-get install -y g++ gcc gfortran
RUN apt-get install -y git cmake


# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lm -lgmp -lgmpxx"


# BOOST
RUN apt-get install -y libboost-dev
ENV BOOST_INCDIR "/usr/include"


# downloading eigen (no compiling needed)
RUN mkdir -p /opt
RUN cd /opt && git clone https://gitlab.com/libeigen/eigen.git
ENV EIGEN_PATH /opt/eigen


# downloading permlib
RUN cd /opt && git clone https://github.com/MathieuDutSik/permlib_mod
ENV PERMLIB_INCDIR /opt/permlib_mod/include


# Threading Building Blocks
RUN apt-get install -y libtbb-dev
ENV TBB_INCDIR=/usr/include
ENV TBB_LINK="-ltbb"


# Gnu Scientific library
RUN apt-get install -y libgsl-dev
ENV GSL_INCDIR=/usr/include
ENV GSL_C_LINK="-lgsl -lgslcblas"


# cloning polyhedral code
RUN mkdir -p GIT
RUN cd GIT &&  git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive


# libshort
ENV LIBSHORT_PATH /opt/libshort
ENV LIBSHORT_INCDIR "$LIBSHORT_PATH/include"
ENV LIBSHORT_LINK "-L$LIBSHORT_PATH/lib -lshort -Wl,-rpath,$LIBSHORT_PATH/lib"
RUN cd GIT/polyhedral_common/ExternalLib/LibShort && make -f Makefile_lib && make -f Makefile_lib install


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
ENV NAUTY_LINK "-L$NAUTY_PATH/lib -lnauty -Wl,-rpath,$NAUTY_PATH/lib"
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty && mkdir build
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build && cmake -DCMAKE_INSTALL_PREFIX:PATH=$NAUTY_PATH ..
RUN cd GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build && make all install


# GLPK
RUN wget https://ftp.gnu.org/gnu/glpk/glpk-4.65.tar.gz
RUN tar -vxzf glpk-4.65.tar.gz && cd glpk-4.65 && ./configure --prefix=/opt/glpk-4.65 && make && make install
ENV GLPK_PATH /opt/glpk-4.65
ENV GLPK_INCLUDE "-I$GLPK_PATH/include"
ENV GLPK_LINK="-L$GLPK_PATH/lib -lglpk -Wl,-rpath,$GLPK_PATH/lib"
RUN ls -la /opt/glpk-4.65/lib


# mpfrc++
RUN apt-get install -y libmpfr-dev
ENV MPFR_INCDIR /usr/include
ENV MPFR_LINK "-lmpfr"
ENV MPREAL_INCLUDE "-I/opt/mpfrc++-3.6.2"
RUN cp -r GIT/polyhedral_common/ExternalLib/mpfrc++-3.6.2 /opt/mpfrc++-3.6.2


# flint
RUN apt-get install -y libflint-dev
ENV FLINT_INCDIR /include
ENV FLINT_LINK -lflint


# CDDLIB
RUN apt-get install -y autoconf autotools-dev libtool
RUN apt-get install -y texlive-latex-base # needed for the manual to be compiled
RUN git clone https://github.com/cddlib/cddlib
RUN cd cddlib && ./bootstrap && ./configure --prefix=/opt/cddlib && make && make install
ENV CDDLIB_PATH /opt/cddlib
ENV CDDLIB_INCLUDE "-I$CDDLIB_PATH/include/cddlib"
ENV CDDLIB_LINK "-L$CDDLIB_PATH/lib  -lcdd -lcddgmp -Wl,-rpath,$CDDLIB_PATH/lib"


# LINBOX
RUN apt-get install -y wget libflint-dev libntl-dev libmpfr-dev
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


# Compiling polyhedral functionality
RUN cd GIT/polyhedral_common/src_linbox && make
RUN cd GIT/polyhedral_common/src_latt && make
RUN cd GIT/polyhedral_common/src_poly && make
RUN cd GIT/polyhedral_common/src_copos && make
RUN cd GIT/polyhedral_common/src_short && make
RUN cd GIT/polyhedral_common/src_perfect && make
# RUN cd GIT/polyhedral_common/src_ctype_mpi && make


