# The Docker environment can be compiled by doing
# 
# docker build -t standard_working_environment:build-env .
#
# Ask Mathieu Dutour Sikiric (mathieu.dutour@gmail.com) for questions.

FROM ubuntu

# Dockerfile maintainer (original author: Mathieu Dutour Sikiric)
MAINTAINER Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>


# Install Ubuntu packages
RUN apt-get update -y
RUN apt-get install -y libgdbm-dev libsnappy-dev make pkg-config rlwrap screen software-properties-common sudo unzip zlib1g-dev zsh
RUN apt-get install -y wget emacs joe
RUN apt-get install -y g++ gcc gfortran
RUN apt-get install -y git mercurial cmake

# GMP
RUN apt-get install -y libgmp-dev
ENV GMP_INCDIR "/usr/include"
ENV GMP_CXX_LINK "-lm -lgmp -lgmpxx"

# BOOST
RUN apt-get install -y libboost-dev
ENV BOOST_INCDIR "/usr/include"

# downloading eigen (no compiling needed)
RUN mkdir -p /opt
RUN cd /opt && hg clone https://bitbucket.org/eigen/eigen/
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
RUN cd GIT && git clone https://github.com/MathieuDutSik/polyhedral_common.git --recursive


# libshort
ENV LIBSHORT_PATH /opt/libshort
ENV LIBSHORT_INCDIR "$LIBSHORT_PATH/include"
ENV LIBSHORT_LINK "-L$LIBSHORT_PATH/lib -lshort -Wl,-rpath,$LIBSHORT_PATH/lib"
RUN cd GIT/polyhedral_common/ExternalLib/LibShort && make -f Makefile_lib && make -f Makefile_lib install


# bliss
ENV LIBBLISS_PATH /opt/libbliss
ENV LIBBLISS_INCDIR "$LIBBLISS_PATH/include"
ENV LIBBLISS_LINK "-L$LIBBLISS_PATH/lib -lbliss -Wl,-rpath,$LIBBLISS_PATH/lib"
RUN cd GIT/polyhedral_common/ExternalLib/bliss && make -f Makefile_lib && make -f Makefile_lib install


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

# Compiling polyhedral functionality
RUN cd GIT/polyhedral_common/src_latt && make
RUN cd GIT/polyhedral_common/src_poly && make
RUN cd GIT/polyhedral_common/src_copos && make
RUN cd GIT/polyhedral_common/src_short && make
RUN cd GIT/polyhedral_common/src_perfect && make


