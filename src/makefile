# Builds JAGWAS in this folder:   make
#
# Needs a C++11 compiler, the Armadillo headers and BLAS/LAPACK (OpenBLAS by default). Armadillo
# calls BLAS/LAPACK directly (ARMA_DONT_USE_WRAPPER), so its runtime library libarmadillo is not
# needed. Set the variables below on the command line when the defaults do not fit, e.g.
#
#   make ARMA_INC=/opt/armadillo/include         Armadillo headers not on the default path
#   make LAPACK_LIBS="-llapack -lblas"           reference BLAS/LAPACK instead of OpenBLAS
#   make STATIC=1 BLAS_DIR=/path/to/OpenBLAS     static binary, like the one in the repository

# folder with the header file 'armadillo', if it is not on the default include path
ARMA_INC =
# folder with the BLAS/LAPACK library, if it is not on the default library path
BLAS_DIR =
LAPACK_LIBS = -lopenblas
# 1 = link statically (needs the static libraries, e.g. libopenblas.a)
STATIC = 0

CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall
# the bundled R pgamma code is included as system headers, so its warnings are not shown
CPPFLAGS = -DARMA_DONT_USE_WRAPPER -isystem ../thirdparty/R_pgamma $(if $(ARMA_INC),-I$(ARMA_INC))
LDFLAGS = $(if $(BLAS_DIR),-L$(BLAS_DIR))
LDLIBS = $(LAPACK_LIBS) -lpthread

ifeq ($(STATIC),1)
LDFLAGS += -static
# a static OpenBLAS also needs the Fortran runtime it was built with
LDLIBS += -lgfortran -lquadmath
endif

all: JAGWAS

JAGWAS: JAGWAS.o
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# JAGWAS.cpp includes the pgamma code, so it is recompiled when that changes
JAGWAS.o: JAGWAS.cpp $(wildcard ../thirdparty/R_pgamma/*.c ../thirdparty/R_pgamma/*.h)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: all clean
clean:
	rm -f JAGWAS.o JAGWAS
