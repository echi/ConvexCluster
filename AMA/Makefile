FC = ifort
FCFLAGS = -arch x86_64 -fno-common -O3 -parallel -c
FLFLAGS = -arch x86_64 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation

MKL_INCLUDE = /opt/intel/mkl/include
MKL_PATH = /opt/intel/mkl/lib
LIBS= \
  $(MKL_PATH)/libmkl_intel_lp64.a \
  $(MKL_PATH)/libmkl_intel_thread.a \
  $(MKL_PATH)/libmkl_core.a  \
  $(MKL_PATH)/libmkl_blas95_lp64.a \
  $(MKL_PATH)/libmkl_lapack95_lp64.a \
  /opt/intel/lib/libiomp5.a \
  -lpthread -lm

all: prox.so test-prox.so ama.so

ama.so: ama.o prox.o hilbert.o sort.o
	$(FC) $(FLFLAGS) $(LIBS) -o ama.so ama.o prox.o hilbert.o sort.o

ama.o: ama.f90 prox.mod constants.mod hilbert.mod
	$(FC) $(FCFLAGS) -I$(MKL_INCLUDE) ama.f90 -o ama.o

test-prox.so: test-prox.o prox.o hilbert.o sort.o
	$(FC) $(FLFLAGS) -o test-prox.so test-prox.o prox.o hilbert.o sort.o

test-prox.o: test-prox.f90 prox.mod sort.mod
	$(FC) $(FCFLAGS) -o test-prox.o test-prox.f90

prox.so: prox.o hilbert.o constants.o sort.o
	$(FC) $(FLFLAGS) -o prox.so prox.o hilbert.o constants.o sort.o

prox.mod: prox.f90 hilbert.mod constants.mod sort.mod
	$(FC) $(FCFLAGS) -o prox.o prox.f90

prox.o: prox.f90 hilbert.mod constants.mod sort.mod
	$(FC) $(FCFLAGS) -o prox.o prox.f90

sort.o: sort.f90 constants.mod
	$(FC) $(FCFLAGS) -o sort.o sort.f90

sort.mod: sort.f90 constants.mod
	$(FC) $(FCFLAGS) -c sort.f90

hilbert.mod: hilbert.f90 constants.mod
	$(FC) $(FCFLAGS) -I$(MKL_INCLUDE) -c hilbert.f90

constants.mod: constants.f90
	$(FC) $(FCFLAGS) -c constants.f90

clean:
	rm -rf *~ *.o *.so *.mod