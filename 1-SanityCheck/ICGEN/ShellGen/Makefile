FC = ifort 
OPT = -O3 -convert big_endian 
NAME=ic_gen
SRC=$(NAME).f
EXEC=$(NAME).x

all:
	$(FC) $(SRC) $(OPT) -I$(TWO_DECOMP_FFT_INC) -I$(FFTW_INC) -o $(EXEC) -L$(TWO_DECOMP_FFT_LIB) -L$(FFTW_LIB) -l2decomp_fft -lfftw3 -lm

clean:
	rm *.mod *.x
