FC = ifort 
FILE=ic_gen
FILEF=$(FILE).f
OUTPUT=$(FILE).x
FFTW_INC=/u/shared/programs/x86_64/fftw/3.3.4/include
FFTW_LIB=/u/shared/programs/x86_64/fftw/3.3.4/lib
OPT = -g
LFLAGS = -lm

all: 
	$(FC) $(OPT) -o $(OUTPUT) -I$(FFTW_INC) $(FFTW_LIB)/libfftw3.a $(FILEF) $(LFLAGS)
	echo "$(FILE)"

clean:
	rm *.mod $(FILE) 
