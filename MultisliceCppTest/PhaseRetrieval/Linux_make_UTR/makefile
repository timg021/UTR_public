IDIR =../include
CC=gcc-9
CPP=g++-9
CFLAGS= -O3 #-g3
CFLAGS1= -Wno-unused-result -Wno-format-overflow -std=gnu++17 -fopenmp
CFLAGSC= -Wno-unused-result -Wno-format-overflow -fopenmp

#User will need to modify these for their own system
#FFTWdir, must contain libfftw3
FFTWdir = /usr/lib/x86_64-linux-gnu

ODIR = obj
LDIR = $(FFTWdir)
XArraydir = ../XArrayLibrary
UTRdir = ../PhaseRetrieval
_DEPS = XA_filt2.h XArray1D.h XA_spln3.h fftwd3drc.h IXAHWave.h XA_data.h XA_head1.h XA_ini.h XArray2D.h XA_tie.h OouraFft.h stdafx.cpp XA_fft2.h XA_head2.h XA_iwfr.h XArray3D.h XA_tiff.h fftwd3frc.h stdafx.h XA_fftr2.h XA_head3.h XA_move2.h XA_move3.h XArray.h IXAHead.h XA_file.h XAHWave.h XA_spln2.h UTR.h
DEPS = $(patsubst %,$(XArraydir)/%,$(_DEPS))

_OBJ =  Ooura_fftsg_flt.o Ooura_fftsg2d_dbl.o Ooura_fftsg_dbl.o XArray.o fftwd2c.o fftwd3drc.o fftwd3frc.o IXAHWave.o XA_ini.o XA_file.o XAHWave.o Ooura_fftsg2d_flt.o FindPeaks.o UTR.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS= -lm -lpthread -lfftw3 -lfftw3_threads -lfftw3f_threads -lfftw3f -lfftw3_omp -lfftw3f_omp

# $(ODIR)/%.o: %.c $(XarrayDEPS)
# 	echo $(ODIR)
# 	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

#UTR

$(ODIR)/Ooura_fftsg2d_flt.o: $(XArraydir)/Ooura_fftsg2d_flt.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS) 
	
$(ODIR)/Ooura_fftsg_flt.o: $(XArraydir)/Ooura_fftsg_flt.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg2d_dbl.o: $(XArraydir)/Ooura_fftsg2d_dbl.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg_dbl.o: $(XArraydir)/Ooura_fftsg_dbl.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/fftwd2c.o: $(XArraydir)/fftwd2c.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/fftwd3drc.o: $(XArraydir)/fftwd3drc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/fftwd3frc.o: $(XArraydir)/fftwd3frc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/IXAHWave.o: $(XArraydir)/IXAHWave.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_ini.o: $(XArraydir)/XA_ini.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XArray.o: $(XArraydir)/XArray.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)
	
$(ODIR)/XA_file.o: $(XArraydir)/XA_file.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XAHWave.o: $(XArraydir)/XAHWave.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/FindPeaks.o: $(UTRdir)/FindPeaks.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary $(CFLAGS)

$(ODIR)/UTR.o: $(UTRdir)/UTR.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary $(CFLAGS)

UTR: $(OBJ)
	$(CPP) $(CFLAGS1) $(OBJ) $(CFLAGS) $(LIBS) -o UTR

#clean:
	rm $(ODIR)/*.o


