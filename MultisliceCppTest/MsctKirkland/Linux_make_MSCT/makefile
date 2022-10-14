CC=gcc-9
CPP=g++-9
CFLAGS= -O3 #-g3
CFLAGS1= -Wno-unused-result -Wno-format-overflow -std=gnu++17 -fopenmp
CFLAGSC= -Wno-unused-result -Wno-format-overflow

#User will need to modify these for their own system
#FFTWdir, must contain libfftw3
FFTWdir = /usr/lib/x86_64-linux-gnu

ODIR = obj
LDIR = -L $(FFTWdir)
XArraydir = ../XArrayLibrary
Temsimdir = ../XArrayLibrary
MSCTdir = ../MsctKirkland

_DEPS = XArray1D.h IXAHWave.h XA_data.h XA_head1.h XA_ini.h XArray2D.h OouraFft.h stdafx.h fftwd2fc.h XA_spln2.h XA_fft2.h XA_head2.h XArray3D.h XA_tiff.h stdafx.h XA_move2.h XA_nrrand.h XArray.h IXAHead.h XA_file.h XAHWave.h
DEPS = $(patsubst %,$(XArraydir)/%,$(_DEPS))

_DEQS = AddIce.h autoslic.hpp autosliccmd.h cfpix.hpp floatTIFF.hpp pdb.h slicelib.hpp
DEQS = $(patsubst %,$(Temsimdir)/%,$(_DEQS))

_OBJ =  Oourafft.o Ooura_fftsg_flt.o Ooura_fftsg2d_dbl.o Ooura_fftsg_dbl.o XArray.o IXAHWave.o XA_ini.o XA_file.o XAHWave.o fftwd2fc.o XA_spln2.o XA_fft2.o Ooura_fftsg2d_flt.o AddIce.o autoslic.o autosliccmd.o cfpix.o floatTIFF.o pdb.o slicelib.o MultisliceT.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS= -lm -lpthread -lfftw3 -lfftw3_threads -lfftw3f_threads -lfftw3f -lfftw3_omp -lfftw3f_omp

# $(ODIR)/%.o: %.c $(DEPS) $(DEQS)
# 	echo $(ODIR)
# 	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

#MsctKirkland

$(ODIR)/Oourafft.o: $(XArraydir)/Oourafft.cpp $(DEPS) $(DEQS)
	$(CC) -c $(CFLAGS1) -o $@ $< $(CFLAGS) 

$(ODIR)/Ooura_fftsg2d_flt.o: $(XArraydir)/Ooura_fftsg2d_flt.c $(DEPS) $(DEQS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS) 
	
$(ODIR)/Ooura_fftsg_flt.o: $(XArraydir)/Ooura_fftsg_flt.c $(DEPS) $(DEQS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg2d_dbl.o: $(XArraydir)/Ooura_fftsg2d_dbl.c $(DEPS) $(DEQS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg_dbl.o: $(XArraydir)/Ooura_fftsg_dbl.c $(DEPS) $(DEQS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/IXAHWave.o: $(XArraydir)/IXAHWave.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_ini.o: $(XArraydir)/XA_ini.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XArray.o: $(XArraydir)/XArray.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)
	
$(ODIR)/XA_file.o: $(XArraydir)/XA_file.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XAHWave.o: $(XArraydir)/XAHWave.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_spln2.o: $(XArraydir)/XA_spln2.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/fftwd2fc.o: $(XArraydir)/fftwd2fc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_fft2.o: $(XArraydir)/XA_fft2.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/AddIce.o: $(Temsimdir)/AddIce.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/autoslic.o: $(Temsimdir)/autoslic.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/autosliccmd.o: $(Temsimdir)/autosliccmd.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/cfpix.o: $(Temsimdir)/cfpix.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/floatTIFF.o: $(Temsimdir)/floatTIFF.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/pdb.o: $(Temsimdir)/pdb.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/slicelib.o: $(Temsimdir)/slicelib.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/MultisliceT.o: $(MSCTdir)/MultisliceT.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary $(CFLAGS)

MsctKirkland: $(OBJ)
	$(CPP) $(CFLAGS1) $(OBJ) $(CFLAGS) $(LIBS) -o MsctKirkland -lpthread

#clean:
	rm $(ODIR)/*.o



