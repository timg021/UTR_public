CC=gcc-9
CPP=g++-9
CFLAGS= -O3 #-g3
CFLAGS1= -Wno-unused-result -Wno-format-overflow -std=gnu++17
CFLAGSC= -Wno-unused-result -Wno-format-overflow

#User will need to modify these for their own system
#FFTWdir, must contain libfftw3
FFTWdir = /usr/lib/x86_64-linux-gnu

ODIR = obj
LDIR = $(FFTWdir)
XArraydir = ../XArrayLibrary
Temsimdir = ../XArrayLibrary
MSCTdir = ../MsctKirkland

_DEPS = XArray1D.h XA_spln2.h IXAHWave.h XA_data.h XA_head1.h XA_ini.h XArray2D.h OouraFft.h stdafx.h XA_fft2.h XA_head2.h XArray3D.h XA_tiff.h stdafx.h XA_move2.h XA_nrrand.h XArray.h IXAHead.h XA_file.h XAHWave.h
DEPS = $(patsubst %,$(XArraydir)/%,$(_DEPS))

_DEQS = AddIce.h autoslic.hpp autosliccmd.h cfpix.hpp floatTIFF.hpp pdb.h slicelib.hpp
DEQS = $(patsubst %,$(Temsimdir)/%,$(_DEQS))

_OBJ =  Ooura_fftsg_flt.o Ooura_fftsg2d_dbl.o Ooura_fftsg_dbl.o XArray.o IXAHWave.o XA_ini.o XA_file.o XAHWave.o Ooura_fftsg2d_flt.o AddIce.o autoslic.o autosliccmd.o cfpix.o floatTIFF.o pdb.o slicelib.o MultisliceK.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS= -lm -lfftw3 -lfftw3f

# $(ODIR)/%.o: %.c $(DEPS) $(DEQS)
# 	echo $(ODIR)
# 	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

#MSCT

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

$(ODIR)/MultisliceK.o: $(MSCTdir)/MultisliceK.cpp $(DEPS) $(DEQS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary $(CFLAGS)

MSCT: $(OBJ)
	$(CPP) $(CFLAGS1) $(OBJ) $(CFLAGS) $(LIBS) -o MsctKirkland -lpthread

#clean:
	rm $(ODIR)/*.o



