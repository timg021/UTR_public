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
LDIR = -L $(FFTWdir)
XArraydir = ../XArrayLibrary
UTRdir = ../PhaseRetrieval
_DEPS = XA_filt2.h XArray1D.h XA_spln3.h XA_fft2.h fftwd3drc.h IXAHWave.h XA_data.h XA_head1.h XA_ini.h XArray2D.h XA_tie.h OouraFft.h stdafx.cpp XA_fft2.h XA_head2.h XA_iwfr.h XArray3D.h XA_tiff.h fftwd3frc.h fftwd2fc.h stdafx.h XA_fftr2.h XA_head3.h XA_move2.h XA_move3.h XArray.h IXAHead.h XA_file.h XAHWave.h XA_spln2.h XA_DICOM.h UTR.h FindPeaks.h
DEPS = $(patsubst %,$(XArraydir)/%,$(_DEPS))

_OBJ =  Oourafft.o Ooura_fftsg_flt.o Ooura_fftsg_dbl.o Ooura_fftsg2d_flt.o Ooura_fftsg2d_dbl.o XArray.o fftwd2fc.o fftwd3drc.o fftwd3frc.o XA_spln2.o XA_fft2.o IXAHWave.o XA_ini.o XA_file.o XAHWave.o FindPeaks.o XA_DICOM.o UTR.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS= -lm -lpthread -lfftw3 -lfftw3_threads -lfftw3f_threads -lfftw3f -lfftw3_omp -lfftw3f_omp -lgdcmCommon -lgdcmDICT -lgdcmDSED -lgdcmIOD -lgdcmMEXD -lgdcmMSFF -lgdcmCommon

# $(ODIR)/%.o: %.c $(XarrayDEPS)
# 	echo $(ODIR)
# 	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

#UTR


$(ODIR)/Oourafft.o: $(XArraydir)/Oourafft.cpp $(DEPS)
	$(CC) -c $(CFLAGS1) -o $@ $< $(CFLAGS) 

$(ODIR)/Ooura_fftsg2d_flt.o: $(XArraydir)/Ooura_fftsg2d_flt.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS) 
	
$(ODIR)/Ooura_fftsg_flt.o: $(XArraydir)/Ooura_fftsg_flt.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg2d_dbl.o: $(XArraydir)/Ooura_fftsg2d_dbl.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/Ooura_fftsg_dbl.o: $(XArraydir)/Ooura_fftsg_dbl.c $(DEPS)
	$(CC) -c $(CFLAGSC) -o $@ $< $(CFLAGS)  

$(ODIR)/fftwd2fc.o: $(XArraydir)/fftwd2fc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/fftwd3drc.o: $(XArraydir)/fftwd3drc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/fftwd3frc.o: $(XArraydir)/fftwd3frc.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_spln2.o: $(XArraydir)/XA_spln2.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< $(CFLAGS)

$(ODIR)/XA_fft2.o: $(XArraydir)/XA_fft2.cpp $(DEPS)
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

$(ODIR)/XA_DICOM.o: $(XArraydir)/XA_DICOM.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I /usr/include/gdcm-3.0 $(CFLAGS)

$(ODIR)/FindPeaks.o: $(UTRdir)/FindPeaks.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary $(CFLAGS)

$(ODIR)/UTR.o: $(UTRdir)/UTR.cpp $(DEPS)
	$(CPP) -c $(CFLAGS1) -o $@ $< -I ../XArrayLibrary -I /usr/include/gdcm-3.0 $(CFLAGS)

UTR: $(OBJ)
	$(CPP) $(CFLAGS1) $(OBJ) $(CFLAGS) $(LDIR) $(LIBS) -o UTR

clean:
	rm $(ODIR)/*.o


