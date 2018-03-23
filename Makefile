
CC=mpicc
FC=ifort
CFLAGS=-I$(INCDIR) -mkl=cluster -Wall -std=c99
DEPS=$(INCDIR)/oropreciplib.h
LIBOUT=lib
BINOUT=bin
SRCDIR=src
INCDIR=inc

all: oroprecip liboroprecip.so liboroprecip.a m_oroprecip.mod
	mv *.o lib/

alllib: liboroprecip.so liboroprecip.a m_oroprecip.mod
	mv *.o lib/

oroprecip: $(DEPS)
	mkdir -p $(BINOUT)
	$(CC) -DOROPRECIP_STANDALONE=1 $(CFLAGS) -o bin/$@ $(SRCDIR)/*.c

liboroprecip.so: $(DEPS)
	mkdir -p $(LIBOUT)
	$(CC) -DOROPRECIP_STANDALONE=0 $(CFLAGS) -shared -fPIC -o lib/$@ $(SRCDIR)/*.c

liboroprecip.a: $(DEPS)
	mkdir -p $(LIBOUT)
	$(CC) -DOROPRECIP_STANDALONE=0 $(CFLAGS) -c $(SRCDIR)/*.c
	ar rcs lib/$@ *.o

m_oroprecip.mod: $(DEPS)
	mkdir -p $(LIBOUT)
	$(FC) -c $(SRCDIR)/oroprecip_iface.f90
	mv *.mod lib/

clean:
	rm -fr bin lib *.o *.mod
