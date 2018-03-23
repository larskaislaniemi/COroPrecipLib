
CC=mpicc
CFLAGS=-I$(INCDIR) -mkl=cluster -Wall -std=c99
DEPS=$(INCDIR)/oropreciplib.h
LIBOUT=lib
BINOUT=bin
SRCDIR=src
INCDIR=inc

oroprecip: $(DEPS)
	mkdir -p $(BINOUT)
	$(CC) -DOROPRECIP_STANDALONE=1 $(CFLAGS) -o bin/$@ $(SRCDIR)/*.c

oroprecip.so: $(DEPS)
	mkdir -p $(LIBOUT)
	$(CC) -DOROPRECIP_STANDALONE=0 $(CFLAGS) -shared -fPIC -o lib/$@ $(SRCDIR)/*.c

oroprecip.a: $(DEPS)
	mkdir -p $(LIBOUT)
	$(CC) -DOROPRECIP_STANDALONE=0 $(CFLAGS) -c $(SRCDIR)/*.c
	ar rcs lib/$@ *.o

clean:
	rm -fr bin lib *.o
