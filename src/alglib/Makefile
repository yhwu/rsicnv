CC=			g++
CFLAGS=		-g -Wall -O2
LOBJS= alglibinterface.o alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o

PROG=		libalglib.a
INCLUDES=	-I.
LIBPATH=

.SUFFIXES:.cpp .o

.cpp.o:
		$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

all:libalglib.a

lib:libalglib.a

libalglib.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

clean:
	rm -fr *.o *.a
