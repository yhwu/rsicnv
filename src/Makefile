EXECUTABLE = rsicnv
SAMTOOLS = samtools-0.1.18
ALGLIB = alglib
LOBJS = ${ALGLIB}/libalglib.a ${SAMTOOLS}/libbam.a
CC = g++
CFLAGS = -g -O2 -Wall
#CFLAGS = -g -Wall
INCLUDE = -I${SAMTOOLS} -I${ALGLIB}
LIBS =  -L${SAMTOOLS} -lbam -L${ALGLIB} -lalglib -lz

MAINCXX = rsi.cpp pairrd.cpp plotcnv.cpp wufunctions.cpp loaddata.cpp gccontent.cpp readref.cpp samfunctions.cpp 
MAINHDR = wu2.h pstream.h $(MAINCXX:.cpp=.h)
MAINOBJ = $(MAINCXX:.cpp=.o)

all : ${EXECUTABLE}

rsicnv : $(LOBJS) ${MAINOBJ}
	${CC} ${CFLAGS} ${MAINOBJ} ${LIBS} -o $@

$(MAINOBJ): %.o: %.cpp %.h wu2.h
	${CC} ${INCLUDE} ${CFLAGS} -c $< -o $@

${ALGLIB}/libalglib.a : 
	cd ${ALGLIB} && make

${SAMTOOLS}/libbam.a :
	cd ${SAMTOOLS} && make libbam.a

clean : 
	cd ${SAMTOOLS} && make clean
	cd ${ALGLIB} && make clean
	rm -fr *.o ${EXECUTABLE}
