# compiler must be compliant to C Standard 99
CC = gcc
CUBATUREPATH = cubature

INCLUDE = -Isrc -I$(CUBATUREPATH)
CFLAGS = -Wall -fno-exceptions -Wfatal-errors -std=c11 -O3 -m64 -fopenmp
LINK = -L. -llapack -lblas -lgfortran -lpthread -lfftw3 -lm -ldl

## Intel MKL version:
#INCLUDE = -Isrc -I$(CUBATUREPATH) -I${MKLROOT}/include
#CFLAGS = -Wall -fno-exceptions -Wfatal-errors -std=c11 -O3 -DMKL_ILP64 -DMKL_Complex16="_Complex double" -m64 -fopenmp
#LINK = -L. -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lfftw3 -lm -ldl

OBJECTS = electrode.o auxiliary.o hcubature.o grid.o linalg.o

.DEFAULT_GOAL := dynamic_library
.PHONY:	clean cleanout dynamic_library

clean	:	cleanout
		rm -f *.exe *.so *.lib
cleanout	:
		rm -f *.o *.csv *.gif *.png

# Objects
hcubature.o	:	$(CUBATUREPATH)/hcubature.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c $(CUBATUREPATH)/hcubature.c $(LINK)

auxiliary.o	:	src/auxiliary.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/auxiliary.c $(LINK)

electrode.o	:	src/electrode.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/electrode.c $(LINK)

grid.o	:	src/grid.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/grid.c $(LINK)

linalg.o	:	src/linalg.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/linalg.c $(LINK)

test	:	$(OBJECTS)
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -Isrc/test -o test.exe src/test/testing.c $(OBJECTS) $(LINK)

# Libraries
dynamic_library	:	libhphem.so

libhphem.so	:	$(OBJECTS)
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhphem.so $(OBJECTS)

# Examples are named based on publication: author_volume_journal_issue
grcev12pwrd01	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev12pwrd01.exe examples/grcev12pwrd01.c $(OBJECTS) $(LINK)

visacro57emc01	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o visacro57emc01.exe examples/visacro57emc01.c $(OBJECTS) $(LINK)

sunjerga173powsys	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o sunjerga173powsys.exe examples/sunjerga173powsys.c $(OBJECTS) $(LINK)

alipio83powsys	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o alipio83powsys.exe examples/alipio83powsys.c $(OBJECTS) $(LINK)
