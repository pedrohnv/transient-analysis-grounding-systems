#MKLROOT is defined by a bash script that is shipped together with intel's MKL.
#Run in the terminal "source mklvars.sh intel64"

#The shared libraries slatec and lapack  are assumed to be in your path.
#If they are not, append their location to the INCLUDE variable or
#build and install them with "make slatec" and make "make lapack"

INTELLINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
INTELFLAGS = -DMKL_ILP64 -m64
CUBATUREPATH = cubature/
SLATECPATH = slatec/
WOLFRAMPATH = /usr/local/Wolfram/Mathematica/11.3/SystemFiles/IncludeFiles/C/
INCLUDE = -Isrc -I$(CUBATUREPATH)
LINK = -L. $(INTELLINK) -lslatec -lgfortran
CFLAGS = -Wall -fno-exceptions -std=c11 $(INTELFLAGS) -Werror -O3
#CFLAGS = -Wall -std=c11 $(INTELFLAGS) -g
OBJECTS = electrode.o auxiliary.o hcubature.o

#Compilers
CC = gcc
FC = gfortran

.DEFAULT_GOAL := dynamicLibrary

.PHONY	:	clean cleanout allclean slatec libhem.so dynamicLibrary wolfram

clean	:
		rm -f *.o *.dat examples/*.o examples/*.dat examples/grcev12pwrd01*.dat

cleanout	:
		rm -f *.a *.so examples/*.a interfaces/mathematica/*.so

cleanall	:
		make clean
		make cleanout
		cd $(SLATECPATH) && $(MAKE) clean
		cd $(CUBATUREPATH) && $(MAKE) clean

hcubature.o	:	$(CUBATUREPATH)hcubature.c
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c $(CUBATUREPATH)hcubature.c $(LINK)

auxiliary.o	:	auxiliary.c
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c auxiliary.c $(LINK)

electrode.o	:	electrode.c auxiliary.o hcubature.o
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c electrode.c $(LINK)

dynamicLibrary	:
		$(CC) -fPIC -shared $(CFLAGS) $(INCLUDE) -o libhem.so electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

test	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o testing.a testing.c $(OBJECTS) $(LINK)

testDLL	:	libhem.so
		$(CC) $(CFLAGS) $(INCLUDE) -o testingDLL.a testing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem

slatec	:
		cd $(SLATECPATH) && $(MAKE) FC=$(FC) all

wolfram	:
		$(CC) -fPIC -shared -std=c11 $(INTELFLAGS) -Werror -O3 $(INCLUDE) -I$(WOLFRAMPATH) -o interfaces/mathematica/libhem_mma.so interface_wolfram.c electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

timing	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o timing.a examples/timing.c $(OBJECTS) $(LINK)
# examples are named based on publication: author_volume_journal_issue
grcev20pwrd02	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev20pwrd02.a examples/grcev20pwrd02.c $(OBJECTS) $(LINK)

grcev51emc03	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev51emc03.a examples/grcev51emc03.c $(OBJECTS) $(LINK)

grcev12pwrd01	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/grcev12pwrd01/grcev12pwrd01.a examples/grcev12pwrd01/grcev12pwrd01.c $(OBJECTS) $(LINK)
