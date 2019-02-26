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
INCLUDE = -I. -I$(CUBATUREPATH)
LINK = -L. $(INTELLINK) -lslatec -lgfortran
CFLAGS = -Wall -fno-exceptions -std=c11 $(INTELFLAGS) -Werror -O3
#CFLAGS = -Wall -std=c11 $(INTELFLAGS) -g
OBJECTS = Electrode.o auxiliary.o hcubature.o

#Compilers
CC = gcc
FC = gfortran

.DEFAULT_GOAL := dynamicLibrary

.PHONY	:	clean cleanout allclean slatec libhem.so dynamicLibrary wolfram

clean	:
		rm -f *.o *.dat examples/*.o examples/*.dat examples/ex12pwrd01grcev*.dat

cleanout	:
		rm -f *.a *.so examples/*.a interfaces/mathematica/*.so

allclean	:
		make clean
		make cleanout
		cd $(SLATECPATH) && $(MAKE) clean
		cd $(CUBATUREPATH) && $(MAKE) clean

hcubature.o	:	$(CUBATUREPATH)hcubature.c
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c $(CUBATUREPATH)hcubature.c $(LINK)

auxiliary.o	:	auxiliary.c
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c auxiliary.c $(LINK)

Electrode.o	:	Electrode.c auxiliary.o hcubature.o
		$(CC) -fPIC $(CFLAGS) $(INCLUDE) -c Electrode.c $(LINK)

libhem.so	:
		make dynamicLibrary

dynamicLibrary	:
		$(CC) -fPIC -shared $(CFLAGS) $(INCLUDE) -o libhem.so Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

test	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o testing.a testing.c $(OBJECTS) $(LINK)

testDLL	:	libhem.so
		$(CC) $(CFLAGS) $(INCLUDE) -o testingDLL.a testing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem

slatec	:
		cd $(SLATECPATH) && $(MAKE) FC=$(FC) all

wolfram	:
		$(CC) -fPIC -shared -std=c11 $(INTELFLAGS) -Werror -O3 $(INCLUDE) -I$(WOLFRAMPATH) -Iinterfaces/mathematica -o interfaces/mathematica/libhem_mma.so interfaces/mathematica/InterfaceWolfram.c Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

timing	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o timing.a examples/timing.c $(OBJECTS) $(LINK)

ex20pwrd02grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o ex20pwrd02grcev.a examples/ex20pwrd02grcev.c $(OBJECTS) $(LINK)

ex51emc03grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o ex51emc03grcev.a examples/ex51emc03grcev.c $(OBJECTS) $(LINK)

ex12pwrd01grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/ex12pwrd01grcev/ex12pwrd01grcev.a examples/ex12pwrd01grcev/ex12pwrd01grcev.c $(OBJECTS) $(LINK)

exAlipioSchroederRCA	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o AlipioSchroederRCA.a examples/AlipioSchroederRCA.c $(OBJECTS) $(LINK)

exMiranda64	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/miranda64.a examples/grounding_miranda64.c $(OBJECTS) $(LINK)

exMalha01	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/malha01.a examples/malha01.c $(OBJECTS) $(LINK)
