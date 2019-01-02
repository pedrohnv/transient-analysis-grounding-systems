#MKLROOT is defined by a bash script that is shipped together with intel's MKL.
#Run in the terminal "source mklvars.sh intel64"

#The shared libraries slatec and lapack  are assumed to be in your path.
#If they are not, append their location to the INCLUDE variable or
#build and install them with "make slatec" and make "make lapack"

INTELLINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
INTELFLAGS = -DMKL_ILP64 -m64
CUBATUREPATH = cubature/
SLATECPATH = slatec/
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
		rm -f *.o *.dat examples/*.o examples/*.dat

cleanout	:
		rm -f *.a *.so examples*.a

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
		$(CC) $(CFLAGS) $(INCLUDE) -o testingDLL.a testing.c '-Wl,-rpath,$$ORIGIN' $(LINK) libhem.so

slatec	:
		cd $(SLATECPATH) && $(MAKE) FC=$(FC) all

wolfram	:
		$(CC) -fPIC -shared -std=c11 $(INTELFLAGS) -Werror -O3 $(INCLUDE) -I$(WOLFRAMPATH) -o libhem_mma.so InterfaceWolfram.c Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

timing	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o timing.a examples/timing.c $(OBJECTS) $(LINK)

example20pwrd02grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o 20pwrd02grcev.a examples/20pwrd02grcev.c $(OBJECTS) $(LINK)

example51emc03grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o 51emc03grcev.a examples/51emc03grcev.c $(OBJECTS) $(LINK)

example12pwrd01grcev	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o 12pwrd01grcev.a examples/ex12pwrd01grcev/12pwrd01grcev.c $(OBJECTS) $(LINK)

exampleAlipioSchroederRCA	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o AlipioSchroederRCA.a examples/AlipioSchroederRCA.c $(OBJECTS) $(LINK)

miranda64	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/miranda64.a examples/grounding_miranda64.c $(OBJECTS) $(LINK)

malha01	:	$(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o examples/malha01.a examples/malha01.c $(OBJECTS) $(LINK)
