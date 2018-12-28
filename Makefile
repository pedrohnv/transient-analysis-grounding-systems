#MKLROOT is defined by a bash script that is shipped together with intel's MKL
#source mklvars.sh intel64
INTELLINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
INTELFLAGS = -DMKL_ILP64 -m64
CUBATUREPATH = cubature/
SLATECPATH = slatec/
LINK = -L. $(INTELLINK) -lslatec -llapack
CFLAGS = -Wall -fno-exceptions -std=c11 $(INTELFLAGS) -Werror -O3
#CFLAGS = -Wall -std=c11 $(INTELFLAGS) -g
OBJECTS = Electrode.o auxiliary.o hcubature.o

.DEFAULT_GOAL := dynamicLibrary

.PHONY	:	clean slatec

clean	:
		rm -f *.[ao] *.so *.dat examples/*.dat examples/*.[ao]
		cd $(SLATECPATH) && $(MAKE) clean
		cd $(CUBATUREPATH) && $(MAKE) clean

hcubature.o	:	$(CUBATUREPATH)hcubature.c
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c $(CUBATUREPATH)hcubature.c $(LINK)

auxiliary.o	:	auxiliary.c
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c auxiliary.c $(LINK)

Electrode.o	:	Electrode.c auxiliary.o hcubature.o
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c Electrode.c $(LINK)

dynamicLibrary	:	Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c
		gcc -fPIC -shared $(CFLAGS) $(INCLUDE) -o libhem.so Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)

test	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o testing.a testing.c $(OBJECTS) $(LINK)

timing	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o timing.a examples/timing.c $(OBJECTS) $(LINK)

example20pwrd02grcev	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o 20pwrd02grcev.a examples/20pwrd02grcev.c $(OBJECTS) $(LINK)

example51emc03grcev	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o 51emc03grcev.a examples/51emc03grcev.c $(OBJECTS) $(LINK)

example12pwrd01grcev	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o 12pwrd01grcev.a examples/ex12pwrd01grcev/12pwrd01grcev.c $(OBJECTS) $(LINK)

exampleAlipioSchroederRCA	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o AlipioSchroederRCA.a examples/AlipioSchroederRCA.c $(OBJECTS) $(LINK)

miranda64	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o examples/miranda64.a examples/grounding_miranda64.c $(OBJECTS) $(LINK)

malha01	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o examples/malha01.a examples/malha01.c $(OBJECTS) $(LINK)

slatec	:
		cd $(SLATECPATH) && $(MAKE) FC=gfortran all
