#MKLROOT is defined by a bash script that is shipped together with intel's MKL
INTELLINK = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
INTELFLAGS = -DMKL_ILP64 -m64
CUBATUREPATH = ./cubature/
INCLUDE = -I. -I$(CUBATUREPATH)
LINK = -L. $(INTELLINK) -lslatec -llapack
CFLAGS = -Wall -fno-exceptions -std=c11 $(INTELFLAGS) -Werror -g
OBJECTS = Electrode.o auxiliary.o hcubature.o

hcubature.o	:	$(CUBATUREPATH)hcubature.c
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c $(CUBATUREPATH)hcubature.c $(LINK)
auxiliary.o	:	auxiliary.c
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c auxiliary.c $(LINK)
Electrode.o	:	Electrode.c auxiliary.o hcubature.o
		gcc -fPIC $(CFLAGS) $(INCLUDE) -c Electrode.c $(LINK)
hem_c.so	:	Electrode.c auxiliary.c hcubature.c
		gcc $(CFLAGS) $(INCLUDE) -shared -o hem.so Electrode.c auxiliary.c $(CUBATUREPATH)hcubature.c $(LINK)
test	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o testing.a testing.c $(OBJECTS) $(LINK)
.PHONY	:	clean
clean	:
		rm -f *.[ao] hem.so $(CUBATUREPATH)*.o *.dat examples/*.dat
example20pwrd02grcev	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o 20pwrd02grcev.a examples/20pwrd02grcev.c $(OBJECTS) $(LINK)
example51emc03grcev	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o 51emc03grcev.a examples/51emc03grcev.c $(OBJECTS) $(LINK)
timing	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o timing.a examples/timing.c $(OBJECTS) $(LINK)
miranda64	:	$(OBJECTS)
		gcc $(CFLAGS) $(INCLUDE) -o examples/miranda64.a examples/grounding_miranda64.c $(OBJECTS) $(LINK)
