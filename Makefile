CUBATUREPATH = cubature
INCLUDE = -Isrc -I$(CUBATUREPATH)
LINK = -L. -llapack -lm
CFLAGS = -Wall -Werror -fno-exceptions -std=c11 -O3
OBJECTS = electrode.o auxiliary.o hcubature.o
CC = gcc # must be compliant to C Standard 99

.DEFAULT_GOAL := dynamic_library

.PHONY	:	clean cleanout dynamic_library

clean	:	cleanout
		rm -f *.exe *.so *.lib
cleanout	:
		rm -f *.o *.dat

hcubature.o	:	$(CUBATUREPATH)/hcubature.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c $(CUBATUREPATH)/hcubature.c $(LINK)

auxiliary.o	:	src/auxiliary.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/auxiliary.c $(LINK)

electrode.o	:	src/electrode.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/electrode.c $(LINK)

linalg.o	:	src/linalg.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src/linalg.c $(LINK)

libhem_electrode.so	:	$(OBJECTS)
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhem_electrode.so $(OBJECTS) $(LINK) -Wl,--out-implib,libhem_electrode.lib

libhem_linalg.so	:	linalg.o $(OBJECTS)
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhem_linalg.so linalg.o $(OBJECTS) $(LINK) -Wl,--out-implib,libhem_linalg.lib

dynamic_library	:	libhem_electrode.so libhem_linalg.so
		make cleanout

test	:	linalg.o $(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o test.exe src/testing.c linalg.o $(OBJECTS) $(LINK)

test_dynlib	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o test_dynlib.exe src/testing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg -lhem_electrode

timing	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o timing.exe examples/timing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg
# examples are named based on publication: author_volume_journal_issue
grcev20pwrd02	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev20pwrd02.exe examples/grcev20pwrd02.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg

grcev51emc03	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev51emc03.exe examples/grcev51emc03.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg

grcev12pwrd01	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev12pwrd01.exe examples/grcev12pwrd01.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg

condition_gs_layered	:	libhem_electrode.so libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o condition_gs_layered.exe examples/condition_gs_layered.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_linalg
