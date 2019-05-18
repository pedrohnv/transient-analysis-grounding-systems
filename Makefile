CUBATUREPATH = cubature
INCLUDE = -Isrc -I$(CUBATUREPATH)
LINK = -L. -llapack -lblas -lm
#CFLAGS = -Wall -Werror -fno-exceptions -std=c11 -O3
CFLAGS = -Wall -fno-exceptions -std=c11 -g
OBJECTS = electrode.o auxiliary.o hcubature.o
CC = gcc # must be compliant to C Standard 99

.DEFAULT_GOAL := program

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
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhem_electrode.so $(OBJECTS) $(LINK)

libhem_linalg.so	:	linalg.o $(OBJECTS)
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhem_linalg.so linalg.o $(OBJECTS) $(LINK)

dynamic_library	:	libhem_electrode.so libhem_linalg.so

test	:	linalg.o $(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -o test.exe src/testing.c linalg.o $(OBJECTS) $(LINK)

test_dynlib	:	libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o test_dynlib.exe src/testing.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg

timing	:	libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o timing.exe examples/timing.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg
# examples are named based on publication: author_volume_journal_issue
grcev20pwrd02	:	libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev20pwrd02.exe examples/grcev20pwrd02.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg

grcev51emc03	:	libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev51emc03.exe examples/grcev51emc03.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg

grcev12pwrd01	:	libhem_linalg.so
		$(CC) $(CFLAGS) $(INCLUDE) -o grcev12pwrd01.exe examples/grcev12pwrd01.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg

program	:	dynamic_library
		$(CC) $(CFLAGS) $(INCLUDE) -o hp_hem.exe src/program.c $(LINK) '-Wl,-rpath,$$ORIGIN' -lhem_linalg
