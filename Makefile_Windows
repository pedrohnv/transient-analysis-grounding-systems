#MKLROOT is defined by a bash script that is shipped together with intel's MKL
#source mklvars.sh intel64
#manual set:
#"C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.1.144\windows\mkl\bin\mklvars.bat" intel64
MKLROOT = "C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries_2019.1.144\\windows\\mkl"
INTELLINK = "C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries_2019.1.144\\windows\\mkl\\lib\\intel64_win\\mkl_rt.lib"
INTELFLAGS = -DMKL_ILP64
#The shared library slatec is assumed to be in your path.
#If it is not, append its location to the INCLUDE variable or
#build and install it with "make slatec"

CUBATUREPATH = cubature
SLATECPATH = slatec
WOLFRAMPATH = "C:\\Program Files\\Wolfram Research\\Mathematica\\11.3\\SystemFiles\\IncludeFiles\\C"
INCLUDE = -Isrc -I$(CUBATUREPATH)
LINK = -L.
CFLAGS = -Wall -Werror -fno-exceptions -std=c11 -g
OBJECTS = electrode.o auxiliary.o hcubature.o
#Compilers
CC = gcc # must be compliant to C Standard 99
FC = gfortran

.DEFAULT_GOAL := dynamic_library

.PHONY	:	clean cleanout dynamicLibrary

clean	:	cleanout
		rm -f *.exe *.dll *.lib
cleanout	:
		rm -f *.o *.dat

hcubature.o	:	$(CUBATUREPATH)\\hcubature.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c $(CUBATUREPATH)\\hcubature.c $(LINK)

auxiliary.o	:	src\\auxiliary.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src\\auxiliary.c $(LINK)

electrode.o	:	src\\electrode.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -c src\\electrode.c $(LINK)

linalg.o	:	src\\linalg.c
		$(CC) $(CFLAGS) -fPIC $(INCLUDE) -I$(MKLROOT)\\include -c src\\linalg.c $(LINK)

libhem_electrode.dll	:	$(OBJECTS)
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -o libhem_electrode.dll $(OBJECTS) $(LINK) -Wl,--out-implib,libhem_electrode.lib

libhem_linalg.dll	:	linalg.o
		$(CC) $(CFLAGS) $(INTELFLAGS) -fPIC -shared $(INCLUDE) -I$(MKLROOT)\\include -o libhem_linalg.dll linalg.o $(OBJECTS) $(LINK) $(INTELLINK) -Wl,--out-implib,libhem_linalg.lib

dynamic_library	:	libhem_electrode.dll libhem_linalg.dll
		make cleanout

test	:	linalg.o $(OBJECTS)
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o test.exe src\\testing.c linalg.o $(OBJECTS) $(LINK) $(INTELLINK)

test_dynlib	:	libhem_electrode.dll libhem_linalg.dll
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o test_dynlib.exe src\\testing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_electrode -lhem_linalg

wolfram_electrode	:
		$(CC) $(CFLAGS) -fPIC -shared $(INCLUDE) -I$(WOLFRAMPATH) -o libhem_electrode_mma.dll interfaces\\mathematica\\interface_wolfram.c src\\electrode.c src\\auxiliary.c $(CUBATUREPATH)\\hcubature.c $(LINK)

slatec	:
		cd $(SLATECPATH) && $(MAKE) FC=$(FC) all

timing	:	libhem_electrode.dll libhem_linalg.dll
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o timing.exe examples\\timing.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_electrode -lhem_linalg
# examples are named based on publication: author_volume_journal_issue
grcev20pwrd02	:	libhem_electrode.dll libhem_linalg.dll
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o grcev20pwrd02.exe examples\\grcev20pwrd02.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_electrode -lhem_linalg

grcev51emc03	:	libhem_electrode.dll libhem_linalg.dll
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o grcev51emc03.exe examples\\grcev51emc03.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_electrode -lhem_linalg

grcev12pwrd01	:	libhem_electrode.dll libhem_linalg.dll
		$(CC) $(CFLAGS) $(INCLUDE) -I$(MKLROOT)\\include -o grcev12pwrd01.exe examples\\grcev12pwrd01.c '-Wl,-rpath,$$ORIGIN' $(LINK) -lhem_electrode -lhem_linalg
