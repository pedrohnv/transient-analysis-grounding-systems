#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
HEMDIR = /home/pedro/codigos/HP_HEM
INCLUDE = -I. -I/usr/include/octave-4.2.2/ -Iswiglib -I$(HEMDIR)
octave	:	InterfaceOctave.cpp
		swig -octave -c++ -o hem_octave.cpp OctaveHEM.i
		mkoctfile $(INCLUDE) hem_octave.cpp InterfaceOctave.cpp '-Wl,-rpath,$$ORIGIN' -L. -lhem
clean	:
		rm -f *.o *.oct hem_octave.cpp
