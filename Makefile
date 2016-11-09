#
# Enable the right compiler options below
# And then run "make" to see her compile
#

CC=g++ # Local Machine
#CC=pgcpp  # HECToR
CFLAGS=-c -g -O3 -Wall -Werror -pedantic -Wno-long-long # -fopenmp Local Machine
#CFLAGS=-c -O3 --pedantic #HECToR
LDFLAGS=-lm
LIBRARIES= # These are mpic++ or g++ flags: -fopenmp 
SOURCES=Functions.cpp \
        Gamess_UK.cpp \
        Genetic.cpp \
        Gradients.cpp \
        History.cpp \
        IO.cpp \
	Linear.cpp \
        Main.cpp \
        Newton_Raphson.cpp \
        Nwchem.cpp \
        Outputs.cpp \
        Powells.cpp \
        Punch.cpp \
        Utils.cpp 


OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fit_my_ecp

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBRARIES)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ 

clean:
	rm -f ${OBJECTS} ${EXECUTABLE}
