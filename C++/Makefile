CC		    = g++
CXXFLAGS    = -std=c++11 -Wall -o3 -flto -g -I/opt/homebrew/opt/boost/include 
LDLIBS		= -lboost_program_options #-lblas
LDFLAGS		= -L/opt/homebrew/opt/boost/lib
TARGET		= main
HDRS 		= Orbit.h

# CXXFLAGS    = -fopenmp -Wall -o3 -flto -g -I/opt/homebrew/opt/libomp/include
# LDFLAGS	  = -fopenmp -L/opt/homebrew/opt/libomp/lib



# Parameters of the test case specified in handout.
TEST_PARAMS 	= --SMA 6878 --e 1.3e-4 --i 58 --RAAN 290.5 --w 110 --theta 250

# Setting the "default" target. 
default: all
all: $(TARGET)

$(TARGET): $(TARGET).o Orbit.o

# A change in the header file requires re-compilation 
# (done as both .cpp include the header file).
%.o: %.cpp $(HDRS)

# These make targets don't create new files but rather run commands.
.PHONY: clean test

# Cleans working directory.
clean:
	rm -f $(TARGET) *.o


test: $(TARGET)
	export OMP_NUM_THREADS=1; ./$(TARGET) $(TEST_PARAMS)
