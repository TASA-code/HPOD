
CC         = g++
CXXFLAGS   = -std=c++17 -Wall -O3 -flto -g -I/opt/homebrew/opt/boost/include -I/opt/homebrew/opt/eigen/include/eigen3

LDLIBS     = -lboost_program_options
LDFLAGS    = -L/opt/homebrew/opt/boost/lib 
TARGET     = main
HDRS       = Orbit.h Coordinate.h

# List of source files
SRCS       = main.cpp Orbit.cpp Coordinate.cpp

# List of object files
OBJS       = $(SRCS:.cpp=.o)

# Parameters of the test case specified in handout.
TEST_PARAMS = --SMA 6945 --e 1.3e-4 --i 60 --RAAN 290.5 --w 110 --theta 250

# Setting the "default" target.
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

# A change in the header file requires re-compilation
# (done as both .cpp include the header file).
%.o: %.cpp $(HDRS)
	$(CC) $(CXXFLAGS) -c $< -o $@

# Clean up intermediate files and the final executable.
clean:
	rm -f $(OBJS) $(TARGET)

# Run the test case
test: $(TARGET)
	./$(TARGET) $(TEST_PARAMS)

