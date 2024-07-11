# Variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGETS = cubic_min.exe quad_min.exe
SRCS_CUBIC = cubic_interp.cpp
SRCS_QUADRATIC = quad_interp.cpp
OBJS_CUBIC = $(SRCS_CUBIC:.cpp=.o)
OBJS_QUADRATIC = $(SRCS_QUADRATIC:.cpp=.o)

# Default target
all: $(TARGETS)

# Linking the cubic interpolation program
cubic_min.exe: $(OBJS_CUBIC)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Linking the quadratic interpolation program
quad_min.exe: $(OBJS_QUADRATIC)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compiling the source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJS_CUBIC) $(OBJS_QUADRATIC) $(TARGETS)

# Phony targets
.PHONY: all clean
