CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra
LDFLAGS = -pthread -lmpfr -lgmp

# Source files
SOURCES = main.cpp kmer_sketch_tool_impl1.cpp kmer_sketch_tool_impl2.cpp MurmurHash3.cpp
HEADERS = kmer_sketch_tool.h sketch.h Newton_method_strong.h MurmurHash3.h mpreal.h
OBJECTS = $(SOURCES:.cpp=.o)

# Target executable
TARGET = repeat_robust_estimator

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJECTS) $(TARGET)

# Install (optional)
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# Test
test: $(TARGET)
	./test.sh

.PHONY: all clean install test