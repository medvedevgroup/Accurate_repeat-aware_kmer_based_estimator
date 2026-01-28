CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
TARGET = mutation_estimator

all: $(TARGET)

$(TARGET): mutation_estimator.cpp
	$(CXX) $(CXXFLAGS) mutation_estimator.cpp -o $(TARGET)

clean:
	rm -f $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

.PHONY: all clean install