CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic

SRC = seqGMM.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = seqGMM

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)