CXX = g++
CXXP = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic
CXXPFLAGS = -std=c++14 -fvisibility=hidden -lpthread -Wall -Wextra

SRC = seqGMM.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = seqGMM

SRC2 = parallelGMM.cpp
OBJ2 = $(SRC2:.cpp=.o)
TARGET2 = parallelGMM

.PHONY: all clean

all: seqGMM parallelGMM

seqGMM: $(OBJ)
	$(CXX) $(CXXFLAGS) seqGMM.o -o $@

parallelGMM:
	$(CXXP) -o $@ $(CXXPFLAGS) parallelGMM.cpp

seqGMM.o: seqGMM.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET) $(OBJ2) $(TARGET2)