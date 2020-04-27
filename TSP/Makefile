CXX = g++
CFLAGS = -std=c++11 -O3

INCLUDES = -I include

NAME := TSP

SRC := src
OBJ := obj

SOURCES := $(wildcard $(SRC)/*.cpp)
OBJECTS := $(patsubst $(SRC)/%.cpp, $(OBJ)/%.o, $(SOURCES))

all: TSP_BUILD

TSP_BUILD: $(OBJECTS)
	$(CXX) $(CFLAGS) $(INCLUDES) -o $(NAME) $(OBJECTS)
	@echo compiled without errors

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) -r $(OBJ)
	$(RM) $(NAME)

$(shell mkdir -p $(OBJ))
