# Makefile adapted from http://www.partow.net/programming/makefile/index.html

CXX		:= -c++
CXXFLAGS	:= -pedantic-errors -Wall -Wextra #-Werror
LDFLAGS		:= -L/usr/lib -lboost_system -lboost_program_options -lboost_filesystem -lboost_date_time -lstdc++ -lm
BUILD		:= ./build
OBJ_DIR		:= $(BUILD)/objects
SRC_DIR		:= ./src
APP_DIR		:= $(BUILD)/apps
# TARGET		:= program
INCLUDE		:= -Iinclude/
#SRC		:= $(wildcard src/*.cpp)
#SRC		:= $(shell find $(SRC_DIR) -type f -name '*.cpp')
#OBJECTS		:= $(SRC:%.cpp=$(OBJ_DIR)/%.o)
#OBJECTS		:= $(patsubst ./src/%, $(BUILD)/%, $(SRC:.cpp=.o))
TF_SRC		:= $(shell find ./src/core ./src/telomereFinder -type f -name '*.cpp') ./src/telomereFinder_driver.cpp
TF_OBJ		:= $(patsubst ./src/%, $(OBJ_DIR)/%, $(TF_SRC:.cpp=.o))
JF_SRC		:= $(shell find ./src/core ./src/junctionFinder -type f -name '*.cpp') ./src/junctionFinder_driver.cpp
JF_OBJ		:= $(patsubst ./src/%, $(OBJ_DIR)/%, $(JF_SRC:.cpp=.o))
SQ_SRC		:= $(shell find ./src/core -type f -name '*.cpp') ./src/sequenceQuality_driver.cpp
SQ_OBJ		:= $(patsubst ./src/%, $(OBJ_DIR)/%, $(SQ_SRC:.cpp=.o))
WCD_SRC		:= $(shell find ./src/core -type f -name '*.cpp') ./src/wcdInterrogate_driver.cpp
WCD_OBJ		:= $(patsubst ./src/%, $(OBJ_DIR)/%, $(WCD_SRC:.cpp=.o))


all: build telomereFinder junctionFinder sequenceQuality wcdInterrogate

telomereFinder: $(APP_DIR)/telomereFinder
junctionFinder: $(APP_DIR)/junctionFinder
sequenceQuality: $(APP_DIR)/sequenceQuality
wcdInterrogate: $(APP_DIR)/wcdInterrogate

$(APP_DIR)/telomereFinder: $(TF_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(APP_DIR)/junctionFinder: $(JF_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(APP_DIR)/sequenceQuality: $(SQ_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(APP_DIR)/wcdInterrogate: $(WCD_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LDFLAGS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
