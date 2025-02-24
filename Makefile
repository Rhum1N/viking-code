#===============================================================
# Makefile for building C Code in command line environement
# using the GCC Open Source Compiler
# Edited by Enzo IGLESIS
#===============================================================

######################################

RM						:= mv
MKDIR					= mkdir -p -- $@

DONE						= [\e[32mDONE\e[0;m]:
NOTE					= [\e[36mNOTE\e[0;m]:

BOLD					=\e[1m
UNDERLINE				=\e[4m
DEFAULT					=\e[0;m

######################################

ifndef TARGET
    $(info Target unspecified, default: TARGET=main.out)
    TARGET 			:= main.out
endif

INC_DIR				:= libraries

SOURCES 			:= $(patsubst src/%.cpp, %.cpp, $(wildcard src/*.cpp))
SOURCES             += $(patsubst $(INC_DIR)/src/%.cpp, %.cpp, $(wildcard $(INC_DIR)/src/*.cpp))

OBJ_DIR 			:= output


######################################

CC              		:= g++
GDB			    		:= gdb
INCLUDES 				:= -I $(INC_DIR)
COMMONFLAGS 			:= -std=c++11 -Wall -fopenmp
CFLAGS 					:= $(COMMONFLAGS) 
LDFLAGS 				:= $(COMMONFLAGS) 

######################################

VPATH := src:$(INC_DIR)/src/

OBJECTS := $(addsuffix .o,$(addprefix $(OBJ_DIR)/, $(basename $(SOURCES))))

######################################

all: $(OBJ_DIR)/$(TARGET)

$(OBJECTS): | $(OBJ_DIR)

$(OBJ_DIR):
	@echo "created directory '$(BOLD)$(UNDERLINE)$(OBJ_DIR)/$(DEFAULT)'"
	@$(MKDIR) $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp
	@echo ======================================================================
	@echo -e "$(NOTE) Generating $(BOLD)$(UNDERLINE)$@$(DEFAULT)"
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)/$(TARGET): $(OBJECTS)
	@echo ======================================================================
	@echo -e "$(NOTE) Linking objects and generating output binary $(BOLD)$(UNDERLINE)$@$(DEFAULT)"
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(INCLUDES)
	@echo "$(DONE) output binary: $(BOLD)$(UNDERLINE)$@$(DEFAULT)"

debug: all
	$(GDB) $(OBJ_DIR)/$(TARGET)

run: all
	$(MKDIR) test
	$(OBJ_DIR)/$(TARGET) test


clean:
	@$(RM) $(OBJ_DIR) trash

rebuild: clean all
