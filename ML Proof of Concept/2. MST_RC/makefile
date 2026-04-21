SYSTEM = x86-64_linux
LIBFORMAT = static_pic

GUROBI_DIR = /home/fakcay/opt/gurobi1103/linux64
CLIB     = -L$(GUROBI_DIR)/lib/ -lgurobi110
INC      = $(GUROBI_DIR)/include/
CPPLIB   = -L$(GUROBI_DIR)/lib/ -lgurobi_c++ $(CLIB)
GUROBI_OPTS = -I$(INC) $(CPPLIB) -lpthread -lm -m64

CC = g++
CFLAGS = -O2 -Wall -DIL_STD -std=c++17

SRC = main.cpp helper_functions.cpp cycleGen_helper_functions.cpp model.cpp

TRUE_VALUES  = 20 40 60 80 100
FALSE_VALUES = 50

TARGETS = $(foreach t,$(TRUE_VALUES),$(foreach f,$(FALSE_VALUES),T$(t)_F$(f)))

all: $(TARGETS)

# Build each target
$(TARGETS):
	$(eval PARTS := $(subst _F, ,$(patsubst T%,%,$@)))
	$(eval T := $(word 1,$(PARTS)))
	$(eval F := $(word 2,$(PARTS)))
	$(CC) $(CFLAGS) -DPRECUT_PERC=$(T) -DPRECUT_MULT=$(F) -o $@ $(SRC) $(GUROBI_OPTS)

clean:
	rm -f $(TARGETS) *.o