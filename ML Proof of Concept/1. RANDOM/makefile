TRUE=100
FALSE=0
SYSTEM = x86-64_linux
LIBFORMAT = static_pic

CPLEXDIR = /home/fakcay/CPLEX_Studio2211/cplex
CPDIR = /home/fakcay/CPLEX_Studio2211/cpoptimizer
CONCERTDIR = /home/fakcay/CPLEX_Studio2211/concert
CPLEXBINDIR = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLIBDIR = $(CPDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CPLEXFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -I$(CPDIR)/include  -L$(CPLIBDIR) -lcp -L$(CPLEXLIBDIR) -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -ldl 

GUROBI_DIR = /home/fakcay/opt/gurobi1103/linux64
CLIB     = -L$(GUROBI_DIR)/lib/ -lgurobi110
INC      = $(GUROBI_DIR)/include/
CPPLIB   = -L$(GUROBI_DIR)/lib/ -lgurobi_c++ $(CLIB)
GUROBI_OPTS = -I$(INC) $(CPPLIB) -lpthread -lm -m64

CC = g++ 
CFLAGS = -O2 -Wall -ansi -pedantic -DIL_STD -std=c++11
DEBUG = -pg -g -Wall -ansi -pedantic -DIL_STD 
OBJECTS = main.o helper_functions.o model.o
 
exec : $(OBJECTS) 
	$(CC) $(CFLAGS) -o R_${TRUE}T_${FALSE}F $(OBJECTS) $(CPLEXFLAGS) $(GUROBI_OPTS)
 
.cpp.o : 
	$(CC) $(CFLAGS) $(CPLEXFLAGS) $(GUROBI_OPTS) -c $< -o $@ 
 
clean : 
	rm -f $(OBJECTS) R_${TRUE}T_${FALSE}F