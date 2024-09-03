# Compiler and flags
MPICC = mpicc
# -DTIME is defined when the main loop needs to be timed.
# Comment out accordingly which ones don't want to be used
# and recompile the code.

CFLAGS = -O3 -Wall -std=c99 $(DEFINE)
LDFLAGS = -lm -lmpi

# Project structure
SRC = src
OBJ = obj
EXE = automaton
VPATH = $(SRC):$(addprefix $(SRC)/, mpi-lib heat-solver-lib util ser-lib par-lib wrap-lib)
INCLUDES = -Iinclude $(addprefix -I, $(subst :, ,$(VPATH)))

# Source files and objects
UTIL_SRCS = mem.c args.c arralloc.c misc.c
AUTOMATON_SRCS = heat_solver.c
MP_SRCS = mpilib.c
VER_SRCS = serlib.c parlib.c wraplib.c
MAIN_SRCS = main.c

UTIL_OBJS = $(UTIL_SRCS:%.c=$(OBJ)/%.o)
AUTOMATON_OBJS = $(AUTOMATON_SRCS:%.c=$(OBJ)/%.o)
MP_OBJS = $(MP_SRCS:%.c=$(OBJ)/%.o)
VER_OBJS = $(VER_SRCS:%.c=$(OBJ)/%.o)
MAIN_OBJS = $(MAIN_SRCS:%.c=$(OBJ)/%.o)

# Compilation rules
COMPILE = $(MPICC) $(CFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: all clean

all: $(OBJ) $(EXE)

$(OBJ):
	mkdir -p $@

$(EXE): $(UTIL_OBJS) $(AUTOMATON_OBJS) $(MP_OBJS) $(VER_OBJS) $(MAIN_OBJS)
	$(MPICC) $(LDFLAGS) -o $@ $^

$(OBJ)/%.o: %.c
	$(COMPILE)

clean:
	rm -rf $(EXE) $(OBJ) core