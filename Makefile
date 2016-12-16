CFLAGS =

PROGRAM=sn

# options
DEBUG = yes
WARNING = yes
OPTIMIZE = yes

ifeq ($(DEBUG), yes)
	CFLAGS += -g
endif

ifeq ($(WARNING), yes)
	CFLAGS += -Wall
endif

all: main.cpp Solver.h Quadrature.h
	g++ $(CFLAGS) -c main.cpp
	g++ main.o -o $(PROGRAM)

.PHONY: clean
clean:
	rm -rf $(PROGRAM) main.o

.PHONY: run
run:
	./$(PROGRAM)
