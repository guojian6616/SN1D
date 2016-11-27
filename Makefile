CFLAGS =

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
	g++ main.o -o tt

.PHONY: clean
clean:
	rm -rf tt main.o

.PHONY: run
run:
	./tt
